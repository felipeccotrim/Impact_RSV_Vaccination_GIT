################################################################################
# REPRODUÇÃO DAS ANÁLISES DO ARTIGO
#
# "Potential impact of maternal vaccination against respiratory syncytial virus
#  on severe acute respiratory infections in infants up to six months of age in
#  Brazil, 2026"
#
# Arquivo único com as análises de ESTIMATIVA DE IMPACTO do artigo, na ordem:
#   0. configuração geral e funções compartilhadas;
#   1. preparação das séries mensais (casos de até 6 meses, por macrorregião);
#   2. validação metodológica: backtesting (2019, 2022-2025) e seleção do
#      algoritmo (XGBoost, LightGBM, CatBoost);
#   3. estimativa de impacto: contrafactual sem vacinação para 2026 (CatBoost),
#      intervalo de previsão empírico de 95%, cenários de vacinação
#      (efetividade x cobertura) e análise de sensibilidade com ajuste por viés;
#   4. comparação exploratória e não causal com os dados observados de 2026 e
#      com a cobertura vacinal real;
#   5. figuras do manuscrito.
#
# ESCOPO E PRÉ-REQUISITOS
#   Este arquivo concentra as análises de estimativa de impacto. Ele NÃO prepara
#   as bases brutas: pressupõe que os arquivos de BD/ foram gerados por scripts
#   de preparação anteriores (fora deste repositório), que aplicam os critérios
#   de elegibilidade aos dados públicos do OpenDATASUS, em especial a seleção de
#   casos de SRAG com RT-PCR positivo para VSR (RES_VSR == 1 no SINAN 2013-2018;
#   PCR_VSR == 1 no SIVEP-Gripe 2019-2025). Para 2026, o filtro PCR_VSR == "1" é
#   aplicado aqui (seção 4). A restrição a crianças de até 6 meses
#   (faixa_idade == "<=6 meses") é feita aqui, na seção 1.
#
# COMO EXECUTAR
#   1. Abra impacto_vacina_vsr.Rproj (os caminhos usam here::here()).
#   2. Coloque as bases em BD/ (ver README.md, seção "Dados"):
#        - BD/base_VSR_2013_2018.RData           (objeto base_VSR_2013_2018)
#        - BD/base_DEF_VSR_2019_2025.RData       (objeto base_DEF_VSR_2019_2025)
#        - BD/INFLUD26-27-07-2026.csv            (SIVEP-Gripe, extração de 27/07/2026)
#        - BD/cobertura_vacinal_vsr_gestantes_uf_2026.xlsx (painel de cobertura)
#   3. Execute este arquivo do início ao fim (source("reproducao_artigo_vsr.R")).
#      A execução completa levou cerca de 30 segundos em um computador pessoal.
#
# SAÍDAS (não versionadas)
#   resultados/bases_processadas/, resultados/tabelas/, resultados/graficos/ e
#   resultados/comparacao_modelos_preditivos/.
#
# SOFTWARE
#   R 4.5.3; catboost 1.2.10; xgboost 3.2.1.1; lightgbm 4.7.0. Semente: 123.
################################################################################


################################################################################
# 0. CONFIGURAÇÃO GERAL E FUNÇÕES COMPARTILHADAS
################################################################################

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

pacotes <- c(
  "dplyr", "tidyr", "lubridate", "slider", "purrr",
  "tibble", "ggplot2", "scales", "readr", "writexl",
  "xgboost", "lightgbm", "catboost", "zoo", "here",
  "stringr", "readxl"
)

instalar_ausentes <- function(pacotes) {
  ausentes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
  if (length(ausentes) > 0) {
    install.packages(ausentes)
  }
}

instalar_ausentes(pacotes)

invisible(lapply(
  pacotes,
  library,
  character.only = TRUE
))

set.seed(123)

# ---- Funções de modelagem (engenharia de atributos, ajuste e previsão recursiva)
################################################################################
# FUNÇÕES COMPARTILHADAS DE MODELAGEM — SRAG POR VSR
#
# Utilizado pelos Scripts 02 e 03 para garantir que a engenharia de atributos,
# o ajuste dos algoritmos, a previsão recursiva e as métricas sejam idênticos.
################################################################################

criar_atributos <- function(y, datas) {
  stopifnot(length(y) == length(datas))

  indice <- seq_along(y)
  mes_num <- lubridate::month(datas)
  ano_num <- lubridate::year(datas)

  media_movel_defasada <- function(x, janela) {
    dplyr::lag(
      zoo::rollapplyr(
        x,
        width = janela,
        FUN = mean,
        fill = NA_real_,
        partial = FALSE,
        na.rm = TRUE
      ),
      1L
    )
  }

  tibble::tibble(
    tendencia = indice,
    ano = ano_num,
    mes_num = mes_num,
    mes_sin = sin(2 * pi * mes_num / 12),
    mes_cos = cos(2 * pi * mes_num / 12),
    periodo_pre_pandemia = as.integer(datas < as.Date("2020-03-01")),
    periodo_pandemia_transicao = as.integer(
      datas >= as.Date("2020-03-01") & datas < as.Date("2023-01-01")
    ),
    periodo_recente = as.integer(datas >= as.Date("2023-01-01")),
    tempo_desde_2023 = pmax(
      0,
      (lubridate::year(datas) - 2023) * 12 + lubridate::month(datas) - 1
    ),
    lag1 = dplyr::lag(y, 1L),
    lag2 = dplyr::lag(y, 2L),
    lag3 = dplyr::lag(y, 3L),
    lag12 = dplyr::lag(y, 12L),
    media_movel3 = media_movel_defasada(y, 3L),
    media_movel6 = media_movel_defasada(y, 6L),
    diferenca_mensal = dplyr::lag(y, 1L) - dplyr::lag(y, 2L),
    razao_lag1_lag12 = dplyr::if_else(
      dplyr::lag(y, 12L) > 0,
      dplyr::lag(y, 1L) / dplyr::lag(y, 12L),
      NA_real_
    )
  )
}

NOMES_ATRIBUTOS <- c(
  "tendencia", "ano", "mes_num", "mes_sin", "mes_cos",
  "periodo_pre_pandemia", "periodo_pandemia_transicao",
  "periodo_recente", "tempo_desde_2023",
  "lag1", "lag2", "lag3", "lag12",
  "media_movel3", "media_movel6",
  "diferenca_mensal", "razao_lag1_lag12"
)

preparar_dados_treino <- function(df) {
  df <- df |>
    dplyr::arrange(.data$mes)

  dplyr::bind_cols(
    tibble::tibble(mes = df$mes, alvo = df$casos),
    criar_atributos(df$casos, df$mes)
  ) |>
    tidyr::drop_na()
}

criar_linha_futura <- function(historico_valores, historico_datas, nova_data) {
  criar_atributos(
    c(historico_valores, NA_real_),
    c(historico_datas, nova_data)
  ) |>
    dplyr::slice_tail(n = 1L)
}

ajustar_modelo <- function(dados, algoritmo, parametros) {
  x <- dados |>
    dplyr::select(dplyr::all_of(NOMES_ATRIBUTOS)) |>
    as.matrix()

  y <- dados$alvo

  switch(
    algoritmo,
    XGBOOST = {
      dtrain <- xgboost::xgb.DMatrix(data = x, label = y)
      xgboost::xgb.train(
        params = parametros$PARAM_XGBOOST,
        data = dtrain,
        nrounds = parametros$NROUND_XGBOOST,
        verbose = 0
      )
    },
    LIGHTGBM = {
      if (!requireNamespace("lightgbm", quietly = TRUE)) {
        stop("O pacote lightgbm não está instalado.")
      }
      dtrain <- lightgbm::lgb.Dataset(data = x, label = y)
      lightgbm::lgb.train(
        params = parametros$PARAM_LIGHTGBM,
        data = dtrain,
        nrounds = parametros$NROUND_LIGHTGBM,
        verbose = -1
      )
    },
    CATBOOST = {
      if (!requireNamespace("catboost", quietly = TRUE)) {
        stop("O pacote catboost não está instalado.")
      }
      pool <- catboost::catboost.load_pool(data = x, label = y)
      catboost::catboost.train(
        learn_pool = pool,
        params = parametros$PARAM_CATBOOST
      )
    },
    stop("Algoritmo desconhecido: ", algoritmo)
  )
}

prever_uma_linha <- function(modelo, nova_linha, algoritmo) {
  x <- nova_linha |>
    dplyr::select(dplyr::all_of(NOMES_ATRIBUTOS)) |>
    as.matrix()

  switch(
    algoritmo,
    XGBOOST = as.numeric(
      predict(modelo, xgboost::xgb.DMatrix(x))
    ),
    LIGHTGBM = as.numeric(predict(modelo, x)),
    CATBOOST = {
      pool <- catboost::catboost.load_pool(data = x)
      as.numeric(catboost::catboost.predict(modelo, pool))
    },
    stop("Algoritmo desconhecido: ", algoritmo)
  )
}

ajustar_prever_recursivo <- function(
  df_treino,
  algoritmo,
  datas_futuras,
  parametros
) {
  dados_modelo <- preparar_dados_treino(df_treino)

  if (nrow(dados_modelo) < 24L) {
    stop(
      "Número insuficiente de linhas após a criação dos atributos: ",
      nrow(dados_modelo)
    )
  }

  modelo <- ajustar_modelo(dados_modelo, algoritmo, parametros)
  historico_valores <- df_treino$casos
  historico_datas <- df_treino$mes
  previsoes <- numeric(length(datas_futuras))

  for (i in seq_along(datas_futuras)) {
    if (i == 1L) {
      valores_anteriores <- numeric()
      datas_anteriores <- as.Date(character())
    } else {
      valores_anteriores <- previsoes[seq_len(i - 1L)]
      datas_anteriores <- datas_futuras[seq_len(i - 1L)]
    }

    nova_linha <- criar_linha_futura(
      historico_valores = c(historico_valores, valores_anteriores),
      historico_datas = c(historico_datas, datas_anteriores),
      nova_data = datas_futuras[i]
    )

    previsoes[i] <- max(
      0,
      prever_uma_linha(modelo, nova_linha, algoritmo)
    )
  }

  list(
    modelo = modelo,
    previsoes = previsoes,
    dados_modelo = dados_modelo
  )
}

calcular_metricas <- function(observado, previsto, treino = NULL) {
  valido <- is.finite(observado) & is.finite(previsto)
  observado <- observado[valido]
  previsto <- previsto[valido]

  if (length(observado) == 0L) {
    return(tibble::tibble(
      n = 0L, MAE = NA_real_, RMSE = NA_real_,
      sMAPE = NA_real_, WAPE = NA_real_, MASE = NA_real_,
      R2 = NA_real_, correlacao = NA_real_,
      vies = NA_real_, razao_total = NA_real_
    ))
  }

  erros <- previsto - observado
  denominador_smape <- abs(observado) + abs(previsto)

  smape <- mean(
    ifelse(
      denominador_smape == 0,
      0,
      200 * abs(erros) / denominador_smape
    ),
    na.rm = TRUE
  )

  wape <- if (sum(abs(observado)) == 0) {
    NA_real_
  } else {
    100 * sum(abs(erros)) / sum(abs(observado))
  }

  mase <- NA_real_
  if (!is.null(treino)) {
    escala_mase <- mean(abs(diff(treino, lag = 12L)), na.rm = TRUE)
    if (is.finite(escala_mase) && escala_mase > 0) {
      mase <- mean(abs(erros), na.rm = TRUE) / escala_mase
    }
  }

  sse <- sum((observado - previsto)^2)
  sst <- sum((observado - mean(observado))^2)

  correlacao <- if (
    length(observado) >= 3L &&
    stats::sd(observado) > 0 &&
    stats::sd(previsto) > 0
  ) {
    stats::cor(observado, previsto)
  } else {
    NA_real_
  }

  tibble::tibble(
    n = length(observado),
    MAE = mean(abs(erros)),
    RMSE = sqrt(mean(erros^2)),
    sMAPE = smape,
    WAPE = wape,
    MASE = mase,
    R2 = if (sst == 0) NA_real_ else 1 - sse / sst,
    correlacao = correlacao,
    vies = mean(erros),
    razao_total = ifelse(
      sum(observado) == 0,
      NA_real_,
      sum(previsto) / sum(observado)
    )
  )
}


# ---- Funções de harmonização (SINAN / SIVEP-Gripe)
################################################################################
# FUNÇÕES COMPARTILHADAS DE HARMONIZAÇÃO — SRAG POR VSR
#
# Utilizado pelos Scripts 01 e 03 para garantir que qualquer extração de dados
# (SINAN Influenza Web, SIVEP-Gripe histórico ou SIVEP-Gripe 2026) seja
# processada exatamente da mesma forma: normalização de UF/região, conversão
# de datas e derivação de idade/faixa etária.
################################################################################

uf_codigo_sigla <- c(
  "11" = "RO", "12" = "AC", "13" = "AM", "14" = "RR",
  "15" = "PA", "16" = "AP", "17" = "TO",
  "21" = "MA", "22" = "PI", "23" = "CE", "24" = "RN",
  "25" = "PB", "26" = "PE", "27" = "AL", "28" = "SE",
  "29" = "BA", "31" = "MG", "32" = "ES", "33" = "RJ",
  "35" = "SP", "41" = "PR", "42" = "SC", "43" = "RS",
  "50" = "MS", "51" = "MT", "52" = "GO", "53" = "DF"
)

uf_regiao <- c(
  AC = "Norte", AP = "Norte", AM = "Norte", PA = "Norte",
  RO = "Norte", RR = "Norte", TO = "Norte",
  AL = "Nordeste", BA = "Nordeste", CE = "Nordeste",
  MA = "Nordeste", PB = "Nordeste", PE = "Nordeste",
  PI = "Nordeste", RN = "Nordeste", SE = "Nordeste",
  DF = "Centro-Oeste", GO = "Centro-Oeste",
  MT = "Centro-Oeste", MS = "Centro-Oeste",
  ES = "Sudeste", MG = "Sudeste", RJ = "Sudeste", SP = "Sudeste",
  PR = "Sul", RS = "Sul", SC = "Sul"
)

normalizar_uf <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- ifelse(z %in% names(uf_codigo_sigla), uf_codigo_sigla[z], z)
  z[!z %in% names(uf_regiao)] <- NA_character_
  unname(z)
}

converter_data <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }

  # Algumas extrações trazem a data com sufixo de horário/fuso
  # (ex.: "2026-01-06T00:00:00.000Z"). O prefixo AAAA-MM-DD é preservado.
  x_chr <- substr(as.character(x), 1, 10)

  suppressWarnings(
    dplyr::coalesce(
      lubridate::ymd(x_chr, quiet = TRUE),
      lubridate::dmy(x_chr, quiet = TRUE),
      as.Date(x_chr)
    )
  )
}

escolher_coluna <- function(df, candidatas, obrigatoria = FALSE) {
  achadas <- intersect(candidatas, names(df))

  if (length(achadas) == 0) {
    if (obrigatoria) {
      stop(
        "Nenhuma das colunas esperadas foi encontrada: ",
        paste(candidatas, collapse = ", ")
      )
    }
    return(NA_character_)
  }

  achadas[[1]]
}

harmonizar_base <- function(df) {
  col_uf <- escolher_coluna(
    df,
    c("SG_UF", "SG_UF_NOT", "SG_UF_INTE"),
    obrigatoria = TRUE
  )

  col_nasc <- escolher_coluna(df, c("DT_NASC"))
  col_tp_idade <- escolher_coluna(df, c("TP_IDADE"))
  col_nu_idade <- escolher_coluna(df, c("NU_IDADE_N"))

  saida <- df |>
    dplyr::mutate(
      DT_SIN_PRI = converter_data(.data[["DT_SIN_PRI"]]),
      SG_UF = normalizar_uf(.data[[col_uf]])
    )

  saida$DT_NASC_PAD <- if (!is.na(col_nasc)) {
    converter_data(saida[[col_nasc]])
  } else {
    as.Date(NA)
  }

  saida$TP_IDADE_PAD <- if (!is.na(col_tp_idade)) {
    suppressWarnings(as.integer(as.numeric(saida[[col_tp_idade]])))
  } else {
    NA_integer_
  }

  saida$NU_IDADE_PAD <- if (!is.na(col_nu_idade)) {
    suppressWarnings(as.numeric(saida[[col_nu_idade]]))
  } else {
    NA_real_
  }

  saida |>
    dplyr::mutate(
      ANO_SIN_PRI = lubridate::year(DT_SIN_PRI),
      REGIAO = unname(uf_regiao[SG_UF])
    )
}

derivar_idade_faixa <- function(df) {
  df |>
    dplyr::mutate(
      idade_dias_nascimento = as.numeric(DT_SIN_PRI - DT_NASC_PAD),

      idade_dias_ficha = dplyr::case_when(
        TP_IDADE_PAD == 1 ~ NU_IDADE_PAD,
        TP_IDADE_PAD == 2 ~ NU_IDADE_PAD * 30.4375,
        TP_IDADE_PAD == 3 ~ NU_IDADE_PAD * 365.25,
        TRUE ~ NA_real_
      ),

      idade_dias = dplyr::coalesce(
        dplyr::if_else(
          idade_dias_nascimento >= 0,
          idade_dias_nascimento,
          NA_real_
        ),
        idade_dias_ficha
      ),

      idade_meses = idade_dias / 30.4375,
      idade_anos = idade_dias / 365.25,

      faixa_idade = dplyr::case_when(
        idade_meses >= 0 & idade_meses <= 6 ~ "<=6 meses",
        idade_meses > 6 & idade_meses <= 12 ~ "7 a 12 meses",
        idade_meses > 12 & idade_meses < 24 ~ "13 a <24 meses",
        idade_anos >= 2 & idade_anos < 6 ~ "2 a 5 anos",
        idade_anos >= 6 & idade_anos < 16 ~ "6 a 15 anos",
        idade_anos >= 16 & idade_anos < 41 ~ "16 a 40 anos",
        idade_anos >= 41 & idade_anos < 60 ~ "41 a 59 anos",
        idade_anos >= 60 ~ "60 anos ou mais",
        TRUE ~ NA_character_
      )
    )
}


################################################################################
# 1. PREPARAÇÃO DAS SÉRIES MENSAIS (casos de até 6 meses, por macrorregião)
################################################################################

#### Parâmetros da análise ####

ARQ_2013_2018 <- here::here("BD", "base_VSR_2013_2018.RData")
ARQ_2019_2025 <- here::here("BD", "base_DEF_VSR_2019_2025.RData")

# Com a base de 2025 completa, o modelo utiliza os dados observados
# até dezembro de 2025 e projeta janeiro a dezembro de 2026.
DATA_CORTE_TREINO <- as.Date("2025-12-01")
DATA_INICIO_PREVISAO <- as.Date("2026-01-01")
DATA_FINAL_FORECAST <- as.Date("2026-12-01")

# Nível do intervalo de previsão empírico (baseado nos resíduos de
# backtesting do modelo final, em escala log1p, por região).
IC_PROB_INFERIOR <- 0.025
IC_PROB_SUPERIOR <- 0.975

FAIXA_ALVO <- "<=6 meses"

REGIOES_ORDEM <- c(
  "Brasil", "Norte", "Nordeste",
  "Centro-Oeste", "Sudeste", "Sul"
)

CENARIOS_VACINA <- tibble::tribble(
  ~cenario,       ~efetividade, ~cobertura,
  "Conservador",          0.60,       0.50,
  "Intermediário",        0.70,       0.70,
  "Otimista",             0.75,       0.85
) |>
  dplyr::mutate(
    reducao_total = efetividade * cobertura
  )

DIR_BASES <- here::here("resultados", "bases_processadas")
DIR_TABELAS <- here::here("resultados", "tabelas")
DIR_GRAFICOS <- here::here("resultados", "graficos")

dir.create(DIR_BASES, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABELAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_GRAFICOS, recursive = TRUE, showWarnings = FALSE)

#### Carregar e harmonizar bases ####

stopifnot(file.exists(ARQ_2013_2018))
stopifnot(file.exists(ARQ_2019_2025))

load(ARQ_2013_2018)
load(ARQ_2019_2025)

if (!exists("base_VSR_2013_2018")) {
  stop("Objeto 'base_VSR_2013_2018' não encontrado no primeiro arquivo.")
}

if (!exists("base_DEF_VSR_2019_2025")) {
  stop("Objeto 'base_DEF_VSR_2019_2025' não encontrado no segundo arquivo.")
}

base_2013_2018_h <- harmonizar_base(base_VSR_2013_2018)
base_2019_2025_h <- harmonizar_base(base_DEF_VSR_2019_2025)

base_vsr <- dplyr::bind_rows(
  base_2013_2018_h,
  base_2019_2025_h
) |>
  dplyr::filter(!is.na(DT_SIN_PRI))

#### Derivar idade ####

base_vsr <- derivar_idade_faixa(base_vsr)

qa_base <- tibble::tibble(
  indicador = c(
    "Registros harmonizados",
    "Data de início dos sintomas ausente",
    "UF ausente ou inválida",
    "Região ausente",
    "Faixa etária ausente",
    "Menores ou iguais a 6 meses"
  ),
  valor = c(
    nrow(base_vsr),
    sum(is.na(base_vsr$DT_SIN_PRI)),
    sum(is.na(base_vsr$SG_UF)),
    sum(is.na(base_vsr$REGIAO)),
    sum(is.na(base_vsr$faixa_idade)),
    sum(base_vsr$faixa_idade == FAIXA_ALVO, na.rm = TRUE)
  )
)

readr::write_csv(
  qa_base,
  file.path(DIR_TABELAS, "01_qa_base.csv")
)

saveRDS(
  base_vsr,
  file.path(DIR_BASES, "01_base_vsr_harmonizada.rds")
)

#### Preparar séries mensais ####

base_vsr <- readRDS(
  file.path(DIR_BASES, "01_base_vsr_harmonizada.rds")
)

serie_regiao <- base_vsr |>
  dplyr::filter(
    faixa_idade == FAIXA_ALVO,
    !is.na(REGIAO),
    !is.na(DT_SIN_PRI)
  ) |>
  dplyr::mutate(
    mes = lubridate::floor_date(DT_SIN_PRI, unit = "month")
  ) |>
  # Utiliza toda a série observada disponível até dezembro de 2025.
  # Registros de 2026 não entram no ajuste contrafactual.
  dplyr::filter(mes <= DATA_CORTE_TREINO) |>
  dplyr::count(REGIAO, mes, name = "casos") |>
  dplyr::group_by(REGIAO) |>
  tidyr::complete(
    mes = seq.Date(min(mes), DATA_CORTE_TREINO, by = "month"),
    fill = list(casos = 0)
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(regiao = REGIAO)

serie_brasil <- serie_regiao |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    casos = sum(casos),
    .groups = "drop"
  ) |>
  dplyr::mutate(regiao = "Brasil")

serie_geo <- dplyr::bind_rows(
  serie_brasil,
  serie_regiao
) |>
  dplyr::mutate(
    regiao = factor(regiao, levels = REGIOES_ORDEM)
  ) |>
  dplyr::arrange(regiao, mes)

readr::write_csv(
  serie_geo,
  file.path(DIR_BASES, "02_serie_mensal_menores_6m.csv")
)

saveRDS(
  serie_geo,
  file.path(DIR_BASES, "02_serie_mensal_menores_6m.rds")
)


# Série mensal completa (2013-01 a 2025-12, todas as regiões com zeros preenchidos
# desde 2013-01), usada como entrada da validação metodológica (seção 2).
serie_completa_2013_2025 <- serie_geo |>
  dplyr::mutate(regiao = as.character(regiao)) |>
  dplyr::group_by(regiao) |>
  tidyr::complete(
    mes = seq.Date(as.Date("2013-01-01"), DATA_CORTE_TREINO, by = "month"),
    fill = list(casos = 0L)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(regiao = factor(regiao, levels = REGIOES_ORDEM)) |>
  dplyr::arrange(regiao, mes)

saveRDS(
  serie_completa_2013_2025,
  file.path(DIR_BASES, "02_serie_mensal_completa_2013_2025.rds")
)

################################################################################
# 2. VALIDAÇÃO METODOLÓGICA: BACKTESTING (2019, 2022-2025) E SELEÇÃO DO ALGORITMO
################################################################################

local({
#### 01. Configuração ####

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(123)

PACOTES_CRAN <- c(
  "dplyr", "tidyr", "purrr", "tibble", "lubridate",
  "ggplot2", "scales", "readr", "writexl", "zoo",
  "here", "stringr", "xgboost"
)

instalar_ausentes <- function(pacotes) {
  ausentes <- pacotes[
    !vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(ausentes) > 0L) {
    install.packages(ausentes)
  }
}

instalar_ausentes(PACOTES_CRAN)
invisible(lapply(PACOTES_CRAN, library, character.only = TRUE))

# (as funções de modelagem já foram definidas na seção 0 deste arquivo)

ARQ_SERIE <- here::here(
  "resultados",
  "bases_processadas",
  "02_serie_mensal_completa_2013_2025.rds"
)

DIR_RESULTADOS <- here::here(
  "resultados",
  "comparacao_modelos_preditivos"
)

DIR_TABELAS <- file.path(DIR_RESULTADOS, "tabelas")
DIR_GRAFICOS <- file.path(DIR_RESULTADOS, "graficos")
DIR_OBJETOS <- file.path(DIR_RESULTADOS, "objetos")

purrr::walk(
  c(DIR_TABELAS, DIR_GRAFICOS, DIR_OBJETOS),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

DATA_INICIO <- as.Date("2013-01-01")
DATA_FIM_DESENVOLVIMENTO <- as.Date("2025-12-01")

ANOS_BACKTEST <- c(2019L, 2022L, 2023L, 2024L, 2025L)
ANOS_PANDEMICOS_EXCLUIDOS_AVALIACAO <- c(2020L, 2021L)

REGIOES_ORDEM <- c(
  "Brasil", "Norte", "Nordeste",
  "Centro-Oeste", "Sudeste", "Sul"
)

ALGORITMO_CALIBRACAO <- "XGBOOST"

# Regra operacional para sinalizar viés sistemático.
# O viés será considerado sistemático quando:
#   a) pelo menos 75% dos anos tiverem a mesma direção de erro; e
#   b) a razão média previsto/observado estiver fora de 0,90–1,10.
TOLERANCIA_RAZAO_INFERIOR <- 0.90
TOLERANCIA_RAZAO_SUPERIOR <- 1.10
PROPORCAO_MINIMA_MESMA_DIRECAO <- 0.75

MODELOS_DISPONIVEIS <- c(
  XGBOOST = requireNamespace("xgboost", quietly = TRUE),
  LIGHTGBM = requireNamespace("lightgbm", quietly = TRUE),
  CATBOOST = requireNamespace("catboost", quietly = TRUE)
)

MODELOS <- names(MODELOS_DISPONIVEIS)[MODELOS_DISPONIVEIS]
MODELOS_AUSENTES <- names(MODELOS_DISPONIVEIS)[!MODELOS_DISPONIVEIS]

if (length(MODELOS) == 0L) {
  stop("Nenhum algoritmo solicitado está disponível.")
}

if (length(MODELOS_AUSENTES) > 0L) {
  warning(
    "Algoritmos não disponíveis nesta execução: ",
    paste(MODELOS_AUSENTES, collapse = ", "),
    call. = FALSE
  )
}

PARAMETROS <- list(
  PARAM_XGBOOST = list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = 0.03,
    max_depth = 3L,
    min_child_weight = 3,
    subsample = 0.85,
    colsample_bytree = 0.85,
    lambda = 2,
    alpha = 0.1,
    nthread = 1L
  ),
  NROUND_XGBOOST = 600L,
  
  PARAM_LIGHTGBM = list(
    objective = "regression",
    metric = "rmse",
    learning_rate = 0.03,
    num_leaves = 7L,
    max_depth = 3L,
    min_data_in_leaf = 8L,
    feature_fraction = 0.85,
    bagging_fraction = 0.85,
    bagging_freq = 1L,
    lambda_l1 = 0.1,
    lambda_l2 = 2,
    verbosity = -1L,
    num_threads = 1L,
    seed = 123L
  ),
  NROUND_LIGHTGBM = 600L,
  
  PARAM_CATBOOST = list(
    loss_function = "RMSE",
    eval_metric = "RMSE",
    iterations = 600L,
    learning_rate = 0.03,
    depth = 3L,
    l2_leaf_reg = 5,
    random_seed = 123L,
    logging_level = "Silent",
    allow_writing_files = FALSE,
    thread_count = 1L
  )
)

#### 02. Validação das entradas ####

arquivos_necessarios <- c(
  ARQ_SERIE
)

arquivos_ausentes <- arquivos_necessarios[
  !file.exists(arquivos_necessarios)
]

if (length(arquivos_ausentes) > 0L) {
  stop(
    "Arquivos não encontrados:\n",
    paste0(" - ", arquivos_ausentes, collapse = "\n")
  )
}

# (funções de modelagem: seção 0)

#### 03. Leitura e preparação da série ####

serie <- readRDS(ARQ_SERIE) |>
  dplyr::transmute(
    regiao = as.character(.data$regiao),
    mes = as.Date(.data$mes),
    casos = pmax(0, as.numeric(.data$casos))
  ) |>
  dplyr::filter(
    .data$mes >= DATA_INICIO,
    .data$mes <= DATA_FIM_DESENVOLVIMENTO,
    .data$regiao %in% REGIOES_ORDEM
  ) |>
  dplyr::group_by(.data$regiao) |>
  tidyr::complete(
    mes = seq.Date(
      DATA_INICIO,
      DATA_FIM_DESENVOLVIMENTO,
      by = "month"
    ),
    fill = list(casos = 0)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(.data$regiao, .data$mes)

regioes_disponiveis <- intersect(
  REGIOES_ORDEM,
  unique(serie$regiao)
)

#### 04. Backtesting temporal dos algoritmos ####

combinacoes <- tidyr::crossing(
  regiao = regioes_disponiveis,
  ano_validacao = ANOS_BACKTEST,
  algoritmo = MODELOS
)

executar_backtest <- function(regiao, ano_validacao, algoritmo) {
  
  regiao_atual <- as.character(regiao)
  ano_atual <- as.integer(ano_validacao)
  algoritmo_atual <- as.character(algoritmo)
  
  inicio_teste <- as.Date(paste0(ano_atual, "-01-01"))
  fim_teste <- as.Date(paste0(ano_atual, "-12-01"))
  
  serie_regiao <- serie |>
    dplyr::filter(.data$regiao == .env$regiao_atual) |>
    dplyr::arrange(.data$mes)
  
  treino <- serie_regiao |>
    dplyr::filter(.data$mes < .env$inicio_teste)
  
  teste <- serie_regiao |>
    dplyr::filter(
      .data$mes >= .env$inicio_teste,
      .data$mes <= .env$fim_teste
    )
  
  if (nrow(treino) == 0L) {
    return(
      tibble::tibble(
        regiao = regiao_atual,
        ano_validacao = ano_atual,
        algoritmo = algoritmo_atual,
        mes = as.Date(character()),
        observado = numeric(),
        previsto = numeric(),
        status = character(),
        mensagem = character()
      ) |>
        dplyr::add_row(
          regiao = regiao_atual,
          ano_validacao = ano_atual,
          algoritmo = algoritmo_atual,
          mes = as.Date(NA),
          observado = NA_real_,
          previsto = NA_real_,
          status = "erro",
          mensagem = "Conjunto de treinamento vazio."
        )
    )
  }
  
  if (nrow(teste) != 12L) {
    return(
      tibble::tibble(
        regiao = regiao_atual,
        ano_validacao = ano_atual,
        algoritmo = algoritmo_atual,
        mes = if (nrow(teste) > 0L) teste$mes else as.Date(NA),
        observado = if (nrow(teste) > 0L) teste$casos else NA_real_,
        previsto = NA_real_,
        status = "erro",
        mensagem = paste(
          "Ano de validação não contém 12 meses:",
          nrow(teste)
        )
      )
    )
  }
  
  tryCatch(
    {
      ajuste <- ajustar_prever_recursivo(
        df_treino = treino,
        algoritmo = algoritmo_atual,
        datas_futuras = teste$mes,
        parametros = PARAMETROS
      )
      
      previsoes <- as.numeric(ajuste$previsoes)
      
      if (length(previsoes) != nrow(teste)) {
        stop(
          "A função ajustar_prever_recursivo retornou ",
          length(previsoes),
          " previsões para ",
          nrow(teste),
          " meses."
        )
      }
      
      if (any(!is.finite(previsoes))) {
        stop("A função ajustar_prever_recursivo retornou valores não finitos.")
      }
      
      tibble::tibble(
        regiao = regiao_atual,
        ano_validacao = ano_atual,
        algoritmo = algoritmo_atual,
        mes = teste$mes,
        observado = teste$casos,
        previsto = pmax(0, previsoes),
        status = "ok",
        mensagem = NA_character_
      )
    },
    error = function(e) {
      tibble::tibble(
        regiao = regiao_atual,
        ano_validacao = ano_atual,
        algoritmo = algoritmo_atual,
        mes = teste$mes,
        observado = teste$casos,
        previsto = NA_real_,
        status = "erro",
        mensagem = conditionMessage(e)
      )
    }
  )
}

message("Executando backtesting temporal...")

previsoes_backtesting <- purrr::pmap_dfr(
  combinacoes,
  \(regiao, ano_validacao, algoritmo) {
    executar_backtest(
      regiao = regiao,
      ano_validacao = ano_validacao,
      algoritmo = algoritmo
    )
  }
)

erros_execucao <- previsoes_backtesting |>
  dplyr::filter(.data$status == "erro") |>
  dplyr::distinct(
    .data$regiao,
    .data$ano_validacao,
    .data$algoritmo,
    .data$mensagem
  )

if (nrow(erros_execucao) > 0L) {
  warning(
    "Ocorreram erros em algumas combinações. Consulte a tabela de erros.",
    call. = FALSE
  )
}

n_combinacoes_ok <- previsoes_backtesting |>
  dplyr::filter(.data$status == "ok") |>
  dplyr::distinct(
    .data$regiao,
    .data$ano_validacao,
    .data$algoritmo
  ) |>
  nrow()

if (n_combinacoes_ok == 0L) {
  resumo_erros <- erros_execucao |>
    dplyr::count(.data$mensagem, sort = TRUE)
  
  print(resumo_erros, n = Inf, width = Inf)
  
  stop(
    paste0(
      "Nenhuma combinação do backtesting foi executada com sucesso. ",
      "A execução foi interrompida antes das métricas, calibração e gráficos. ",
      "Consulte a tabela impressa acima."
    ),
    call. = FALSE
  )
}

#### 05. Pergunta 1 — Qual algoritmo prevê melhor? ####

metricas_por_ano <- previsoes_backtesting |>
  dplyr::filter(.data$status == "ok") |>
  dplyr::distinct(
    .data$regiao,
    .data$ano_validacao,
    .data$algoritmo
  ) |>
  purrr::pmap_dfr(function(regiao, ano_validacao, algoritmo) {
    
    ano_atual <- as.integer(ano_validacao)
    regiao_atual <- as.character(regiao)
    algoritmo_atual <- as.character(algoritmo)
    
    previsoes_grupo <- previsoes_backtesting |>
      dplyr::filter(
        .data$status == "ok",
        .data$regiao == .env$regiao_atual,
        .data$ano_validacao == .env$ano_atual,
        .data$algoritmo == .env$algoritmo_atual
      )
    
    treino <- serie |>
      dplyr::filter(
        .data$regiao == .env$regiao_atual,
        lubridate::year(.data$mes) < .env$ano_atual
      ) |>
      dplyr::arrange(.data$mes) |>
      dplyr::pull(.data$casos)
    
    calcular_metricas(
      observado = previsoes_grupo$observado,
      previsto = previsoes_grupo$previsto,
      treino = treino
    ) |>
      dplyr::mutate(
        regiao = regiao_atual,
        ano_validacao = ano_atual,
        algoritmo = algoritmo_atual,
        .before = 1
      )
  })

grade_metricas_agregadas <- previsoes_backtesting |>
  dplyr::filter(.data$status == "ok") |>
  dplyr::distinct(.data$regiao, .data$algoritmo)

metricas_agregadas <- purrr::pmap_dfr(
  grade_metricas_agregadas,
  \(regiao, algoritmo) {
    
    regiao_atual <- as.character(regiao)
    algoritmo_atual <- as.character(algoritmo)
    
    dados_grupo <- previsoes_backtesting |>
      dplyr::filter(
        .data$status == "ok",
        .data$regiao == .env$regiao_atual,
        .data$algoritmo == .env$algoritmo_atual
      )
    
    calcular_metricas(
      observado = dados_grupo$observado,
      previsto = dados_grupo$previsto
    ) |>
      dplyr::mutate(
        regiao = regiao_atual,
        algoritmo = algoritmo_atual,
        .before = 1
      )
  }
)

ranking_modelos <- metricas_agregadas |>
  dplyr::group_by(.data$regiao) |>
  dplyr::mutate(
    rank_RMSE = rank(.data$RMSE, ties.method = "average"),
    rank_MAE = rank(.data$MAE, ties.method = "average"),
    rank_sMAPE = rank(.data$sMAPE, ties.method = "average"),
    rank_WAPE = rank(.data$WAPE, ties.method = "average"),
    rank_MASE = rank(.data$MASE, ties.method = "average"),
    desvio_razao = abs(.data$razao_total - 1),
    rank_razao = rank(
      .data$desvio_razao,
      ties.method = "average"
    ),
    posicao_media = rowMeans(
      cbind(
        .data$rank_RMSE,
        .data$rank_MAE,
        .data$rank_sMAPE,
        .data$rank_WAPE,
        .data$rank_MASE,
        .data$rank_razao
      ),
      na.rm = TRUE
    )
  ) |>
  dplyr::arrange(
    .data$regiao,
    .data$posicao_media,
    .data$RMSE
  ) |>
  dplyr::mutate(
    posicao = dplyr::row_number(),
    selecionado = .data$posicao == 1L
  ) |>
  dplyr::ungroup()

modelos_selecionados <- ranking_modelos |>
  dplyr::filter(.data$selecionado) |>
  dplyr::select(
    .data$regiao,
    algoritmo_selecionado = .data$algoritmo,
    .data$RMSE,
    .data$MAE,
    .data$sMAPE,
    .data$WAPE,
    .data$MASE,
    .data$R2,
    .data$correlacao,
    .data$vies,
    .data$razao_total,
    .data$posicao_media
  )



#### 06. Pergunta 2 — Qual modelo será adotado? ####

# Para manter um único pipeline nacional, o modelo final é definido pelo
# algoritmo mais vezes selecionado (posicao_media == 1, ranking composto de
# RMSE/MAE/sMAPE/WAPE/MASE/razão) entre as macrorregiões — em vez de um valor
# fixo no código, para que a decisão reflita sempre o resultado real do
# backtesting quando o script for reexecutado.

contagem_selecionados <- modelos_selecionados |>
  dplyr::count(.data$algoritmo_selecionado, sort = TRUE, name = "regioes_vencidas")

MODELO_FINAL <- contagem_selecionados$algoritmo_selecionado[[1]]

decisao_modelo_final <- tibble::tibble(
  modelo_final = MODELO_FINAL,
  criterio_decisao = paste0(
    "Algoritmo com melhor ranking composto (RMSE, MAE, sMAPE, WAPE, MASE, ",
    "razão previsto/observado) em ", contagem_selecionados$regioes_vencidas[[1]],
    " das ", nrow(modelos_selecionados), " macrorregiões avaliadas."
  ),
  anos_backtesting = paste(ANOS_BACKTEST, collapse = ", "),
  calibracao_adicional = "Não adotada",
  justificativa_calibracao = paste(
    "As estratégias de calibração histórica avaliadas anteriormente",
    "aumentaram os erros preditivos em relação à previsão bruta."
  )
)

#### 07. Pergunta 3 — O XGBoost apresenta viés histórico relevante? ####

diagnostico_vies_anual_xgboost <- previsoes_backtesting |>
  dplyr::filter(
    .data$status == "ok",
    toupper(.data$algoritmo) == toupper(.env$MODELO_FINAL)
  ) |>
  dplyr::group_by(
    .data$regiao,
    .data$ano_validacao,
    .data$algoritmo
  ) |>
  dplyr::summarise(
    total_observado = sum(.data$observado, na.rm = TRUE),
    total_previsto = sum(.data$previsto, na.rm = TRUE),
    diferenca_casos = .data$total_previsto - .data$total_observado,
    diferenca_absoluta_casos = abs(.data$diferenca_casos),
    vies_medio_mensal = mean(
      .data$previsto - .data$observado,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    razao_previsto_observado = dplyr::if_else(
      .data$total_observado > 0,
      .data$total_previsto / .data$total_observado,
      NA_real_
    ),
    erro_percentual_total = dplyr::if_else(
      .data$total_observado > 0,
      100 * .data$diferenca_casos / .data$total_observado,
      NA_real_
    ),
    direcao_vies = dplyr::case_when(
      .data$razao_previsto_observado >
        TOLERANCIA_RAZAO_SUPERIOR ~ "Superestimação",
      .data$razao_previsto_observado <
        TOLERANCIA_RAZAO_INFERIOR ~ "Subestimação",
      TRUE ~ "Dentro da tolerância"
    )
  ) |>
  dplyr::arrange(.data$regiao, .data$ano_validacao)

if (nrow(diagnostico_vies_anual_xgboost) == 0L) {
  stop(
    "Não há previsões válidas do XGBoost para avaliar o viés histórico.",
    call. = FALSE
  )
}

resumo_agregado_xgboost <- diagnostico_vies_anual_xgboost |>
  dplyr::group_by(.data$regiao, .data$algoritmo) |>
  dplyr::summarise(
    anos_avaliados = dplyr::n(),
    anos_usados = paste(
      sort(unique(.data$ano_validacao)),
      collapse = ", "
    ),
    total_observado = sum(.data$total_observado, na.rm = TRUE),
    total_previsto = sum(.data$total_previsto, na.rm = TRUE),
    diferenca_casos = .data$total_previsto - .data$total_observado,
    diferenca_absoluta_casos = abs(.data$diferenca_casos),
    razao_agregada_previsto_observado = dplyr::if_else(
      .data$total_observado > 0,
      .data$total_previsto / .data$total_observado,
      NA_real_
    ),
    erro_percentual_agregado = dplyr::if_else(
      .data$total_observado > 0,
      100 * .data$diferenca_casos / .data$total_observado,
      NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    direcao_agregada = dplyr::case_when(
      .data$razao_agregada_previsto_observado >
        TOLERANCIA_RAZAO_SUPERIOR ~ "Superestimação",
      .data$razao_agregada_previsto_observado <
        TOLERANCIA_RAZAO_INFERIOR ~ "Subestimação",
      TRUE ~ "Dentro da tolerância"
    )
  )

resumo_vies_xgboost <- diagnostico_vies_anual_xgboost |>
  dplyr::group_by(.data$regiao, .data$algoritmo) |>
  dplyr::summarise(
    anos_avaliados = dplyr::n(),
    anos_superestimacao = sum(
      .data$direcao_vies == "Superestimação",
      na.rm = TRUE
    ),
    anos_subestimacao = sum(
      .data$direcao_vies == "Subestimação",
      na.rm = TRUE
    ),
    anos_tolerancia = sum(
      .data$direcao_vies == "Dentro da tolerância",
      na.rm = TRUE
    ),
    proporcao_maior_direcao = max(
      .data$anos_superestimacao,
      .data$anos_subestimacao
    ) / .data$anos_avaliados,
    direcao_predominante = dplyr::case_when(
      .data$anos_superestimacao >
        .data$anos_subestimacao ~ "Superestimação",
      .data$anos_subestimacao >
        .data$anos_superestimacao ~ "Subestimação",
      TRUE ~ "Sem direção predominante"
    ),
    media_razao_anual = mean(
      .data$razao_previsto_observado,
      na.rm = TRUE
    ),
    mediana_razao_anual = stats::median(
      .data$razao_previsto_observado,
      na.rm = TRUE
    ),
    dp_razao_anual = stats::sd(
      .data$razao_previsto_observado,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    resumo_agregado_xgboost |>
      dplyr::select(
        .data$regiao,
        .data$razao_agregada_previsto_observado,
        .data$erro_percentual_agregado,
        .data$direcao_agregada
      ),
    by = "regiao"
  ) |>
  dplyr::mutate(
    razao_agregada_fora_tolerancia =
      .data$razao_agregada_previsto_observado <
      TOLERANCIA_RAZAO_INFERIOR |
      .data$razao_agregada_previsto_observado >
      TOLERANCIA_RAZAO_SUPERIOR,
    vies_sistematico =
      .data$proporcao_maior_direcao >=
      PROPORCAO_MINIMA_MESMA_DIRECAO &
      .data$razao_agregada_fora_tolerancia,
    conclusao_vies = dplyr::case_when(
      .data$vies_sistematico &
        .data$direcao_predominante == "Superestimação" ~
        "Evidência de viés sistemático de superestimação",
      .data$vies_sistematico &
        .data$direcao_predominante == "Subestimação" ~
        "Evidência de viés sistemático de subestimação",
      TRUE ~
        "Sem evidência operacional de viés sistemático"
    )
  )

#### 08. Pergunta 4 — Qual é a distribuição mensal dos resíduos? ####

# Convenção:
# residuo = observado - previsto
# residuo > 0: o modelo subestimou.
# residuo < 0: o modelo superestimou.

residuos_mensais_xgboost <- previsoes_backtesting |>
  dplyr::filter(
    .data$status == "ok",
    toupper(.data$algoritmo) == toupper(.env$MODELO_FINAL)
  ) |>
  dplyr::transmute(
    regiao = .data$regiao,
    ano_validacao = as.integer(.data$ano_validacao),
    mes = as.Date(.data$mes),
    mes_numero = lubridate::month(.data$mes),
    mes_nome = factor(
      lubridate::month(
        .data$mes,
        label = TRUE,
        abbr = TRUE,
        locale = "pt_BR"
      ),
      ordered = TRUE
    ),
    observado = as.numeric(.data$observado),
    previsto = as.numeric(.data$previsto),
    residuo = .data$observado - .data$previsto,
    erro_modelo = .data$previsto - .data$observado,
    erro_absoluto = abs(.data$previsto - .data$observado),
    erro_quadratico = (.data$previsto - .data$observado)^2,
    erro_percentual = dplyr::if_else(
      .data$observado > 0,
      100 * (.data$previsto - .data$observado) /
        .data$observado,
      NA_real_
    ),
    razao_previsto_observado = dplyr::if_else(
      .data$observado > 0,
      .data$previsto / .data$observado,
      NA_real_
    ),
    razao_observado_previsto = dplyr::if_else(
      .data$previsto > 0,
      .data$observado / .data$previsto,
      NA_real_
    ),
    residuo_log1p =
      log1p(.data$observado) - log1p(.data$previsto)
  ) |>
  dplyr::arrange(.data$regiao, .data$mes)

resumo_residuos_mensais <- residuos_mensais_xgboost |>
  dplyr::group_by(.data$regiao) |>
  dplyr::summarise(
    n_meses = dplyr::n(),
    media_residuo = mean(.data$residuo, na.rm = TRUE),
    mediana_residuo = stats::median(
      .data$residuo,
      na.rm = TRUE
    ),
    desvio_padrao_residuo = stats::sd(
      .data$residuo,
      na.rm = TRUE
    ),
    p025_residuo = stats::quantile(
      .data$residuo,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    p05_residuo = stats::quantile(
      .data$residuo,
      probs = 0.05,
      na.rm = TRUE,
      names = FALSE
    ),
    p10_residuo = stats::quantile(
      .data$residuo,
      probs = 0.10,
      na.rm = TRUE,
      names = FALSE
    ),
    p90_residuo = stats::quantile(
      .data$residuo,
      probs = 0.90,
      na.rm = TRUE,
      names = FALSE
    ),
    p95_residuo = stats::quantile(
      .data$residuo,
      probs = 0.95,
      na.rm = TRUE,
      names = FALSE
    ),
    p975_residuo = stats::quantile(
      .data$residuo,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    MAE_mensal = mean(.data$erro_absoluto, na.rm = TRUE),
    RMSE_mensal = sqrt(
      mean(.data$erro_quadratico, na.rm = TRUE)
    ),
    media_residuo_log1p = mean(
      .data$residuo_log1p,
      na.rm = TRUE
    ),
    dp_residuo_log1p = stats::sd(
      .data$residuo_log1p,
      na.rm = TRUE
    ),
    p025_residuo_log1p = stats::quantile(
      .data$residuo_log1p,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    p975_residuo_log1p = stats::quantile(
      .data$residuo_log1p,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

resumo_residuos_por_mes_calendario <- residuos_mensais_xgboost |>
  dplyr::group_by(
    .data$regiao,
    .data$mes_numero,
    .data$mes_nome
  ) |>
  dplyr::summarise(
    n_observacoes = dplyr::n(),
    media_residuo = mean(.data$residuo, na.rm = TRUE),
    mediana_residuo = stats::median(
      .data$residuo,
      na.rm = TRUE
    ),
    desvio_padrao_residuo = stats::sd(
      .data$residuo,
      na.rm = TRUE
    ),
    p025_residuo = stats::quantile(
      .data$residuo,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    p975_residuo = stats::quantile(
      .data$residuo,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$regiao, .data$mes_numero)

#### 09. Síntese metodológica ####

sintese_metodologica <- modelos_selecionados |>
  dplyr::full_join(
    resumo_vies_xgboost |>
      dplyr::select(
        .data$regiao,
        .data$conclusao_vies,
        .data$direcao_predominante,
        .data$razao_agregada_previsto_observado,
        .data$erro_percentual_agregado
      ),
    by = "regiao"
  ) |>
  dplyr::mutate(
    modelo_final_estudo = MODELO_FINAL,
    tipo_previsao_final = "Previsão bruta",
    calibracao_adotada = FALSE,
    justificativa_final = paste(
      "O XGBoost foi adotado como modelo único.",
      "A calibração histórica não foi incorporada por apresentar",
      "desempenho inferior à previsão bruta."
    )
  ) |>
  dplyr::rename(
    melhor_algoritmo_por_regiao =
      .data$algoritmo_selecionado
  )

#### 10. Exportações ####

configuracao <- tibble::tibble(
  item = c(
    "inicio_serie",
    "fim_desenvolvimento",
    "anos_backtesting",
    "anos_pandemicos_excluidos_avaliacao",
    "algoritmos_avaliados",
    "algoritmos_ausentes",
    "modelo_final",
    "tipo_previsao_final",
    "calibracao_adotada",
    "base_incerteza_futura",
    "regra_vies_tolerancia",
    "regra_vies_proporcao_direcao"
  ),
  valor = c(
    as.character(DATA_INICIO),
    as.character(DATA_FIM_DESENVOLVIMENTO),
    paste(ANOS_BACKTEST, collapse = ", "),
    paste(
      ANOS_PANDEMICOS_EXCLUIDOS_AVALIACAO,
      collapse = ", "
    ),
    paste(MODELOS, collapse = ", "),
    paste(MODELOS_AUSENTES, collapse = ", "),
    MODELO_FINAL,
    "Previsão bruta",
    "Não",
    "Distribuição empírica dos resíduos mensais históricos",
    paste0(
      TOLERANCIA_RAZAO_INFERIOR,
      " a ",
      TOLERANCIA_RAZAO_SUPERIOR
    ),
    as.character(PROPORCAO_MINIMA_MESMA_DIRECAO)
  )
)

tabelas_csv <- list(
  "01_previsoes_backtesting_algoritmos.csv" =
    previsoes_backtesting,
  "02_metricas_por_ano_algoritmo.csv" =
    metricas_por_ano,
  "03_metricas_agregadas_algoritmos.csv" =
    metricas_agregadas,
  "04_ranking_algoritmos_por_regiao.csv" =
    ranking_modelos,
  "05_modelos_selecionados_por_regiao.csv" =
    modelos_selecionados,
  "06_decisao_modelo_final.csv" =
    decisao_modelo_final,
  "07_diagnostico_vies_anual_xgboost.csv" =
    diagnostico_vies_anual_xgboost,
  "08_resumo_agregado_xgboost.csv" =
    resumo_agregado_xgboost,
  "09_resumo_vies_xgboost.csv" =
    resumo_vies_xgboost,
  "10_residuos_mensais_xgboost.csv" =
    residuos_mensais_xgboost,
  "11_resumo_residuos_mensais.csv" =
    resumo_residuos_mensais,
  "12_resumo_residuos_por_mes_calendario.csv" =
    resumo_residuos_por_mes_calendario,
  "13_sintese_metodologica.csv" =
    sintese_metodologica
)

purrr::iwalk(
  tabelas_csv,
  ~ readr::write_csv(
    .x,
    file.path(DIR_TABELAS, .y)
  )
)

if (nrow(erros_execucao) > 0L) {
  readr::write_csv(
    erros_execucao,
    file.path(
      DIR_TABELAS,
      "99_erros_execucao.csv"
    )
  )
}

atributos_objeto <- if (exists("NOMES_ATRIBUTOS")) {
  NOMES_ATRIBUTOS
} else {
  NULL
}

objeto_script02 <- list(
  configuracao = configuracao,
  parametros = PARAMETROS,
  atributos = atributos_objeto,
  anos_backtest = ANOS_BACKTEST,
  modelo_final = MODELO_FINAL,
  decisao_modelo_final = decisao_modelo_final,
  previsoes_backtesting = previsoes_backtesting,
  metricas_por_ano = metricas_por_ano,
  metricas_agregadas = metricas_agregadas,
  ranking_modelos = ranking_modelos,
  modelos_selecionados = modelos_selecionados,
  diagnostico_vies_anual_xgboost =
    diagnostico_vies_anual_xgboost,
  resumo_agregado_xgboost =
    resumo_agregado_xgboost,
  resumo_vies_xgboost =
    resumo_vies_xgboost,
  residuos_mensais_xgboost =
    residuos_mensais_xgboost,
  resumo_residuos_mensais =
    resumo_residuos_mensais,
  resumo_residuos_por_mes_calendario =
    resumo_residuos_por_mes_calendario,
  sintese_metodologica =
    sintese_metodologica
)

saveRDS(
  objeto_script02,
  file.path(
    DIR_OBJETOS,
    "resultados_selecao_validacao_modelo_script02.rds"
  )
)

writexl::write_xlsx(
  list(
    configuracao = configuracao,
    decisao_modelo = decisao_modelo_final,
    sintese = sintese_metodologica,
    ranking_algoritmos = ranking_modelos,
    modelos_selecionados = modelos_selecionados,
    metricas_por_ano = metricas_por_ano,
    metricas_agregadas = metricas_agregadas,
    vies_anual_xgboost =
      diagnostico_vies_anual_xgboost,
    resumo_agregado_xgboost =
      resumo_agregado_xgboost,
    resumo_vies_xgboost =
      resumo_vies_xgboost,
    residuos_mensais =
      residuos_mensais_xgboost,
    resumo_residuos =
      resumo_residuos_mensais,
    residuos_por_mes =
      resumo_residuos_por_mes_calendario
  ),
  file.path(
    DIR_TABELAS,
    "selecao_validacao_modelo_srag_vsr_V6.xlsx"
  )
)

#### 11. Gráficos ####

if (nrow(ranking_modelos) > 0L) {
  
  grafico_ranking <- ranking_modelos |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = reorder(
          .data$algoritmo,
          .data$posicao_media
        ),
        y = .data$posicao_media
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(
      ~regiao,
      scales = "free_y"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Ranking dos algoritmos por macrorregião",
      subtitle = paste(
        "Backtesting nos anos",
        paste(ANOS_BACKTEST, collapse = ", ")
      ),
      x = "Algoritmo",
      y = "Posição média"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  ggplot2::ggsave(
    file.path(
      DIR_GRAFICOS,
      "01_ranking_algoritmos.png"
    ),
    grafico_ranking,
    width = 14,
    height = 9,
    dpi = 300
  )
}

if (nrow(diagnostico_vies_anual_xgboost) > 0L) {
  
  grafico_vies <- diagnostico_vies_anual_xgboost |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = factor(.data$ano_validacao),
        y = .data$razao_previsto_observado,
        group = 1
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(
      ~regiao,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = "Viés anual do XGBoost",
      subtitle = "Razão entre o total previsto e o total observado",
      x = "Ano de validação",
      y = "Previsto/observado"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  ggplot2::ggsave(
    file.path(
      DIR_GRAFICOS,
      "02_vies_anual_xgboost.png"
    ),
    grafico_vies,
    width = 14,
    height = 9,
    dpi = 300
  )
}

dados_series_xgboost <- previsoes_backtesting |>
  dplyr::filter(
    .data$status == "ok",
    toupper(.data$algoritmo) == toupper(.env$MODELO_FINAL)
  )

if (nrow(dados_series_xgboost) > 0L) {
  
  grafico_series <- dados_series_xgboost |>
    dplyr::select(
      .data$regiao,
      .data$mes,
      Observado = .data$observado,
      `Previsão bruta` = .data$previsto
    ) |>
    tidyr::pivot_longer(
      cols = -c(.data$regiao, .data$mes),
      names_to = "serie",
      values_to = "casos"
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$mes,
        y = .data$casos,
        linetype = .data$serie
      )
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(
      ~regiao,
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = "Validação histórica da previsão bruta do XGBoost",
      subtitle = paste(
        "Anos avaliados:",
        paste(ANOS_BACKTEST, collapse = ", ")
      ),
      x = "Mês",
      y = "Casos",
      linetype = "Série"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom"
    )
  
  ggplot2::ggsave(
    file.path(
      DIR_GRAFICOS,
      "03_observado_previsao_bruta_xgboost.png"
    ),
    grafico_series,
    width = 16,
    height = 10,
    dpi = 300
  )
}

if (nrow(resumo_agregado_xgboost) > 0L) {
  
  dados_totais <- resumo_agregado_xgboost |>
    dplyr::select(
      .data$regiao,
      Observado = .data$total_observado,
      `Previsão bruta` = .data$total_previsto
    ) |>
    tidyr::pivot_longer(
      cols = - .data$regiao,
      names_to = "serie",
      values_to = "total"
    )
  
  grafico_totais <- dados_totais |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$regiao,
        y = .data$total,
        fill = .data$serie
      )
    ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      title = "Totais observados e previstos pelo XGBoost",
      subtitle = paste(
        "Soma dos anos",
        paste(ANOS_BACKTEST, collapse = ", ")
      ),
      x = NULL,
      y = "Casos",
      fill = "Série"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
    )
  
  ggplot2::ggsave(
    file.path(
      DIR_GRAFICOS,
      "04_totais_observados_previsao_bruta.png"
    ),
    grafico_totais,
    width = 12,
    height = 7,
    dpi = 300
  )
}

if (nrow(residuos_mensais_xgboost) > 0L) {
  
  grafico_residuos <- residuos_mensais_xgboost |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$regiao,
        y = .data$residuo
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_boxplot(
      outlier.alpha = 0.5
    ) +
    ggplot2::labs(
      title = "Distribuição dos resíduos mensais do XGBoost",
      subtitle = paste(
        "Resíduo = observado − previsto;",
        "valores positivos indicam subestimação"
      ),
      x = NULL,
      y = "Resíduo mensal"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
    )
  
  ggplot2::ggsave(
    file.path(
      DIR_GRAFICOS,
      "05_distribuicao_residuos_mensais.png"
    ),
    grafico_residuos,
    width = 12,
    height = 7,
    dpi = 300
  )
}

message("\nSeleção e validação histórica concluídas.")
message("Modelo final: ", MODELO_FINAL)
message("Tipo de previsão adotado: previsão bruta")
message("Calibração adicional: não adotada")
message("Resultados salvos em: ", DIR_RESULTADOS)
message(
  "Objeto consolidado do Script 02: ",
  file.path(
    DIR_OBJETOS,
    "resultados_selecao_validacao_modelo_script02.rds"
  )
)

print(decisao_modelo_final)
print(sintese_metodologica)

})

################################################################################
# 3. ESTIMATIVA DE IMPACTO: CONTRAFACTUAL 2026, INTERVALO DE PREVISÃO, CENÁRIOS E AJUSTE POR VIÉS
################################################################################

ARQ_RESULTADOS_SCRIPT02 <- here::here(
  "resultados", "comparacao_modelos_preditivos", "objetos",
  "resultados_selecao_validacao_modelo_script02.rds"
)

stopifnot(file.exists(ARQ_RESULTADOS_SCRIPT02))

resultados_script02 <- readRDS(ARQ_RESULTADOS_SCRIPT02)

# O algoritmo, os parâmetros e a engenharia de atributos usados na estimativa
# oficial de 2026 são exatamente os validados (e escolhidos) na seção 2 —
# evita que o modelo "testado" e o modelo "usado no artigo" divirjam.
MODELO_ESCOLHIDO <- resultados_script02$modelo_final
PARAMETROS_MODELO <- resultados_script02$parametros
RESIDUOS_BACKTESTING <- resultados_script02$residuos_mensais_xgboost

message("Modelo final adotado (seção 2): ", MODELO_ESCOLHIDO)

# Reinicia a semente antes da modelagem final (mesma semente das demais etapas).
set.seed(123)

#### Modelar cenário contrafactual sem vacinação ####

serie_geo <- readRDS(
  file.path(DIR_BASES, "02_serie_mensal_menores_6m.rds")
)

meses_futuros <- seq.Date(
  from = DATA_INICIO_PREVISAO,
  to = DATA_FINAL_FORECAST,
  by = "month"
)

# Engenharia de atributos, ajuste do modelo e previsão recursiva usam as
# MESMAS funções (seção 0 deste arquivo) e o MESMO algoritmo
# (MODELO_ESCOLHIDO) validados por backtesting no Script 02, garantindo que
# a estimativa oficial de 2026 não divirja do que foi de fato testado.
ajustar_prever_regiao <- function(df_regiao) {
  regiao_atual <- as.character(unique(df_regiao$regiao))

  df_treino <- df_regiao |>
    dplyr::select(mes, casos) |>
    dplyr::arrange(mes)

  ultima_data <- max(df_treino$mes, na.rm = TRUE)

  if (ultima_data != DATA_CORTE_TREINO) {
    stop(
      "A série de ", regiao_atual,
      " termina em ", format(ultima_data, "%Y-%m"),
      ", mas deveria terminar em ",
      format(DATA_CORTE_TREINO, "%Y-%m"), "."
    )
  }

  ajuste <- ajustar_prever_recursivo(
    df_treino = df_treino,
    algoritmo = MODELO_ESCOLHIDO,
    datas_futuras = meses_futuros,
    parametros = PARAMETROS_MODELO
  )

  tibble::tibble(
    mes = meses_futuros,
    casos_forecast = ajuste$previsoes,
    regiao = regiao_atual
  )
}

lista_forecast <- serie_geo |>
  dplyr::filter(regiao != "Brasil") |>
  dplyr::group_split(regiao) |>
  purrr::map(ajustar_prever_regiao)

forecast_regioes <- dplyr::bind_rows(lista_forecast)

# O total Brasil é construído pela soma das regiões, evitando um modelo
# nacional independente que não feche com as estimativas regionais.
forecast_brasil <- forecast_regioes |>
  dplyr::group_by(mes) |>
  dplyr::summarise(
    casos_forecast = sum(casos_forecast),
    .groups = "drop"
  ) |>
  dplyr::mutate(regiao = "Brasil")

forecast_geo <- dplyr::bind_rows(
  forecast_brasil,
  forecast_regioes
) |>
  dplyr::filter(
    mes >= as.Date("2026-01-01"),
    mes <= as.Date("2026-12-01")
  ) |>
  dplyr::mutate(
    regiao = factor(regiao, levels = REGIOES_ORDEM)
  ) |>
  dplyr::arrange(regiao, mes)

#### Intervalo de previsão empírico (95%), por região ####

# Construído a partir dos TOTAIS ANUAIS de backtesting do modelo final
# (Script 02: 2019, 2022, 2023, 2024, 2025) — não dos resíduos mensais.
#
# Um teste inicial com resíduos mensais (log1p(observado) - log1p(previsto))
# produzia intervalos irrealistas (ex.: limite superior > 3x o ponto central
# para o Brasil), porque meses de baixa sazonalidade têm contagens próximas
# de zero e qualquer erro absoluto pequeno vira uma razão multiplicativa
# enorme, contaminando o intervalo quando aplicada aos meses de pico.
# Agregando no nível anual (onde o volume nunca é próximo de zero) evita essa
# distorção. A mesma razão anual (escala log1p) é aplicada a todos os meses
# de 2026 da região, preservando o formato sazonal e mantendo os limites
# mensais coerentes com o limite anual.
#
# Limitação: com apenas 5 anos de backtesting por região, os quantis 2,5%/
# 97,5% se aproximam do mínimo/máximo observado — o intervalo deve ser lido
# como uma faixa indicativa da variação histórica do erro, não como um IC
# paramétrico no sentido estrito. Isso é reportado como limitação no artigo.
residuo_anual_log1p <- RESIDUOS_BACKTESTING |>
  dplyr::group_by(regiao, ano_validacao) |>
  dplyr::summarise(
    observado_anual = sum(observado, na.rm = TRUE),
    previsto_anual = sum(previsto, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    residuo_log1p_anual = log1p(observado_anual) - log1p(previsto_anual)
  )

readr::write_csv(
  residuo_anual_log1p,
  file.path(DIR_TABELAS, "03_residuo_anual_log1p_por_regiao.csv")
)

quantis_residuo_log1p <- residuo_anual_log1p |>
  dplyr::group_by(regiao) |>
  dplyr::summarise(
    n_anos = dplyr::n(),
    quantil_inferior = stats::quantile(
      residuo_log1p_anual, probs = IC_PROB_INFERIOR, na.rm = TRUE, names = FALSE
    ),
    quantil_superior = stats::quantile(
      residuo_log1p_anual, probs = IC_PROB_SUPERIOR, na.rm = TRUE, names = FALSE
    ),
    .groups = "drop"
  )

readr::write_csv(
  quantis_residuo_log1p,
  file.path(DIR_TABELAS, "03_quantis_residuo_log1p_por_regiao.csv")
)

forecast_geo <- forecast_geo |>
  dplyr::mutate(regiao_chr = as.character(regiao)) |>
  dplyr::left_join(
    quantis_residuo_log1p,
    by = c("regiao_chr" = "regiao")
  ) |>
  dplyr::mutate(
    casos_forecast_li95 = pmax(
      0, expm1(log1p(casos_forecast) + quantil_inferior)
    ),
    casos_forecast_ls95 = expm1(log1p(casos_forecast) + quantil_superior),
    casos_forecast_li95 = pmin(casos_forecast_li95, casos_forecast),
    casos_forecast_ls95 = pmax(casos_forecast_ls95, casos_forecast)
  ) |>
  dplyr::select(-regiao_chr, -quantil_inferior, -quantil_superior)

#### Cenário contrafactual ajustado por viés histórico (análise de sensibilidade) ####

# O backtesting (Script 02) mostrou subestimação sistemática nos anos de
# crescimento acelerado (2022, 2023, 2024, 2025), sendo 2025 o maior valor
# de toda a série e agora incluído como ano-alvo de backtesting. Para não
# depender de uma única leitura (o contrafactual bruto), esta seção soma um
# cenário ANÁLOGO ajustado pelo viés histórico médio, rodando em paralelo ao
# cenário oficial (bruto).
#
# Fator de viés por região = média geométrica da razão observado/previsto
# nos 5 anos de backtesting = exp(média do resíduo log1p anual). Aplicado de
# forma constante a todos os meses de 2026 da região, preservando o formato
# sazonal — mesma lógica usada no intervalo de previsão.
#
# Este cenário é uma ANÁLISE DE SENSIBILIDADE, não substitui a estimativa
# oficial: usa apenas 4 pontos por região para estimar o fator, então deve
# ser interpretado como indicativo da direção e magnitude plausível do viés,
# não como uma correção definitiva.
fator_vies_regional <- residuo_anual_log1p |>
  dplyr::group_by(regiao) |>
  dplyr::summarise(
    n_anos = dplyr::n(),
    fator_vies = exp(mean(residuo_log1p_anual, na.rm = TRUE)),
    .groups = "drop"
  )

readr::write_csv(
  fator_vies_regional,
  file.path(DIR_TABELAS, "03_fator_vies_regional.csv")
)

forecast_geo <- forecast_geo |>
  dplyr::mutate(regiao_chr = as.character(regiao)) |>
  dplyr::left_join(
    fator_vies_regional |> dplyr::select(regiao, fator_vies),
    by = c("regiao_chr" = "regiao")
  ) |>
  dplyr::mutate(
    casos_forecast_ajustado = casos_forecast * fator_vies
  ) |>
  dplyr::select(-regiao_chr, -fator_vies)

#### Validar horizonte de janeiro a dezembro de 2026 ####

qa_forecast <- forecast_geo |>
  dplyr::mutate(
    regiao = as.character(regiao)
  ) |>
  dplyr::group_by(regiao) |>
  dplyr::summarise(
    primeiro_mes = min(mes),
    ultimo_mes = max(mes),
    numero_meses = dplyr::n(),
    .groups = "drop"
  )

if (
  any(qa_forecast$primeiro_mes != as.Date("2026-01-01")) ||
  any(qa_forecast$ultimo_mes != as.Date("2026-12-01")) ||
  any(qa_forecast$numero_meses != 12)
) {
  print(qa_forecast)
  stop(
    "O forecast não contém exatamente janeiro a dezembro de 2026 ",
    "para todas as regiões."
  )
}

readr::write_csv(
  qa_forecast,
  file.path(DIR_TABELAS, "03_qa_horizonte_forecast.csv")
)

readr::write_csv(
  forecast_geo,
  file.path(DIR_BASES, "03_forecast_sem_vacina_2026.csv")
)

saveRDS(
  forecast_geo,
  file.path(DIR_BASES, "03_forecast_sem_vacina_2026.rds")
)


#### Série histórica e previsão para 2026 ####

# Série observada utilizada no gráfico.
# O recorte a partir de 2020 melhora a leitura visual, mas não altera
# o período usado no treinamento do modelo.
serie_observada_grafico <- serie_geo |>
  dplyr::filter(
    regiao != "Brasil",
    mes >= as.Date("2020-01-01"),
    mes <= DATA_CORTE_TREINO
  ) |>
  dplyr::transmute(
    regiao = as.character(regiao),
    mes,
    casos,
    casos_li95 = NA_real_,
    casos_ls95 = NA_real_,
    tipo = "Observado"
  )

serie_prevista_grafico <- forecast_geo |>
  dplyr::filter(regiao != "Brasil") |>
  dplyr::transmute(
    regiao = as.character(regiao),
    mes,
    casos = casos_forecast,
    casos_li95 = casos_forecast_li95,
    casos_ls95 = casos_forecast_ls95,
    tipo = "Forecasting 2026"
  )

base_historico_forecast <- dplyr::bind_rows(
  serie_observada_grafico,
  serie_prevista_grafico
) |>
  dplyr::filter(
    !is.na(regiao),
    !is.na(mes),
    !is.na(casos),
    !is.na(tipo)
  ) |>
  dplyr::mutate(
    regiao = factor(
      regiao,
      levels = c(
        "Centro-Oeste",
        "Nordeste",
        "Norte",
        "Sudeste",
        "Sul"
      )
    ),
    tipo = factor(
      tipo,
      levels = c(
        "Observado",
        "Forecasting 2026"
      )
    )
  ) |>
  dplyr::arrange(regiao, mes)

# Verificação: cada região deve ter 12 meses previstos em 2026.
qa_historico_forecast <- base_historico_forecast |>
  dplyr::filter(tipo == "Forecasting 2026") |>
  dplyr::group_by(regiao) |>
  dplyr::summarise(
    primeiro_mes = min(mes),
    ultimo_mes = max(mes),
    meses_previstos = dplyr::n(),
    .groups = "drop"
  )

if (
  any(qa_historico_forecast$primeiro_mes != as.Date("2026-01-01")) ||
  any(qa_historico_forecast$ultimo_mes != as.Date("2026-12-01")) ||
  any(qa_historico_forecast$meses_previstos != 12)
) {
  print(qa_historico_forecast)
  stop(
    "A base histórica + forecast não contém os 12 meses de 2026 ",
    "para todas as regiões."
  )
}

grafico_historico_forecast <- ggplot2::ggplot(
  base_historico_forecast,
  ggplot2::aes(
    x = mes,
    y = casos,
    color = tipo,
    group = tipo
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = casos_li95, ymax = casos_ls95),
    fill = "#E74C3C",
    alpha = 0.15,
    color = NA,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(
    ~ regiao,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Observado" = "#2C3E50",
      "Forecasting 2026" = "#E74C3C"
    ),
    drop = TRUE
  ) +
  ggplot2::scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    limits = c(
      as.Date("2020-01-01"),
      as.Date("2026-12-01")
    ),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    ),
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  ggplot2::labs(
    title = "SRAG por VSR em menores de 6 meses",
    subtitle = paste0(
      "Série histórica observada até dezembro de 2025 e ",
      "previsão para 2026, por região geográfica"
    ),
    x = "Ano",
    y = "Casos",
    color = "Tipo de série",
    caption = paste0(
      "Forecast contrafactual para 2026 estimado por ", MODELO_ESCOLHIDO,
      " (modelo validado por backtesting no Script 02). Faixa sombreada: ",
      "intervalo de previsão empírico de 95%, por região. ",
      "O gráfico exibe a série observada a partir de 2020."
    )
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.minor = ggplot2::element_blank()
  )

readr::write_csv(
  base_historico_forecast,
  file.path(
    DIR_TABELAS,
    "03_serie_historica_e_forecast_2026.csv"
  )
)

readr::write_csv(
  qa_historico_forecast,
  file.path(
    DIR_TABELAS,
    "03_qa_serie_historica_e_forecast.csv"
  )
)

ggplot2::ggsave(
  filename = file.path(
    DIR_GRAFICOS,
    "00_serie_historica_e_forecast_2026.png"
  ),
  plot = grafico_historico_forecast,
  width = 13,
  height = 7.5,
  dpi = 300
)

saveRDS(
  list(
    base = base_historico_forecast,
    grafico = grafico_historico_forecast
  ),
  file.path(
    DIR_BASES,
    "03_serie_historica_e_forecast_2026.rds"
  )
)

#### Estimar impacto da vacinação ####

forecast_geo <- readRDS(
  file.path(DIR_BASES, "03_forecast_sem_vacina_2026.rds")
)

impacto_mensal <- forecast_geo |>
  tidyr::crossing(CENARIOS_VACINA) |>
  dplyr::mutate(
    casos_pos_vacina = casos_forecast * (1 - reducao_total),
    casos_evitados = casos_forecast * reducao_total,
    # Os limites do intervalo de previsão do contrafactual são propagados
    # linearmente pelos mesmos fatores de redução — os cenários de VE e
    # cobertura não têm incerteza própria estimada (são valores determinísticos
    # informados pela literatura, ver Tabela 1), então o intervalo final
    # reflete apenas a incerteza da previsão do contrafactual sem vacina.
    casos_pos_vacina_li95 = casos_forecast_li95 * (1 - reducao_total),
    casos_pos_vacina_ls95 = casos_forecast_ls95 * (1 - reducao_total),
    casos_evitados_li95 = casos_forecast_li95 * reducao_total,
    casos_evitados_ls95 = casos_forecast_ls95 * reducao_total
  ) |>
  dplyr::arrange(regiao, cenario, mes)

resumo_anual <- impacto_mensal |>
  dplyr::group_by(
    regiao,
    cenario,
    efetividade,
    cobertura,
    reducao_total
  ) |>
  dplyr::summarise(
    casos_sem_vacina = sum(casos_forecast),
    casos_sem_vacina_li95 = sum(casos_forecast_li95),
    casos_sem_vacina_ls95 = sum(casos_forecast_ls95),
    casos_com_vacina = sum(casos_pos_vacina),
    casos_com_vacina_li95 = sum(casos_pos_vacina_li95),
    casos_com_vacina_ls95 = sum(casos_pos_vacina_ls95),
    casos_evitados = sum(casos_evitados),
    casos_evitados_li95 = sum(casos_evitados_li95),
    casos_evitados_ls95 = sum(casos_evitados_ls95),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    reducao_percentual = 100 * casos_evitados / casos_sem_vacina,
    dplyr::across(
      dplyr::starts_with("casos_"),
      round
    ),
    reducao_percentual = round(reducao_percentual, 1),
    regiao = factor(regiao, levels = REGIOES_ORDEM),
    cenario = factor(
      cenario,
      levels = c("Conservador", "Intermediário", "Otimista")
    )
  ) |>
  dplyr::arrange(regiao, cenario)

readr::write_csv(
  impacto_mensal,
  file.path(DIR_TABELAS, "04_impacto_mensal_2026.csv")
)

readr::write_csv(
  resumo_anual,
  file.path(DIR_TABELAS, "05_resumo_impacto_anual_2026.csv")
)

writexl::write_xlsx(
  list(
    "Cenarios" = CENARIOS_VACINA,
    "Impacto_mensal" = impacto_mensal,
    "Resumo_anual" = resumo_anual
  ),
  path = file.path(DIR_TABELAS, "impacto_vacina_vsr_2026.xlsx")
)

saveRDS(
  impacto_mensal,
  file.path(DIR_BASES, "04_impacto_mensal_2026.rds")
)

#### Análise de sensibilidade: impacto sob contrafactual ajustado por viés ####

# Mesma lógica de impacto (VE x cobertura) aplicada ao contrafactual
# ajustado pelo viés histórico (ver seção acima), em vez do contrafactual
# bruto. Não é combinado com a tabela oficial — fica em arquivos separados
# para deixar explícito que se trata de uma checagem de robustez, não da
# estimativa principal do estudo.
impacto_mensal_ajustado <- forecast_geo |>
  tidyr::crossing(CENARIOS_VACINA) |>
  dplyr::mutate(
    casos_pos_vacina_ajustado = casos_forecast_ajustado * (1 - reducao_total),
    casos_evitados_ajustado = casos_forecast_ajustado * reducao_total
  ) |>
  dplyr::arrange(regiao, cenario, mes)

resumo_anual_ajustado <- impacto_mensal_ajustado |>
  dplyr::group_by(regiao, cenario, efetividade, cobertura, reducao_total) |>
  dplyr::summarise(
    casos_sem_vacina_bruto = sum(casos_forecast),
    casos_sem_vacina_ajustado = sum(casos_forecast_ajustado),
    casos_com_vacina_ajustado = sum(casos_pos_vacina_ajustado),
    casos_evitados_ajustado = sum(casos_evitados_ajustado),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    fator_vies_medio = round(casos_sem_vacina_ajustado / casos_sem_vacina_bruto, 3),
    reducao_percentual_ajustada = round(
      100 * casos_evitados_ajustado / casos_sem_vacina_ajustado, 1
    ),
    dplyr::across(
      c(
        casos_sem_vacina_bruto, casos_sem_vacina_ajustado,
        casos_com_vacina_ajustado, casos_evitados_ajustado
      ),
      round
    ),
    regiao = factor(regiao, levels = REGIOES_ORDEM),
    cenario = factor(
      cenario, levels = c("Conservador", "Intermediário", "Otimista")
    )
  ) |>
  dplyr::arrange(regiao, cenario)

readr::write_csv(
  resumo_anual_ajustado,
  file.path(DIR_TABELAS, "05_resumo_impacto_anual_2026_ajustado_por_vies.csv")
)

message(
  "\nAnálise de sensibilidade (contrafactual ajustado por viés histórico):\n"
)
print(
  resumo_anual_ajustado |>
    dplyr::filter(regiao == "Brasil") |>
    dplyr::select(
      cenario, casos_sem_vacina_bruto, casos_sem_vacina_ajustado,
      fator_vies_medio, casos_evitados_ajustado, reducao_percentual_ajustada
    )
)

#### Gráficos e produtos finais ####

serie_geo <- readRDS(
  file.path(DIR_BASES, "02_serie_mensal_menores_6m.rds")
)

forecast_geo <- readRDS(
  file.path(DIR_BASES, "03_forecast_sem_vacina_2026.rds")
)

impacto_mensal <- readRDS(
  file.path(DIR_BASES, "04_impacto_mensal_2026.rds")
)

cores <- c(
  "Sem vacina" = "#222222",
  "Conservador" = "#D97904",
  "Intermediário" = "#1B7F79",
  "Otimista" = "#1F4E79"
)

base_sem_vacina <- forecast_geo |>
  dplyr::transmute(
    regiao,
    mes,
    cenario = "Sem vacina",
    casos = casos_forecast
  )

base_pos_vacina <- impacto_mensal |>
  dplyr::transmute(
    regiao,
    mes,
    cenario,
    casos = casos_pos_vacina
  )

base_grafico <- dplyr::bind_rows(
  base_sem_vacina,
  base_pos_vacina
) |>
  dplyr::mutate(
    regiao = factor(regiao, levels = REGIOES_ORDEM),
    cenario = factor(
      cenario,
      levels = c(
        "Sem vacina",
        "Conservador",
        "Intermediário",
        "Otimista"
      )
    )
  )

grafico_regioes <- base_grafico |>
  dplyr::filter(regiao != "Brasil") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = mes,
      y = casos,
      color = cenario,
      linetype = cenario
    )
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::facet_wrap(
    ~ regiao,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::scale_color_manual(values = cores) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Sem vacina" = "solid",
      "Conservador" = "longdash",
      "Intermediário" = "dashed",
      "Otimista" = "dotted"
    )
  ) +
  ggplot2::scale_x_date(
    breaks = seq.Date(
      as.Date("2026-01-01"),
      as.Date("2026-12-01"),
      by = "month"
    ),
    date_labels = "%b/%Y",
    limits = c(as.Date("2026-01-01"), as.Date("2026-12-01")),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    )
  ) +
  ggplot2::labs(
    title = "Impacto estimado da vacinação materna contra o VSR",
    subtitle = "Casos de SRAG por VSR em crianças de até 6 meses — projeções para 2026",
    x = "Mês",
    y = "Casos estimados",
    color = "Cenário",
    linetype = "Cenário",
    caption = paste0(
      "Cenários: redução = efetividade × cobertura. ",
      "Estimativas baseadas em cenário contrafactual sem vacinação."
    )
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

grafico_brasil <- base_grafico |>
  dplyr::filter(regiao == "Brasil") |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = mes,
      y = casos,
      color = cenario,
      linetype = cenario
    )
  ) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::scale_color_manual(values = cores) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Sem vacina" = "solid",
      "Conservador" = "longdash",
      "Intermediário" = "dashed",
      "Otimista" = "dotted"
    )
  ) +
  ggplot2::scale_x_date(
    breaks = seq.Date(
      as.Date("2026-01-01"),
      as.Date("2026-12-01"),
      by = "month"
    ),
    date_labels = "%b/%Y",
    limits = c(as.Date("2026-01-01"), as.Date("2026-12-01")),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      big.mark = ".",
      decimal.mark = ","
    )
  ) +
  ggplot2::labs(
    title = "Brasil — impacto estimado da vacinação materna contra o VSR",
    subtitle = "SRAG por VSR em crianças de até 6 meses — projeções para 2026",
    x = "Mês",
    y = "Casos estimados",
    color = "Cenário",
    linetype = "Cenário",
    caption = paste0(
      "Brasil calculado pela soma das cinco macrorregiões. ",
      "Cenários: redução = efetividade × cobertura."
    )
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(DIR_GRAFICOS, "01_impacto_regioes_2026.png"),
  plot = grafico_regioes,
  width = 12,
  height = 9,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(DIR_GRAFICOS, "02_impacto_brasil_2026.png"),
  plot = grafico_brasil,
  width = 11,
  height = 6.5,
  dpi = 300
)

saveRDS(
  list(
    grafico_regioes = grafico_regioes,
    grafico_brasil = grafico_brasil
  ),
  file.path(DIR_BASES, "05_graficos.rds")
)

################################################################################
# FINALIZAÇÃO
################################################################################

message(
  "\nAnálises concluídas com sucesso.\n\n",
  "Arquivos gerados em:\n",
  "  - resultados/bases_processadas\n",
  "  - resultados/tabelas\n",
  "  - resultados/graficos\n"
)

################################################################################
# 4. COMPARAÇÃO EXPLORATÓRIA E NÃO CAUSAL COM OS DADOS OBSERVADOS DE 2026
################################################################################

################################################################################
# 03. COMPARAÇÃO OBSERVADO x CONTRAFACTUAL — SRAG POR VSR EM 2026 (pop<=6meses)
#
# IMPORTANTE — o que este script NÃO faz:
#   Não recalibra nem substitui a estimativa oficial de impacto (Script 01).
#   2026 é o primeiro ano de implementação da vacina (RSVpreF a partir de
#   dez/2025; nirsevimabe a partir de fev/2026): os casos observados em 2026
#   já estão sob efeito parcial da vacinação real. Por isso os dados
#   observados de 2026 NÃO podem ser usados para "validar" ou recalibrar o
#   contrafactual sem-vacina — isso penalizaria o modelo exatamente na
#   direção em que a vacina deveria atuar (observado mais baixo que o
#   previsto), confundindo erro de modelo com efeito real da vacina.
#
# O que este script faz:
#   1. Compara, mês a mês, o contrafactual oficial sem vacina (Script 01,
#      CatBoost + IC95 empírico) com os casos observados reais de 2026,
#      apresentando a diferença como um SINAL PRECOCE E EXPLORATÓRIO do
#      possível impacto da vacinação já em curso — não como uma reestimativa
#      causal do impacto.
#   2. Usa a cobertura vacinal REAL (RNDS/SINASC, painel do Ministério da
#      Saúde, por UF/mês) para checar se a magnitude do sinal acima é
#      compatível com VE x cobertura real — em vez de depender só dos
#      cenários determinísticos fixos (50/70/85%) do Script 01.
#
# Corte de completude dos dados observados de CASOS:
#   A extração de 2026 (BD/INFLUD26-27-07-2026.csv) traz notificações até 26/07/2026 (extração de 27/07/2026) e
#   está sujeita a atraso de digitação do SIVEP-Gripe. Por isso a comparação
#   quantitativa é restrita a jan-mai/2026 (DATA_LIMITE_COMPLETUDE); jun-jul
#   são exibidos no gráfico apenas como "dado preliminar", sem uso nas
#   conclusões.
#
# Sobre a cobertura vacinal real:
#   "Cobertura Vacinal Mensal Acumulada" no painel do MS é medida por coorte
#   de nascimento (mês de nascimento do bebê), não é uma série "ao vivo" —
#   coortes mais recentes (maio/junho) têm menos tempo para a dose da mãe
#   ser vinculada no RNDS até o corte de extração e tendem a subir com o
#   tempo, mesmo viés de atraso que afeta os dados de caso. Alguns UFs
#   mostram cobertura >100% (denominador SINASC pode não bater exatamente
#   com o numerador RNDS) — tratado como limitação, não corrigido.
#
# ESCOPO E PRÉ-REQUISITOS:
#   Este script concentra a comparação exploratória com 2026. O arquivo bruto
#   de 2026 (SRAG de qualquer etiologia) é filtrado aqui para RT-PCR positivo
#   para VSR (PCR_VSR == "1"), replicando o critério já aplicado, em scripts de
#   preparação anteriores (fora deste repositório), às bases 2013-2018 e
#   2019-2025 usadas no Script 01. Os resultados do Script 01 são necessários.
#
# Entradas:
#   funções de harmonização (seção 0 deste arquivo)
#   BD/INFLUD26-27-07-2026.csv (extração SIVEP-Gripe de 2026)
#   BD/cobertura_vacinal_vsr_gestantes_uf_2026.xlsx (painel infoms.saude.gov.br)
#   resultados/bases_processadas/03_forecast_sem_vacina_2026.rds (Script 01)
################################################################################

local({
#### 01. Configuração ####

options(stringsAsFactors = FALSE, scipen = 999)

PACOTES <- c(
  "dplyr", "tidyr", "purrr", "tibble", "lubridate",
  "ggplot2", "scales", "readr", "writexl", "readxl", "here"
)

instalar_ausentes <- function(pacotes) {
  ausentes <- pacotes[
    !vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(ausentes) > 0L) {
    install.packages(ausentes)
  }
}

instalar_ausentes(PACOTES)
invisible(lapply(PACOTES, library, character.only = TRUE))

ARQ_2026 <- here::here("BD", "INFLUD26-27-07-2026.csv")
ARQ_COBERTURA <- here::here(
  "BD", "cobertura_vacinal_vsr_gestantes_uf_2026.xlsx"
)
ARQ_CONTRAFACTUAL <- here::here(
  "resultados", "bases_processadas", "03_forecast_sem_vacina_2026.rds"
)

arquivos_necessarios <- c(
  ARQ_2026, ARQ_COBERTURA, ARQ_CONTRAFACTUAL
)
arquivos_ausentes <- arquivos_necessarios[!file.exists(arquivos_necessarios)]

if (length(arquivos_ausentes) > 0L) {
  stop(
    "Arquivos não encontrados:\n",
    paste0(" - ", arquivos_ausentes, collapse = "\n"),
    "\n\nConfirme se a extração de 2026 está em BD/ e se o Script 01 já ",
    "foi executado (gera o contrafactual oficial)."
  )
}

# (funções de harmonização: seção 0)

FAIXA_ALVO <- "<=6 meses"
REGIOES_ORDEM <- c(
  "Brasil", "Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul"
)

# Último mês de 2026 considerado suficientemente completo para uso
# quantitativo (ver nota de cabeçalho sobre atraso de digitação).
DATA_LIMITE_COMPLETUDE <- as.Date("2026-05-01")
DATA_INICIO_2026 <- as.Date("2026-01-01")
DATA_FIM_EXTRACAO <- as.Date("2026-12-01")

DIR_TABELAS <- here::here("resultados", "tabelas")
DIR_GRAFICOS <- here::here("resultados", "graficos")
dir.create(DIR_TABELAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_GRAFICOS, recursive = TRUE, showWarnings = FALSE)

#### 02. Carregar e harmonizar os dados observados de 2026 ####

# A extração de 2026 traz TODOS os casos de SRAG (qualquer etiologia); o
# filtro por PCR_VSR positivo replica o pré-filtro já aplicado nas bases
# 2013-2018 e 2019-2025 usadas para o treinamento do modelo (ver Script 01).
bruto_2026 <- readr::read_delim(
  ARQ_2026,
  delim = ";",
  locale = readr::locale(encoding = "latin1"),
  col_types = readr::cols(.default = readr::col_character()),
  progress = FALSE
) |>
  dplyr::filter(.data$PCR_VSR == "1")

base_2026_h <- harmonizar_base(bruto_2026) |>
  dplyr::filter(!is.na(.data$DT_SIN_PRI)) |>
  derivar_idade_faixa()

serie_observada_2026_regiao <- base_2026_h |>
  dplyr::filter(
    .data$faixa_idade == FAIXA_ALVO,
    !is.na(.data$REGIAO),
    .data$DT_SIN_PRI >= DATA_INICIO_2026,
    .data$DT_SIN_PRI <= DATA_FIM_EXTRACAO
  ) |>
  dplyr::mutate(mes = lubridate::floor_date(.data$DT_SIN_PRI, "month")) |>
  dplyr::count(REGIAO, mes, name = "casos_observados") |>
  dplyr::rename(regiao = REGIAO)

serie_observada_2026_brasil <- serie_observada_2026_regiao |>
  dplyr::group_by(mes) |>
  dplyr::summarise(casos_observados = sum(casos_observados), .groups = "drop") |>
  dplyr::mutate(regiao = "Brasil")

serie_observada_2026 <- dplyr::bind_rows(
  serie_observada_2026_brasil,
  serie_observada_2026_regiao
)

readr::write_csv(
  serie_observada_2026,
  file.path(DIR_TABELAS, "07_observado_2026_por_regiao.csv")
)

#### 03. Comparar observado x contrafactual oficial (Script 01) ####

contrafactual <- readRDS(ARQ_CONTRAFACTUAL) |>
  dplyr::mutate(regiao = as.character(.data$regiao))

comparacao_mensal <- contrafactual |>
  dplyr::full_join(serie_observada_2026, by = c("regiao", "mes")) |>
  dplyr::mutate(
    regiao = factor(regiao, levels = REGIOES_ORDEM),
    completo = mes <= DATA_LIMITE_COMPLETUDE,
    diferenca = casos_observados - casos_forecast,
    diferenca_percentual = dplyr::if_else(
      casos_forecast > 0,
      100 * diferenca / casos_forecast,
      NA_real_
    ),
    dentro_do_ic95 = dplyr::case_when(
      is.na(casos_observados) ~ NA,
      casos_observados >= casos_forecast_li95 &
        casos_observados <= casos_forecast_ls95 ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  dplyr::arrange(regiao, mes)

readr::write_csv(
  comparacao_mensal,
  file.path(DIR_TABELAS, "07_comparacao_observado_contrafactual_2026.csv")
)

# Resumo quantitativo restrito aos meses considerados completos (jan-mai).
# Meses posteriores (jun-jul) entram apenas no gráfico, como referência
# visual preliminar, e são explicitamente excluídos deste resumo.
resumo_completo <- comparacao_mensal |>
  dplyr::filter(.data$completo, !is.na(.data$casos_observados)) |>
  dplyr::group_by(regiao) |>
  dplyr::summarise(
    n_meses = dplyr::n(),
    casos_observados = sum(casos_observados),
    casos_contrafactual = sum(casos_forecast),
    casos_contrafactual_li95 = sum(casos_forecast_li95),
    casos_contrafactual_ls95 = sum(casos_forecast_ls95),
    diferenca = casos_observados - casos_contrafactual,
    diferenca_percentual = round(
      100 * diferenca / casos_contrafactual, 1
    ),
    meses_dentro_ic95 = sum(dentro_do_ic95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(regiao = factor(regiao, levels = REGIOES_ORDEM)) |>
  dplyr::arrange(regiao)

readr::write_csv(
  resumo_completo,
  file.path(DIR_TABELAS, "07_resumo_observado_vs_contrafactual_jan_mai_2026.csv")
)

writexl::write_xlsx(
  list(
    "Observado_2026_mensal" = serie_observada_2026,
    "Comparacao_mensal_completa" = comparacao_mensal,
    "Resumo_jan_mai_2026" = resumo_completo
  ),
  file.path(DIR_TABELAS, "comparacao_observado_contrafactual_2026.xlsx")
)

#### 04. Gráfico — Brasil e regiões ####

grafico_dados <- comparacao_mensal |>
  dplyr::mutate(
    status = dplyr::case_when(
      is.na(casos_observados) ~ "Sem dado",
      completo ~ "Observado (completo)",
      TRUE ~ "Observado (preliminar, sujeito a revisão)"
    )
  )

grafico_comparacao <- ggplot2::ggplot(
  grafico_dados,
  ggplot2::aes(x = mes)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = casos_forecast_li95, ymax = casos_forecast_ls95),
    fill = "#1B7F79", alpha = 0.15
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = casos_forecast, color = "Contrafactual sem vacina"),
    linewidth = 1
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = casos_observados, color = status, group = 1),
    linewidth = 1
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = casos_observados, color = status),
    size = 1.8, na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~regiao, scales = "free_y", ncol = 3) +
  ggplot2::scale_color_manual(
    values = c(
      "Contrafactual sem vacina" = "#1B7F79",
      "Observado (completo)" = "#1a1a1a",
      "Observado (preliminar, sujeito a revisão)" = "#B0B0B0"
    )
  ) +
  ggplot2::scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(big.mark = ".", decimal.mark = ",")
  ) +
  ggplot2::labs(
    title = "SRAG por VSR em menores de 6 meses — observado x contrafactual, 2026",
    subtitle = paste0(
      "Comparação exploratória e preliminar. NÃO é uma reestimativa do ",
      "impacto vacinal nem uma recalibração do modelo — ver observações ",
      "no cabeçalho do script."
    ),
    x = "Mês de 2026", y = "Casos", color = NULL,
    caption = paste0(
      "Faixa sombreada: IC95 empírico do contrafactual (Script 01). ",
      "Dados completos até ", format(DATA_LIMITE_COMPLETUDE, "%m/%Y"),
      "; jun-jul/2026 exibidos como preliminares (atraso de notificação)."
    )
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(
    DIR_GRAFICOS, "04_observado_vs_contrafactual_2026.png"
  ),
  plot = grafico_comparacao,
  width = 13, height = 8, dpi = 300
)

#### 05. Cobertura vacinal real — checar se VE x cobertura explica o sinal ####

# Painel do Ministério da Saúde (infoms.saude.gov.br): cobertura vacinal
# mensal acumulada de VSR em gestantes, por UF, jan-jun/2026. Exportado
# manualmente ("Baixar Dados (Estados)") e salvo em BD/.
raw_cobertura <- readxl::read_excel(ARQ_COBERTURA, col_names = FALSE)

meses_cobertura <- as.Date(
  as.numeric(raw_cobertura[2, seq(4, 19, 3)]),
  origin = "1899-12-30"
)

metricas <- c("cobertura", "doses", "populacao")
nomes_colunas <- c("regiao_label", "data_extracao", "uf_label")
for (m in seq_along(meses_cobertura)) {
  nomes_colunas <- c(
    nomes_colunas,
    paste0(metricas, "_", format(meses_cobertura[m], "%Y_%m"))
  )
}

dados_cobertura <- raw_cobertura[-(1:3), ]
names(dados_cobertura) <- nomes_colunas

dados_cobertura <- dados_cobertura |>
  dplyr::mutate(
    uf_sigla = ifelse(
      is.na(uf_label), "BRASIL", trimws(sub("^[0-9]+ - ", "", uf_label))
    )
  ) |>
  dplyr::filter(!(uf_sigla == "BRASIL" & dplyr::row_number() > 1))

cobertura_longa <- purrr::map_dfr(seq_along(meses_cobertura), function(i) {
  suf <- format(meses_cobertura[i], "%Y_%m")
  dados_cobertura |>
    dplyr::transmute(
      uf_sigla,
      mes = meses_cobertura[i],
      doses = as.numeric(.data[[paste0("doses_", suf)]]),
      populacao = as.numeric(.data[[paste0("populacao_", suf)]])
    )
})

cobertura_brasil <- cobertura_longa |>
  dplyr::filter(uf_sigla == "BRASIL") |>
  dplyr::transmute(regiao = "Brasil", mes, doses, populacao)

cobertura_regional <- cobertura_longa |>
  dplyr::filter(uf_sigla != "BRASIL") |>
  dplyr::mutate(regiao = unname(uf_regiao[uf_sigla])) |>
  dplyr::filter(!is.na(regiao)) |>
  dplyr::group_by(regiao, mes) |>
  dplyr::summarise(
    doses = sum(doses, na.rm = TRUE),
    populacao = sum(populacao, na.rm = TRUE),
    .groups = "drop"
  )

cobertura_regiao_mes <- dplyr::bind_rows(cobertura_brasil, cobertura_regional) |>
  dplyr::mutate(
    cobertura_pct_bruta = 100 * doses / populacao,
    # Alguns UFs somam >100% (denominador SINASC não bate exatamente com o
    # numerador RNDS); limitado a 100% só para uso na estimativa abaixo.
    cobertura_frac = pmin(1, doses / populacao),
    regiao = factor(regiao, levels = REGIOES_ORDEM)
  ) |>
  dplyr::arrange(regiao, mes)

readr::write_csv(
  cobertura_regiao_mes,
  file.path(DIR_TABELAS, "08_cobertura_vacinal_vsr_por_regiao_mes.csv")
)

# VE ancorada nos mesmos três cenários do Script 01 (60/70/75%), aplicada à
# cobertura REAL observada em cada região-mês (em vez da cobertura fixa dos
# cenários) — checa se a MAGNITUDE do gap observado-vs-contrafactual (seção
# 03 acima) é compatível com o que a cobertura real do programa já entrega.
VE_CENARIOS <- c(Conservador = 0.60, Intermediario = 0.70, Otimista = 0.75)

comparacao_cobertura <- contrafactual |>
  dplyr::inner_join(cobertura_regiao_mes, by = c("regiao", "mes")) |>
  dplyr::filter(mes <= DATA_LIMITE_COMPLETUDE) |>
  dplyr::left_join(
    serie_observada_2026, by = c("regiao", "mes")
  ) |>
  tidyr::crossing(tibble::enframe(VE_CENARIOS, "cenario_ve", "efetividade")) |>
  dplyr::mutate(
    reducao_estimada = efetividade * cobertura_frac,
    casos_com_vacina_estimado = casos_forecast * (1 - reducao_estimada),
    casos_evitados_estimado = casos_forecast * reducao_estimada
  )

resumo_cobertura_real <- comparacao_cobertura |>
  dplyr::group_by(regiao, cenario_ve, efetividade) |>
  dplyr::summarise(
    n_meses = dplyr::n(),
    cobertura_media_pct = round(100 * mean(cobertura_frac), 1),
    casos_contrafactual = round(sum(casos_forecast)),
    casos_evitados_estimado = round(sum(casos_evitados_estimado)),
    casos_com_vacina_estimado = round(sum(casos_com_vacina_estimado)),
    casos_observados = sum(casos_observados),
    reducao_percentual_estimada = round(
      100 * sum(casos_evitados_estimado) / sum(casos_forecast), 1
    ),
    reducao_percentual_observada = round(
      100 * (sum(casos_forecast) - sum(casos_observados)) / sum(casos_forecast), 1
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(regiao, cenario_ve)

readr::write_csv(
  resumo_cobertura_real,
  file.path(
    DIR_TABELAS, "08_resumo_ve_x_cobertura_real_jan_mai_2026.csv"
  )
)

#### 06. Síntese no console ####

message(
  "\nComparação observado x contrafactual concluída (jan-mai/2026, dados completos).\n"
)
print(resumo_completo)

message(
  "\nRedução estimada por VE x cobertura VACINAL REAL (RNDS/SINASC) vs. ",
  "redução realmente observada — Brasil, jan-mai/2026:\n"
)
print(
  resumo_cobertura_real |>
    dplyr::filter(regiao == "Brasil") |>
    dplyr::select(
      cenario_ve, efetividade, cobertura_media_pct,
      reducao_percentual_estimada, reducao_percentual_observada
    )
)

message(
  "\nAVISO: esta comparação é exploratória. Os casos observados de 2026 já ",
  "ocorrem sob vacinação real parcial (RSVpreF desde dez/2025, nirsevimabe ",
  "desde fev/2026); a diferença observada NÃO deve ser interpretada como ",
  "uma reestimativa causal do impacto vacinal. A cobertura vacinal real por ",
  "coorte de nascimento tende a subir com o tempo (atraso de vinculação no ",
  "RNDS), então meses mais recentes (maio) provavelmente estão subestimados."
)

})

################################################################################
# 5. FIGURAS DO MANUSCRITO (padrão da revista: 107 mm de largura, 300 dpi)
################################################################################

# Série observada completa para as figuras: 2013-2025 (seção 1) e 2026 (seção 4).
serie_observada_completa <- dplyr::bind_rows(
  readRDS(file.path(DIR_BASES, "02_serie_mensal_completa_2013_2025.rds")) |>
    dplyr::mutate(regiao = as.character(regiao)) |>
    dplyr::select(regiao, mes, casos),
  readr::read_csv(
    here::here("resultados", "tabelas", "07_observado_2026_por_regiao.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::transmute(regiao, mes = as.Date(mes), casos = casos_observados)
)

# ---- Figure 1 — validação (backtesting) e contrafactual 2026
local({
###############################################################################
# Figure 1 (artigo, em ingles; padrao The Lancet Regional Health - Americas)
# Backtesting 2019 e 2022-2025 (CatBoost) + contrafactual 2026 com IP 95%, por regiao.
# Replota saidas ja geradas (sem rodar nenhuma analise nova):
#   - resultados/bases_processadas/02_serie_mensal_completa_2013_2026.rds (observado)
#   - resultados/comparacao_modelos_preditivos/tabelas/01_previsoes_backtesting_algoritmos.csv
#   - resultados/bases_processadas/03_forecast_sem_vacina_2026.rds (contrafactual + IP 95%)
# Padrao da revista: 107 mm de largura, 300 dpi, sem titulo/subtitulo/rodape na arte,
# fonte uniforme, sem moldura nos paineis (titulo e nota vao na legenda do manuscrito).
# Saidas: resultados/graficos/figura1_validacao_lancet_en.png (300 dpi) e .pdf (vetorial)
###############################################################################
suppressMessages({library(dplyr); library(ggplot2); library(scales); library(readr); library(here)})

EN <- c("Brasil"="Brazil","Norte"="North","Nordeste"="Northeast",
        "Centro-Oeste"="Centre-West","Sudeste"="Southeast","Sul"="South")
ORDEM <- c("Brazil","North","Northeast","Centre-West","Southeast","South")
reg <- function(x) factor(EN[as.character(x)], levels = ORDEM)

obs <- serie_observada_completa |>
  transmute(regiao = reg(regiao), mes, casos)

bt <- read_csv(here("resultados","comparacao_modelos_preditivos","tabelas",
                    "01_previsoes_backtesting_algoritmos.csv"), show_col_types = FALSE) |>
  filter(algoritmo == "CATBOOST", status == "ok") |>
  transmute(regiao = reg(regiao), ano_validacao, mes = as.Date(mes), previsto)

fc <- readRDS(here("resultados","bases_processadas","03_forecast_sem_vacina_2026.rds")) |>
  transmute(regiao = reg(regiao), mes = as.Date(mes), casos_forecast,
            li = casos_forecast_li95, ls = casos_forecast_ls95)

bandas <- tibble(ano = c(2019, 2022:2025)) |>
  mutate(xmin = as.Date(paste0(ano, "-01-01")), xmax = as.Date(paste0(ano + 1, "-01-01")),
         tom = c("a", "a", "b", "a", "b"))

COR <- c("Observed" = "#1a1a1a", "CatBoost (backtest forecast)" = "#D62728",
         "2026 counterfactual (no vaccination)" = "#1F4E79")

p <- ggplot() +
  geom_rect(data = bandas, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = tom),
            alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c(a = "grey88", b = "grey78")) +
  geom_ribbon(data = fc, aes(mes, ymin = li, ymax = ls), fill = "#1F4E79", alpha = 0.18) +
  geom_line(data = obs, aes(mes, casos, color = "Observed"), linewidth = 0.35) +
  geom_line(data = bt, aes(mes, previsto, color = "CatBoost (backtest forecast)",
                           group = ano_validacao), linewidth = 0.35) +
  geom_line(data = fc, aes(mes, casos_forecast, color = "2026 counterfactual (no vaccination)"),
            linewidth = 0.5) +
  facet_wrap(~regiao, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_color_manual(values = COR, breaks = names(COR)) +
  scale_x_date(breaks = as.Date(paste0(2013:2027, "-01-01")), date_labels = "%Y",
               expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Year", y = "Monthly cases (up to 6 months)", color = NULL) +
  theme_minimal(base_size = 6.5, base_family = "Helvetica") +
  theme(legend.position = "bottom", legend.key.width = unit(8, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0),
        strip.text = element_text(face = "bold", size = 6.5),
        axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 5.5),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), plot.background = element_rect(fill = "white", color = NA),
        panel.spacing.y = unit(1.2, "mm"))

dir.create(here("resultados","graficos"), showWarnings = FALSE, recursive = TRUE)
W <- 107 / 25.4; H <- 175 / 25.4
ggsave(here("resultados","graficos","figura1_validacao_lancet_en.png"), p,
       width = W, height = H, dpi = 300, bg = "white")
ggsave(here("resultados","graficos","figura1_validacao_lancet_en.pdf"), p,
       width = W, height = H, device = "pdf")
cat("salvo\n")

})

# ---- Figure 2 — cenários de vacinação por região
local({
###############################################################################
# Figura de cenarios de vacinacao (artigo, em ingles) - eixo Y fixo, Norte por ultimo
# Sugestao do coautor (W. Araujo): eixo Y comum a todos os paineis. Replota
# resultados/tabelas/04_impacto_mensal_2026.csv, sem rodar nenhuma analise nova.
# Saida: resultados/graficos/figura_cenarios_eixoY_fixo_en.png
###############################################################################
Sys.setlocale("LC_TIME", "en_US.UTF-8")
suppressMessages({library(dplyr); library(ggplot2); library(scales); library(readr); library(here)})
EN <- c("Nordeste"="Northeast","Centro-Oeste"="Centre-West","Sudeste"="Southeast","Sul"="South","Norte"="North")
d <- read_csv(here("resultados","tabelas","04_impacto_mensal_2026.csv"), show_col_types = FALSE) |>
  filter(regiao != "Brasil") |>
  mutate(regiao = factor(EN[regiao], levels = c("Northeast","Centre-West","Southeast","South","North")),
         cenario = factor(cenario, levels = c("Conservador","Intermediário","Otimista"),
                          labels = c("Conservative","Intermediate","Optimistic")))
base <- d |> distinct(regiao, mes, casos_forecast)
p <- ggplot() +
  geom_line(data = base, aes(mes, casos_forecast), color = "black", linewidth = 0.9) +
  geom_line(data = d, aes(mes, casos_pos_vacina, color = cenario, linetype = cenario), linewidth = 0.8) +
  facet_wrap(~regiao, ncol = 2) +
  scale_y_continuous(limits = c(0, 2200), breaks = seq(0, 2000, 500), labels = label_number(big.mark = ",")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Conservative"="#D98E04","Intermediate"="#1F7A6E","Optimistic"="#1F4E79")) +
  scale_linetype_manual(values = c("Conservative"="longdash","Intermediate"="dashed","Optimistic"="dotted")) +
  labs(x = "Month of 2026", y = "Estimated hospitalizations", color = "Scenario", linetype = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
dir.create(here("resultados","graficos"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("resultados","graficos","figura_cenarios_eixoY_fixo_en.png"), p, width = 9, height = 9, dpi = 300)
cat("salvo\n")

})

# ---- Figure 3 — observado vs contrafactual 2026
local({
###############################################################################
# Figure 3 (artigo, em ingles; padrao The Lancet Regional Health - Americas)
# Observado 2026 vs contrafactual (IP 95%), por regiao. Comparacao exploratoria.
# Replota saidas ja geradas (sem rodar nenhuma analise nova):
#   - resultados/bases_processadas/03_forecast_sem_vacina_2026.rds
#   - resultados/bases_processadas/02_serie_mensal_completa_2013_2026.rds
# Padrao da revista: 107 mm de largura, 300 dpi, sem titulo/rodape na arte.
# Saidas: resultados/graficos/figura3_observado_vs_contrafactual_lancet_en.png (.pdf)
###############################################################################
suppressMessages({library(dplyr); library(ggplot2); library(scales); library(here)})
Sys.setlocale("LC_TIME", "en_US.UTF-8")

EN <- c("Brasil"="Brazil","Norte"="North","Nordeste"="Northeast",
        "Centro-Oeste"="Centre-West","Sudeste"="Southeast","Sul"="South")
ORDEM <- c("Brazil","North","Northeast","Centre-West","Southeast","South")
reg <- function(x) factor(EN[as.character(x)], levels = ORDEM)

fc <- readRDS(here("resultados","bases_processadas","03_forecast_sem_vacina_2026.rds")) |>
  transmute(regiao = reg(regiao), mes = as.Date(mes), casos_forecast,
            li = casos_forecast_li95, ls = casos_forecast_ls95)

obs <- serie_observada_completa |>
  filter(mes >= as.Date("2026-01-01")) |>
  transmute(regiao = reg(regiao), mes = as.Date(mes), casos,
            status = ifelse(mes <= as.Date("2026-05-01"), "Observed (complete)",
                            "Observed (preliminary, subject to revision)"))

COR <- c("Counterfactual (no vaccination)" = "#1B7F73", "Observed (complete)" = "#1a1a1a",
         "Observed (preliminary, subject to revision)" = "grey60")

p <- ggplot() +
  geom_ribbon(data = fc, aes(mes, ymin = li, ymax = ls), fill = "#1B7F73", alpha = 0.18) +
  geom_line(data = fc, aes(mes, casos_forecast, color = "Counterfactual (no vaccination)"), linewidth = 0.45) +
  geom_line(data = filter(obs, status == "Observed (complete)"),
            aes(mes, casos, color = "Observed (complete)"), linewidth = 0.4) +
  geom_line(data = filter(obs, mes >= as.Date("2026-05-01"), mes <= as.Date("2026-06-01")),
            aes(mes, casos, color = "Observed (complete)"), linewidth = 0.4) +
  geom_point(data = filter(obs, status == "Observed (complete)"),
             aes(mes, casos, color = "Observed (complete)"), size = 0.8) +
  geom_line(data = filter(obs, status != "Observed (complete)"),
            aes(mes, casos, color = "Observed (preliminary, subject to revision)"), linewidth = 0.4) +
  geom_point(data = filter(obs, status != "Observed (complete)"),
             aes(mes, casos, color = "Observed (preliminary, subject to revision)"), size = 0.8) +
  facet_wrap(~regiao, ncol = 2, scales = "free_y") +
  expand_limits(y = 0) +
  scale_color_manual(values = COR, breaks = names(COR),
                     guide = guide_legend(nrow = 3, override.aes = list(
                       linetype = "solid", shape = c(NA, 16, 16)))) +
  scale_x_date(breaks = as.Date(paste0("2026-", c("01","04","07","10"), "-01")),
               date_labels = "%b", expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "Month of 2026", y = "Monthly hospitalizations (up to 6 months)", color = NULL) +
  theme_minimal(base_size = 6.5, base_family = "Helvetica") +
  theme(legend.position = "bottom", legend.key.width = unit(7, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-2, 0, 0, 0),
        legend.text = element_text(size = 6), legend.key.height = unit(3, "mm"),
        strip.text = element_text(face = "bold", size = 6.5),
        axis.text = element_text(size = 5.5),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.spacing = unit(2.5, "mm"),
        plot.background = element_rect(fill = "white", color = NA))

dir.create(here("resultados","graficos"), showWarnings = FALSE, recursive = TRUE)
W <- 107 / 25.4; H <- 130 / 25.4
ggsave(here("resultados","graficos","figura3_observado_vs_contrafactual_lancet_en.png"), p,
       width = W, height = H, dpi = 300, bg = "white")
ggsave(here("resultados","graficos","figura3_observado_vs_contrafactual_lancet_en.pdf"), p,
       width = W, height = H, device = "pdf")
cat("salvo\n")

})

# ---- Figure 4 — redução observada vs esperada (VE × cobertura real)
local({
###############################################################################
# Figure 4 (artigo, em ingles; padrao The Lancet Regional Health - Americas)
# Reducao observada vs reducao esperada (VE 60-75% x cobertura real), jan-mai/2026, por regiao.
# Replota saida ja gerada (sem rodar nenhuma analise nova):
#   - resultados/tabelas/08_resumo_ve_x_cobertura_real_jan_mai_2026.csv
# Padrao da revista: 107 mm de largura, 300 dpi, sem titulo/rodape na arte.
# Valores negativos = casos observados ACIMA do contrafactual.
# Saidas: resultados/graficos/figura4_observado_vs_esperado_lancet_en.png (.pdf)
###############################################################################
suppressMessages({library(dplyr); library(ggplot2); library(readr); library(here)})

EN <- c("Brasil"="Brazil","Norte"="North","Nordeste"="Northeast",
        "Centro-Oeste"="Centre-West","Sudeste"="Southeast","Sul"="South")
ORDEM <- c("Brazil","North","Northeast","Centre-West","Southeast","South")

d <- read_csv(here("resultados","tabelas","08_resumo_ve_x_cobertura_real_jan_mai_2026.csv"),
              show_col_types = FALSE) |>
  mutate(regiao = factor(EN[regiao], levels = ORDEM))

faixa <- d |>
  group_by(regiao) |>
  summarise(cobertura = first(cobertura_media_pct),
            ymin = min(reducao_percentual_estimada), ymax = max(reducao_percentual_estimada),
            obs = first(reducao_percentual_observada), .groups = "drop") |>
  mutate(x = as.integer(regiao),
         rotulo = sprintf("%s\n%.1f%%", regiao, cobertura),
         obs_txt = sprintf("%s%.0f%%", ifelse(obs < 0, "−", ""), abs(obs)))

p <- ggplot(faixa) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  geom_rect(aes(xmin = x - 0.15, xmax = x + 0.15, ymin = ymin, ymax = ymax,
                fill = "Expected reduction (VE 60–75% × real coverage)")) +
  geom_point(aes(x, obs, shape = "Observed reduction"), size = 1.8, color = "black") +
  geom_text(aes(x + 0.2, obs, label = obs_txt), size = 2.1, hjust = 0, fontface = "bold") +
  scale_fill_manual(values = c("Expected reduction (VE 60–75% × real coverage)" = "#2A9D8F")) +
  scale_shape_manual(values = c("Observed reduction" = 16)) +
  scale_x_continuous(breaks = faixa$x, labels = faixa$rotulo, expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(-95, 95), breaks = seq(-80, 80, 40),
                     labels = function(v) paste0(ifelse(v < 0, "−", ""), abs(v), "%")) +
  labs(x = "Region (real vaccination coverage, January–May 2026)",
       y = "Reduction in hospitalizations vs. counterfactual", fill = NULL, shape = NULL) +
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_minimal(base_size = 6.5, base_family = "Helvetica") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"),
        axis.text.x = element_text(size = 5.8, lineheight = 1.05),
        axis.text.y = element_text(size = 5.8), axis.title.x = element_text(size = 6),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_blank(), plot.background = element_rect(fill = "white", color = NA))

dir.create(here("resultados","graficos"), showWarnings = FALSE, recursive = TRUE)
W <- 107 / 25.4; H <- 85 / 25.4
ggsave(here("resultados","graficos","figura4_observado_vs_esperado_lancet_en.png"), p,
       width = W, height = H, dpi = 300, bg = "white")
ggsave(here("resultados","graficos","figura4_observado_vs_esperado_lancet_en.pdf"), p,
       width = W, height = H, device = "pdf")
cat("salvo\n")

})

# ---- Supplementary Figure S1 — contrafactual bruto vs ajustado por viés
local({
###############################################################################
# Supplementary Figure S1 (artigo, em ingles; padrao The Lancet Regional Health - Americas)
# Contrafactual bruto vs ajustado por vies vs observado em 2026 - Brasil.
# Replota saidas ja geradas (sem rodar nenhuma analise nova):
#   - resultados/bases_processadas/03_forecast_sem_vacina_2026.rds (bruto, IP 95%, ajustado)
#   - resultados/bases_processadas/02_serie_mensal_completa_2013_2026.rds (observado)
# Padrao da revista: 107 mm de largura, 300 dpi, sem titulo/subtitulo/rodape na arte.
# Saidas: resultados/graficos/figuraS1_vies_lancet_en.png (300 dpi) e .pdf
###############################################################################
suppressMessages({library(dplyr); library(ggplot2); library(scales); library(here)})
Sys.setlocale("LC_TIME", "en_US.UTF-8")

fc <- readRDS(here("resultados","bases_processadas","03_forecast_sem_vacina_2026.rds")) |>
  filter(regiao == "Brasil") |>
  transmute(mes = as.Date(mes), bruto = casos_forecast, li = casos_forecast_li95,
            ls = casos_forecast_ls95, ajustado = casos_forecast_ajustado)

obs <- serie_observada_completa |>
  filter(regiao == "Brasil", mes >= as.Date("2026-01-01")) |>
  transmute(mes = as.Date(mes), casos,
            status = ifelse(mes <= as.Date("2026-05-01"), "Observed (complete)",
                            "Observed (preliminary, subject to revision)"))

COR <- c("Raw counterfactual" = "#1B9E8F", "Bias-adjusted counterfactual" = "#0B4F4A",
         "Observed (complete)" = "#1a1a1a", "Observed (preliminary, subject to revision)" = "grey60")

p <- ggplot() +
  geom_ribbon(data = fc, aes(mes, ymin = li, ymax = ls), fill = "#1B9E8F", alpha = 0.18) +
  geom_line(data = fc, aes(mes, bruto, color = "Raw counterfactual"), linewidth = 0.6) +
  geom_line(data = fc, aes(mes, ajustado, color = "Bias-adjusted counterfactual"),
            linewidth = 0.6, linetype = "dashed") +
  geom_line(data = filter(obs, status == "Observed (complete)"),
            aes(mes, casos, color = "Observed (complete)"), linewidth = 0.5) +
  geom_point(data = filter(obs, status == "Observed (complete)"),
             aes(mes, casos, color = "Observed (complete)"), size = 1.2) +
  geom_line(data = filter(obs, status != "Observed (complete)"),
            aes(mes, casos, color = "Observed (preliminary, subject to revision)"), linewidth = 0.5) +
  geom_point(data = filter(obs, status != "Observed (complete)"),
             aes(mes, casos, color = "Observed (preliminary, subject to revision)"), size = 1.2) +
  scale_color_manual(values = COR, breaks = names(COR),
                     guide = guide_legend(nrow = 2, override.aes = list(
                       linetype = c("solid", "dashed", "solid", "solid"),
                       shape = c(NA, NA, 16, 16)))) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "Month of 2026", y = "Monthly hospitalizations (up to 6 months)", color = NULL) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = "bottom", legend.key.width = unit(8, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 6),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        plot.background = element_rect(fill = "white", color = NA))

dir.create(here("resultados","graficos"), showWarnings = FALSE, recursive = TRUE)
W <- 107 / 25.4; H <- 85 / 25.4
ggsave(here("resultados","graficos","figuraS1_vies_lancet_en.png"), p,
       width = W, height = H, dpi = 300, bg = "white")
ggsave(here("resultados","graficos","figuraS1_vies_lancet_en.pdf"), p, width = W, height = H, device = "pdf")
cat("salvo\n")

})

message("\nReprodução concluída. Figuras em resultados/graficos/; tabelas em resultados/tabelas/.")
