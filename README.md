# IMPACTO_VACINA_VSR_BRA — Código de reprodução do artigo

Código de reprodução do estudo **"Potential impact of maternal vaccination against respiratory syncytial virus on severe acute respiratory infections in infants up to six months of age in Brazil, 2026"** (manuscrito em preparação).

*English summary: see the [last section](#english-summary).*

## Escopo

Este repositório publica **um único arquivo**, [`reproducao_artigo_vsr.R`](reproducao_artigo_vsr.R), com as **análises de estimativa de impacto** do artigo: validação dos algoritmos, projeção contrafactual de 2026, cenários de vacinação, análise de sensibilidade, comparação exploratória com 2026 e figuras do manuscrito.

Ele **não inclui a preparação das bases brutas**. Os arquivos de `BD/` foram gerados por scripts de preparação anteriores a esta análise, mantidos fora deste repositório, que aplicam os critérios de elegibilidade aos dados públicos (em especial, a seleção de casos de SRAG com RT-PCR positivo para VSR). Esses filtros prévios são pré-requisito para reproduzir os resultados (ver "Dados").

## Objetivo

Estimar o impacto potencial da vacinação materna contra o vírus sincicial respiratório (VSR) sobre as hospitalizações por síndrome respiratória aguda grave (SRAG) com VSR confirmado por RT-PCR em crianças de até seis meses no Brasil em 2026, a partir de um modelo preditivo ex-ante treinado com dados de vigilância de 2013 a dezembro de 2025, e comparar de forma exploratória e **não causal** as projeções com os dados observados de janeiro a maio de 2026 e com a cobertura vacinal real.

## Estrutura

```text
IMPACTO_VACINA_VSR_BRA_ARTIGO/
├── reproducao_artigo_vsr.R      # arquivo único com todas as análises e figuras
├── BD/                          # bases de entrada (NÃO versionadas)
├── resultados/                  # produtos gerados (NÃO versionados)
├── impacto_vacina_vsr.Rproj
└── README.md
```

## Como executar

1. Abra `impacto_vacina_vsr.Rproj` (o arquivo usa `here::here()`; a raiz do projeto é a pasta do `.Rproj`).
2. Coloque as bases em `BD/` (ver "Dados").
3. Execute `source("reproducao_artigo_vsr.R")` do início ao fim. Os pacotes ausentes são instalados automaticamente (`install.packages`). A execução completa levou cerca de 30 segundos em um computador pessoal (macOS, Apple Silicon).

O arquivo é organizado em seções, executadas nesta ordem:

| Seção | Conteúdo |
|---|---|
| 0 | Configuração geral, pacotes e funções compartilhadas (modelagem e harmonização) |
| 1 | Preparação das séries mensais: casos de até 6 meses por macrorregião, 2013–2025 |
| 2 | Validação metodológica: backtesting (2019, 2022–2025) e seleção do algoritmo |
| 3 | Estimativa de impacto: contrafactual 2026, intervalo de previsão de 95%, cenários e ajuste por viés |
| 4 | Comparação exploratória e não causal com 2026 e com a cobertura vacinal real |
| 5 | Figuras do manuscrito |

**Verificação.** Uma execução completa a partir do zero (pasta vazia, apenas as quatro bases de `BD/`) reproduziu **byte a byte** as 30 tabelas CSV geradas nas análises do artigo (série mensal, previsões do backtesting, métricas e ranking dos algoritmos, contrafactual de 2026 com intervalos, cenários de impacto, ajuste por viés, comparação com 2026 e cobertura real).

## Métodos em síntese (conforme o manuscrito)

- **Dados:** SINAN Influenza Web (2009–2018) e SIVEP-Gripe (2019 em diante); série mensal contínua de 2013 a dezembro de 2025. Casos elegíveis: SRAG com RT-PCR positivo para VSR em crianças de até 6 meses no início dos sintomas.
- **Algoritmos comparados:** XGBoost, LightGBM e CatBoost, por backtesting temporal de janela expansível nos anos **2019 e 2022–2025** (2020–2021 ficam no treino, mas fora da avaliação). Ranking composto por RMSE, MAE, sMAPE, WAPE, MASE e razão previsto/observado.
- **Seleção:** o CatBoost teve o melhor ranking em 4 de 6 unidades (Norte, Nordeste, Sudeste e Sul); o XGBoost teve melhor ranking no Brasil e no Centro-Oeste. O CatBoost foi adotado como **algoritmo único para todas as regiões**, por consistência operacional. A previsão é bruta, sem fatores de calibração.
- **Modelagem regional:** um modelo por macrorregião; o total do Brasil (estimativa pontual) é obtido pela soma das cinco regiões. Atributos: tendência, ano, termos sazonais trigonométricos, indicadores de período, defasagens (1–3 e 12), médias móveis de 3 e 6 meses, diferença mensal e razão entre defasagens.
- **Incerteza:** intervalo de previsão empírico de 95% por região, a partir dos percentis 2,5 e 97,5 dos resíduos anuais `log1p(observado) − log1p(previsto)` do backtesting, aplicados multiplicativamente a todos os meses de 2026.
- **Cenários determinísticos (efetividade × cobertura):** conservador 60% × 50% = 30%; intermediário 70% × 70% = 49%; otimista 75% × 85% = 63,75%.
- **Análise de sensibilidade:** contrafactual ajustado por viés, com fator igual à média geométrica da razão observado/previsto nos cinco anos de backtesting.
- **Comparação exploratória com 2026:** janeiro–maio de 2026 (período com notificação considerada completa), cobertura vacinal real por região (painel do Ministério da Saúde) combinada com efetividade de 60–75%. Não recalibra nem valida o modelo.

## Dados (`BD/`, não versionada)

Todas as bases são **públicas** e os registros individuais são anonimizados antes da divulgação. Nenhum dado individual além dos publicados foi usado. As bases não são versionadas neste repositório; para reproduzir, baixe-as nas fontes abaixo e coloque-as em `BD/`.

### Fontes

| Base | Arquivo em `BD/` | Usada em | Onde obter |
|---|---|---|---|
| SINAN Influenza Web (SRAG), 2013–2018 | `base_VSR_2013_2018.RData` (objeto `base_VSR_2013_2018`) | Seção 1 | Portal de Dados Abertos do SUS (https://dadosabertos.saude.gov.br, busca por "SRAG") |
| SIVEP-Gripe (SRAG), 2019–2025, versão consolidada | `base_DEF_VSR_2019_2025.RData` (objeto `base_DEF_VSR_2019_2025`) | Seção 1 | OpenDATASUS, conjunto "SRAG 2019 a 2026" (https://dadosabertos.saude.gov.br/dataset/srag-2019-a-2026) |
| SIVEP-Gripe (SRAG), 2026 | `INFLUD26-27-07-2026.csv` | Seção 4 | Mesmo conjunto acima; extração de **27/07/2026** (notificações até 26/07/2026; separador `;`, codificação latin1) |
| Cobertura vacinal contra VSR em gestantes, por UF de residência, janeiro–junho de 2026 | `cobertura_vacinal_vsr_gestantes_uf_2026.xlsx` | Seção 4 | Painel de cobertura vacinal do Ministério da Saúde (https://infoms.saude.gov.br), "Cobertura Vacinal — Calendário Nacional — Gestantes — Residência" |

### Pré-filtro aplicado às bases do treino

As duas bases `.RData` contêm **apenas casos com RT-PCR positivo para VSR**, selecionados em scripts de preparação anteriores a esta análise, que não fazem parte deste repositório (a seleção é pré-requisito, e não é refeita aqui para as bases de treino):

- 2013–2018 (SINAN Influenza Web): `RES_VSR == 1` (16.374 registros, todas as idades);
- 2019–2025 (SIVEP-Gripe): `PCR_VSR == 1` (128.313 registros, todas as idades).

A seção 4 reaplica o mesmo critério ao arquivo bruto de 2026, que traz SRAG de qualquer etiologia: `PCR_VSR == "1"`. A restrição a crianças de até 6 meses é feita **dentro do arquivo** (faixa etária `<=6 meses`, ver abaixo). Para reconstruir as bases a partir dos arquivos públicos, aplique esses mesmos filtros de RT-PCR positivo (e os demais critérios de elegibilidade da etapa de preparação) e salve os objetos com os nomes acima. Os registros de SRAG nesses sistemas são, por definição, de pacientes hospitalizados; mais de 99% dos registros das duas bases têm `HOSPITAL = 1`.

### Variáveis das bases usadas

A harmonização (seção 0, funções `harmonizar_base` e `derivar_idade_faixa`) usa o mesmo procedimento para as três extrações.

| Variável | Significado | Como é usada |
|---|---|---|
| `DT_SIN_PRI` | Data de início dos sintomas | Obrigatória. Define o mês e o ano do caso; registros sem data são excluídos. A série mensal usa o mês do início dos sintomas. |
| `SG_UF` (alternativas: `SG_UF_NOT`, `SG_UF_INTE`) | UF (código IBGE ou sigla) | Usa-se a primeira coluna disponível nessa ordem; a UF é normalizada e agrupada nas cinco macrorregiões (Norte, Nordeste, Centro-Oeste, Sudeste, Sul). Registros sem UF válida são excluídos da série regional. |
| `DT_NASC` | Data de nascimento | Idade em dias = `DT_SIN_PRI − DT_NASC` (usada quando ≥ 0). |
| `TP_IDADE` e `NU_IDADE_N` | Tipo da idade (1 = dias, 2 = meses, 3 = anos) e valor | Idade de reserva quando a data de nascimento está ausente ou inconsistente. A base de 2013–2018 não traz `TP_IDADE`, então depende de `DT_NASC`. |
| `PCR_VSR` (SIVEP-Gripe) / `RES_VSR` (SINAN) | RT-PCR positivo para VSR | Pré-filtro de elegibilidade (ver acima). |

**Faixa etária analisada.** `FAIXA_ALVO = "<=6 meses"`: crianças de até 6 meses, isto é, idade em meses (dias / 30,4375) entre 0 e 6, inclusive.

**Séries e variáveis derivadas.** Os casos são contados por macrorregião e mês (jan/2013 a dez/2025); o total do Brasil é a soma das regiões. Os atributos do modelo (função `criar_atributos`, seção 0) são: tendência, ano, mês, seno e cosseno sazonais, indicadores de período (pré-pandemia, transição pandêmica, a partir de 2023), tempo desde 2023, defasagens 1, 2, 3 e 12, médias móveis de 3 e 6 meses, diferença mensal e razão entre as defasagens 1 e 12.

**Principais colunas das tabelas de saída** (em `resultados/tabelas/`): `casos_forecast` (contrafactual sem vacina), `casos_forecast_li95` / `casos_forecast_ls95` (limites do intervalo de previsão de 95%), `casos_forecast_ajustado` (contrafactual ajustado por viés), `efetividade`, `cobertura` e `reducao_total` (cenário), `casos_pos_vacina`, `casos_evitados` e seus limites, `fator_vies` (por região), `casos_observados` e `diferenca_percentual` (comparação com 2026).

### Cobertura vacinal (seção 4)

O arquivo do painel traz, por UF e por mês de janeiro a junho de 2026, três colunas: cobertura mensal acumulada, número de doses (numerador: pessoas com ao menos uma dose registrada na Rede Nacional de Dados em Saúde, RNDS) e população (denominador: nascidos vivos do Sistema de Informações sobre Nascidos Vivos, SINASC). O arquivo soma doses e população das UFs de cada região e recalcula a cobertura como doses ÷ população, limitada a 100%. Essa cobertura é combinada com a efetividade de 60%, 70% e 75% para obter a redução esperada de janeiro a maio.

## Figuras do manuscrito (seção 5)

As figuras replotam saídas das seções 1–4, sem nova análise. Os arquivos saem em `resultados/graficos/`. Os dados observados de 2026 vêm da extração de 27/07/2026 (a mesma das tabelas).

| Figura no manuscrito | Arquivo gerado |
|---|---|
| Figure 1 (validação e contrafactual) | `figura1_validacao_lancet_en.png` / `.pdf` |
| Figure 2 (cenários por região) | `figura_cenarios_eixoY_fixo_en.png` |
| Figure 3 (observado vs contrafactual) | `figura3_observado_vs_contrafactual_lancet_en.png` / `.pdf` |
| Figure 4 (redução observada vs esperada) | `figura4_observado_vs_esperado_lancet_en.png` / `.pdf` |
| Supplementary Figure S1 (ajuste por viés) | `figuraS1_vies_lancet_en.png` / `.pdf` |

Observação: a Figure 2 ainda não foi regerada no padrão da revista (107 mm de largura, 300 dpi); as demais já estão nesse padrão.

## Software

R 4.5.3; `catboost` 1.2.10, `xgboost` 3.2.1.1, `lightgbm` 4.7.0. Semente aleatória: `set.seed(123)` (e `random_seed = 123` nos modelos).

## Licença e citação

A definir pelos autores antes da publicação do repositório.

---

## English summary

Reproduction code for the manuscript *"Potential impact of maternal vaccination against respiratory syncytial virus on severe acute respiratory infections in infants up to six months of age in Brazil, 2026"* (in preparation). The repository publishes a single file, `reproducao_artigo_vsr.R`.

**Scope.** The file contains the impact-estimation analyses only. Input files in `BD/` come from earlier data-preparation scripts (kept outside this repository) that apply the eligibility criteria to the public data, notably selection of RT-PCR-positive RSV cases; these prior filters are a prerequisite.

**Design.** Ecological time-series study using national surveillance data (SINAN Influenza Web 2009–2018; SIVEP-Gripe 2019 onwards; continuous monthly series 2013–December 2025) on RT-PCR-confirmed RSV-associated SARI hospitalizations in children up to six months of age.

**Methods.** Three gradient-boosting algorithms (XGBoost, LightGBM, CatBoost) were compared by expanding-window backtesting (2019 and 2022–2025). CatBoost, which ranked best in four of six units, was adopted as the single algorithm for all regions and projected 2026 hospitalizations under a no-vaccination counterfactual by region; the national point estimate is the sum of the five regions. An empirical 95% prediction interval was built from annual backtesting residuals. Deterministic scenarios combined vaccine effectiveness and coverage (60%×50%=30%, 70%×70%=49%, 75%×85%=63.75%). A bias-adjusted counterfactual was run as a sensitivity analysis. A non-causal, exploratory comparison used observed January–May 2026 hospitalizations and real vaccination coverage.

**How to run.** Open `impacto_vacina_vsr.Rproj` (paths use `here::here()`), place the four input files in `BD/`, and run `source("reproducao_artigo_vsr.R")`. Sections: 0 setup and functions, 1 monthly series, 2 backtesting and algorithm selection, 3 counterfactual, scenarios and bias-adjusted sensitivity analysis, 4 exploratory comparison with 2026 data and real coverage, 5 figures. Data and outputs (`BD/`, `resultados/`) are git-ignored. A full from-scratch run reproduced all 30 result tables of the analyses byte for byte.

**Data.** Public: OpenDATASUS (https://opendatasus.saude.gov.br) and the Ministry of Health vaccination coverage panel (https://infoms.saude.gov.br). Individual-level records are anonymized.

**Software.** R 4.5.3; catboost 1.2.10; xgboost 3.2.1.1; lightgbm 4.7.0; seed 123.
