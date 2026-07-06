# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                         BLOQUE 4 — Inferencia causal                         
#                           Código R · Fichas 11–14                            
# ==============================================================================

# ═══════════════════════════════════════════════════════════════════════════════
# IDENTIDAD VISUAL DEL LIBRO — paleta y tema (cargar una vez por sesión)
# ═══════════════════════════════════════════════════════════════════════════════
# Paleta coherente con el esquema del libro (navy/ochre + apoyo).
col_libro <- c(
  navy     = "#1A365D",  # primario: elementos principales, líneas, puntos
  ochre    = "#B7791F",  # secundario: contraste, segundo grupo
  teal     = "#2C7A7B",  # terciario categórico
  rojo     = "#C53030",  # énfasis / valores nulos / alertas
  verde    = "#2F855A",  # categórico adicional
  gris     = "#718096",  # referencia, líneas guía, texto secundario
  azul_med = "#2C5282"   # azul intermedio (continuidad con figuras previas)
)
# Escala categórica ordenada para grupos discretos
pal_libro <- unname(col_libro[c("navy", "ochre", "teal", "verde", "rojo", "gris")])

# Tema ggplot2 con acabado profesional (ready-to-publish)
theme_libro <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = base_size + 3,
                                               colour = col_libro["navy"]),
      plot.subtitle   = ggplot2::element_text(size = base_size, colour = col_libro["gris"]),
      plot.caption    = ggplot2::element_text(size = base_size - 3, colour = col_libro["gris"]),
      axis.title      = ggplot2::element_text(face = "bold", size = base_size),
      axis.text       = ggplot2::element_text(colour = "grey20"),
      legend.title    = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      plot.margin      = ggplot2::margin(10, 14, 10, 10)
    )
}
# Uso: p + scale_colour_manual(values = pal_libro) + theme_libro()
# Gráficas base (metafor/rms/rpart): usar col.* y palette con col_libro directamente.
# ═══════════════════════════════════════════════════════════════════════════════

#
# Fichas incluidas en este bloque:
#   Ficha 11 · Variables instrumentales (IV)
#   Ficha 12 · Aleatorización mendeliana (MR)
#   Ficha 13 · Puntuación de propensión (PS)
#   Ficha 14 · G-methods
#
# Contenido: 4 fichas · 22 fragmentos de código · 513 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v94.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 11 · Variables instrumentales (IV)
# ║  Método: AER::ivreg()/ivDiag · estimación IV, instrumentos débiles, sobreidentificación
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estimación clásica con AER::ivreg
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(AER) # ivreg() — gold standard
library(ivmodel) # IC Anderson-Rubin, LIML, k-class
library(sandwich) # vcovHC, vcovCL
library(lmtest) # coeftest
library(dplyr)
library(broom)
# Datos de ejemplo (CardKrueger): efecto de educación sobre salario
data(CollegeDistance, package = "AER")
dat <- CollegeDistance %>%
 filter(complete.cases(.))
# Sintaxis de AER::ivreg
# y ~ x1 + x2 + ... | z1 + z2 + ... + x2 + ... (después del |, instrumentos + exógenas)
# Modelo IV: efecto de 'education' sobre 'score'
# Z = distance to college (válido por argumento sustantivo: distancia al campus afecta
# la decisión educativa pero no las habilidades cognitivas si controlamos por factores
# socioeconómicos)
fit_iv <- ivreg(score ~ education + gender + ethnicity + income |
 distance + gender + ethnicity + income,
 data = dat)
summary(fit_iv, diagnostics = TRUE)
# La opción diagnostics=TRUE reporta:
# - Weak instruments: F de la 1.ª etapa (queremos F > 10)
# - Wu-Hausman: compara IV con OLS (si p < 0.05, IV preferible)
# - Sargan: test de validez (solo aplica si over-identified)

# ──────────────────────────────────────────────────────────────────────────────
# ■ SE robustos y cluster-robust
# ──────────────────────────────────────────────────────────────────────────────

# SE robustos a heterocedasticidad
coeftest(fit_iv, vcov = vcovHC(fit_iv, type = "HC1"))
# SE cluster-robust (si hay correlación dentro de centros, regiones, etc.)
coeftest(fit_iv, vcov = vcovCL(fit_iv, cluster = ~ region))
# IC robustos al instrumento débil con ivmodel
library(ivmodel)
iv_model <- ivmodel(Y = dat$score,
 D = dat$education,
 Z = dat$distance,
 X = model.matrix(~ gender + ethnicity + income, dat)[, -1])
summary(iv_model)
# Reporta: 2SLS, LIML, Fuller, Anderson-Rubin (robusto a inst. débil)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación OLS vs IV
# ──────────────────────────────────────────────────────────────────────────────

# Estimaciones lado a lado
fit_ols <- lm(score ~ education + gender + ethnicity + income, data = dat)
# Tabla comparativa
library(modelsummary)
modelsummary(list("OLS" = fit_ols, "IV (2SLS)" = fit_iv),
 vcov = c("classical", "robust"),
 statistic = "conf.int",
 stars = TRUE)
# Test de Hausman-Wu para comparar formalmente
# Disponible directamente en summary(fit_iv, diagnostics=TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnóstico del instrumento débil
# ──────────────────────────────────────────────────────────────────────────────

# Primera etapa explícita (para diagnóstico, NO para inferencia)
first_stage <- lm(education ~ distance + gender + ethnicity + income, data = dat)
summary(first_stage)
# F-test del coeficiente de 'distance'
library(car)
linearHypothesis(first_stage, "distance = 0")
# F > 10 → aceptable
# F < 10 → instrumento débil; usar Anderson-Rubin
# F efectivo de Olea-Pflueger (más estricto)
future::plan("sequential")  # ejecución secuencial — evita error de workers en paralelo
library(ivDiag)
ivDiag::ivDiag(data = dat,
 Y = "score",
 D = "education",
 Z = "distance",
 controls = c("gender", "ethnicity", "income"),
 cl = NULL)
# Reporta F_eff y umbrales actualizados

# ──────────────────────────────────────────────────────────────────────────────
# ■ Múltiples instrumentos (over-identified)
# ──────────────────────────────────────────────────────────────────────────────

# Cuando hay más instrumentos que endógenas, se puede testar la validez conjunta
# Ejemplo con 2 instrumentos
fit_iv_multi <- ivreg(score ~ education + gender + ethnicity + income |
 distance + tuition + gender + ethnicity + income,
 data = dat)
summary(fit_iv_multi, diagnostics = TRUE)
# El test de Sargan-Hansen aparece en diagnostics:
# H0: todos los instrumentos son válidos
# Si p < 0.05 → algún instrumento viola la exclusión

# ──────────────────────────────────────────────────────────────────────────────
# ■ IV para datos binarios y de supervivencia
# ──────────────────────────────────────────────────────────────────────────────

# CUIDADO: el 'forbidden regression' (logit + 2SLS naive) está SESGADO
# Soluciones correctas:
# 1) Bivariate probit (recpoisson + probit)
# Datos simulados: Y desenlace binario, X exposición endógena (binaria),
# W covariable exógena, Z instrumento
set.seed(42); n <- 500
W <- rnorm(n); Z <- rnorm(n)
X <- rbinom(n, 1, plogis(0.5 * Z + 0.3 * W))
Y <- rbinom(n, 1, plogis(-0.5 + 1.2 * X + 0.4 * W))
dat_gjrm <- data.frame(Y, X, W, Z)
library(GJRM)
biv_probit <- gjrm(list(Y ~ X + W, X ~ Z + W),
 data = dat_gjrm,
 model = "B",
 margins = c("probit", "probit"))
# 2) Two-stage residual inclusion (2SRI) - apropiado para no-lineales
# Etapa 1: residuos
first <- lm(X ~ Z + W, data = dat_gjrm)
dat_gjrm$resid_1 <- residuals(first)
# Etapa 2: GLM con los residuos como covariable
second <- glm(Y ~ X + resid_1 + W, data = dat_gjrm, family = binomial)
# El coeficiente de X es el efecto causal estimado
# 3) Para Cox / supervivencia: paquete ivtools (G-estimation)
library(ivtools)
library(survival)
library(survival)
# Simulación coherente: tiempos de fallo por inversión de la función de riesgo (Bender)
set.seed(42)
n_cox <- 500
W_cox <- rnorm(n_cox); Z_cox <- rnorm(n_cox)
X_cox <- rbinom(n_cox, 1, plogis(0.6 * Z_cox + 0.4 * W_cox))
dat_iv_cox <- data.frame(W = W_cox, Z = Z_cox, X = X_cox)
dat_iv_cox$time  <- -log(runif(n_cox)) / (0.1 * exp(0.7 * X_cox + 0.3 * W_cox))
dat_iv_cox$event <- rbinom(n_cox, 1, 0.8)
# G-estimation: fitZ.L (glm) + fitT.LZX (Cox con Z incluida)
fit_Z_L   <- glm(Z ~ W, data = dat_iv_cox)
fit_T_LZX <- coxph(Surv(time, event) ~ X + Z + W, data = dat_iv_cox)
# Simulación coherente: tiempos de fallo (método de Bender)
set.seed(42)
n_cox <- 500
W_cox <- rnorm(n_cox); Z_cox <- rnorm(n_cox)
X_cox <- rbinom(n_cox, 1, plogis(0.6 * Z_cox + 0.4 * W_cox))
dat_iv_cox <- data.frame(W = W_cox, Z = Z_cox, X = X_cox)
dat_iv_cox$time  <- -log(runif(n_cox)) / (0.1 * exp(0.7 * X_cox + 0.3 * W_cox))
dat_iv_cox$event <- rbinom(n_cox, 1, 0.8)
# G-estimation: fitZ.L (glm) + fitT.LZX (Cox con Z incluida)
fit_Z_L   <- glm(Z ~ W, data = dat_iv_cox)
fit_T_LZX <- coxph(Surv(time, event) ~ X + Z + W, data = dat_iv_cox)
fit_iv_cox <- ivcoxph(estmethod = "g", X = "X",
 fitZ.L   = fit_Z_L,
 fitT.LZX = fit_T_LZX,
 data     = dat_iv_cox)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 12 · Aleatorización mendeliana (MR)
# ║  Método: TwoSampleMR · IVW, MR-Egger, weighted median; MR-PRESSO; MVMR
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Paquete TwoSampleMR del MRC IEU
# ──────────────────────────────────────────────────────────────────────────────

# Instalación
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("rondolab/MR-PRESSO")
library(TwoSampleMR)
library(MRPRESSO)
library(MendelianRandomization)
library(ieugwasr) # acceso a OpenGWAS
# Token OpenGWAS: añadir OPENGWAS_JWT="eyJhb..." en ~/.Renviron y reiniciar R
Sys.getenv("OPENGWAS_JWT")  # devuelve "" si no está configurado
library(dplyr)
# ---------------- 1) Obtener SNPs instrumento para la exposición --------------
# Ejemplo: HDL-c como exposición, CHD como desenlace
# Buscar el ID del GWAS de HDL
ao <- available_outcomes()
head(ao %>% filter(grepl("HDL", trait, ignore.case = TRUE)),
 n = 5)
# IDs (ejemplos: ieu-b-109 = HDL, ieu-a-7 = CHD)
exp_dat <- extract_instruments(outcomes = "ieu-b-109",
 p1 = 5e-08, # umbral GWS
 clump = TRUE, # LD clumping
 r2 = 0.001,
 kb = 10000)
# Verificar nº de SNPs y F-statistic
nrow(exp_dat)
exp_dat <- exp_dat %>%
 mutate(F_stat = beta.exposure^2 / se.exposure^2)
summary(exp_dat$F_stat) # mediana > 10 ideal

# ──────────────────────────────────────────────────────────────────────────────
# ■ Extracción de efectos sobre el desenlace y harmonización
# ──────────────────────────────────────────────────────────────────────────────

# ---------------- 2) Extraer efectos de los SNPs sobre el desenlace ----------
out_dat <- extract_outcome_data(snps = exp_dat$SNP,
 outcomes = "ieu-a-7") # CHD
# ---------------- 3) Harmonización -------------------------------------------
# Asegurar que el alelo de efecto es el mismo en ambos GWAS
dat <- harmonise_data(exposure_dat = exp_dat,
 outcome_dat = out_dat,
 action = 2) # 2 = inferir strand basado en frecuencia alélica
# Verificar harmonización
table(dat$mr_keep) # FALSE indica SNPs problemáticos (paliándromos ambiguos)
dat <- dat %>% filter(mr_keep)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis principal (IVW) y métodos de sensibilidad
# ──────────────────────────────────────────────────────────────────────────────

# ---------------- 4) Análisis principal --------------------------------------
res <- mr(dat,
 method_list = c("mr_ivw", # principal
 "mr_egger_regression", # sensibilidad
 "mr_weighted_median",
 "mr_weighted_mode",
 "mr_simple_mode"))
# Convertir a OR si el desenlace es binario (HDL → CHD)
res_or <- generate_odds_ratios(res)
print(res_or %>% dplyr::select(method, nsnp, b, or, or_lci95, or_uci95, pval))
# ---------------- 5) Diagnósticos --------------------------------------------
# Test de intercepto de MR-Egger (pleiotropía direccional)
egger_int <- mr_pleiotropy_test(dat)
print(egger_int)
# Si pval < 0.05 → evidencia de pleiotropía direccional
# Heterogeneidad entre SNPs
het <- mr_heterogeneity(dat)
print(het)
# Q_IVW alto + I² > 50% → heterogeneidad relevante
# Leave-one-out: descartar un SNP cada vez
loo <- mr_leaveoneout(dat)
mr_leaveoneout_plot(loo) +
 ggtitle("Análisis leave-one-out — sensibilidad a SNPs individuales")

# ──────────────────────────────────────────────────────────────────────────────
# ■ MR-PRESSO para detectar outliers pleiotrópicos
# ──────────────────────────────────────────────────────────────────────────────

# ---------------- 6) MR-PRESSO -----------------------------------------------
presso <- mr_presso(
 BetaOutcome = "beta.outcome",
 BetaExposure = "beta.exposure",
 SdOutcome = "se.outcome",
 SdExposure = "se.exposure",
 OUTLIERtest = TRUE,
 DISTORTIONtest = TRUE,
 data = dat,
 NbDistribution = 5000,
 SignifThreshold = 0.05
)
# El test global indica si HAY pleiotropía;
# OutlierTest identifica SNPs específicos;
# DistortionTest evalúa si eliminarlos cambia las conclusiones
print(presso$`Main MR results`)
print(presso$`MR-PRESSO results`$`Global Test`)
print(presso$`MR-PRESSO results`$`Outlier Test`)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización: scatter plot y forest plot
# ──────────────────────────────────────────────────────────────────────────────

# ---------------- 7) Visualizaciones -----------------------------------------
# Scatter plot con las distintas estimaciones
mr_scatter_plot(res, dat) +
 ggtitle("Estimaciones de aleatorización mendeliana por método")
# Forest plot de ratios de Wald individuales por SNP
mr_forest_plot(mr_singlesnp(dat)) +
 ggtitle("Ratios de Wald por SNP individual")
# Funnel plot para detectar asimetría (pleiotropía direccional)
mr_funnel_plot(mr_singlesnp(dat)) +
 ggtitle("Funnel plot — asimetría y pleiotropía direccional")
# Todos los gráficos en un panel
plots <- mr_singlesnp(dat) %>% mr_forest_plot()
library(ggplot2)  # requerido por ggsave()
ggsave("mr_forest.png", plots[[1]], width = 8, height = 10)

# · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·

# PSEUDOCÓDIGO ILUSTRATIVO — no ejecutable sin datos propios del usuario.
# Requiere columnas: snp1..snpN, weight1..weightN, HDL, CHD, age, sex, pc1-pc4.
# Sustituir por los nombres reales del dataset propio.
library(AER) # ivreg() — 2SLS clásico (Ficha 11)
library(dplyr)
# Datos individuales con G, X, Y y covariables
# G puede ser: un SNP único, dosis alélica, o un PRS (Polygenic Risk Score)
# Construcción de un PRS instrumento
# (sumas ponderadas de alelos según efectos de un GWAS externo)
dat$PRS_hdl <- with(dat, snp1 * weight1 + snp2 * weight2 + snp3 * weight3 + ...)
# F-statistic de la primera etapa
fs <- lm(HDL ~ PRS_hdl + age + sex + pc1 + pc2 + pc3 + pc4, data = dat)
summary(fs)$fstatistic # F > 10 = aceptable
# MR como IV con 2SLS
fit_mr <- ivreg(
 CHD ~ HDL + age + sex + pc1 + pc2 + pc3 + pc4 |
 PRS_hdl + age + sex + pc1 + pc2 + pc3 + pc4,
 data = dat
)
summary(fit_mr, diagnostics = TRUE)
# SE robustos (recomendable)
library(sandwich); library(lmtest)
coeftest(fit_mr, vcov = vcovHC(fit_mr, type = "HC1"))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Multivariable MR (MVMR)
# ──────────────────────────────────────────────────────────────────────────────

library(MVMR)
# Datos: tres exposiciones correlacionadas (HDL, LDL, TG) sobre CHD
# Datos simulados: K SNPs instrumentando 3 exposiciones (HDL, LDL, TG) sobre CHD
set.seed(2024); K <- 100
dat_mvmr <- data.frame(
 SNP      = paste0("rs", 1:K),
 beta_HDL = rnorm(K, 0, 0.05), se_HDL = runif(K, 0.01, 0.03),
 beta_LDL = rnorm(K, 0, 0.05), se_LDL = runif(K, 0.01, 0.03),
 beta_TG  = rnorm(K, 0, 0.05), se_TG  = runif(K, 0.01, 0.03),
 beta_CHD = rnorm(K, 0, 0.10), se_CHD = runif(K, 0.02, 0.05)
)
mvmr_input <- format_mvmr(
 BXGs = cbind(dat_mvmr$beta_HDL, dat_mvmr$beta_LDL, dat_mvmr$beta_TG),
 BYG = dat_mvmr$beta_CHD,
 seBXGs = cbind(dat_mvmr$se_HDL, dat_mvmr$se_LDL, dat_mvmr$se_TG),
 seBYG = dat_mvmr$se_CHD,
 RSID = dat_mvmr$SNP
)
# Test de instrumentos débiles condicionales
strength_mvmr(mvmr_input)
# Estimación MVMR-IVW
mv_result <- ivw_mvmr(mvmr_input)
print(mv_result)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 13 · Puntuación de propensión (PS)
# ║  Método: MatchIt/IPTW/estratificación · cuatro estrategias de ajuste por PS y balance
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y estimación del PS
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(MatchIt) # matching y estimación
library(WeightIt) # IPTW
library(cobalt) # love plots y balance
library(survey) # análisis ponderado
library(survival) # Cox y Kaplan-Meier
library(dplyr)
library(ggplot2)
# Datos: estudio de efectividad de un fármaco (ejemplo simulado)
set.seed(2024)
N <- 5000
dat <- data.frame(
 id = 1:N,
 tratado = rbinom(N, 1, 0.3), # 30% tratados (no aleatorio)
 edad = NA, sexo = NA, IMC = NA, HTA = NA, DM = NA,
 fum = NA, alcohol = NA, educacion = NA
)
# Generar covariables con desbalance (más severo en los tratados)
dat$edad <- ifelse(dat$tratado == 1,
 rnorm(N, 62, 8), rnorm(N, 55, 10))
dat$IMC <- ifelse(dat$tratado == 1,
 rnorm(N, 28, 3.5), rnorm(N, 26, 4))
dat$HTA <- ifelse(dat$tratado == 1,
 rbinom(N, 1, 0.55), rbinom(N, 1, 0.30))
dat$DM <- ifelse(dat$tratado == 1,
 rbinom(N, 1, 0.30), rbinom(N, 1, 0.15))
dat$fum <- rbinom(N, 1, 0.30)
dat$sexo <- rbinom(N, 1, 0.50)
dat$alcohol <- rnorm(N, 15, 12)
dat$educacion <- sample(c("baja", "media", "alta"), N, replace = TRUE)
# Generar tiempo y evento con HR_verdadero = 0.65 (tratamiento protector)
hazard <- 0.02 * exp(-0.43 * dat$tratado +
 0.03 * (dat$edad - 55) +
 0.05 * dat$HTA + 0.06 * dat$DM)
dat$tiempo <- rexp(N, hazard)
dat$tiempo_cens <- runif(N, 5, 15)
dat$evento <- as.numeric(dat$tiempo <= dat$tiempo_cens)
dat$tiempo <- pmin(dat$tiempo, dat$tiempo_cens)
# Definición ÚNICA del conjunto de covariables del PS
ps_vars <- c("edad", "sexo", "IMC", "HTA", "DM", "fum",
 "alcohol", "educacion")
formula_ps <- as.formula(
 paste("tratado ~", paste(ps_vars, collapse = " + "))
)
# Estimación del PS por regresión logística
fit_ps <- glm(formula_ps, data = dat, family = binomial)
dat$ps <- predict(fit_ps, type = "response")
# Verificación de solapamiento (common support)
ggplot(dat, aes(x = ps, fill = factor(tratado))) +
 geom_density(alpha = 0.5) +
 scale_fill_manual(values = c("0" = unname(col_libro["navy"]),
 "1" = unname(col_libro["ochre"])),
 labels = c("Control", "Tratado")) +
 labs(title    = "Solapamiento del propensity score entre grupos",
 subtitle = "Distribución de la probabilidad de tratamiento (supuesto de positividad)",
 x        = "Propensity score",
 y        = "Densidad",
 fill     = "Grupo",
 caption  = "El solapamiento adecuado es condición necesaria para la comparabilidad") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estrategia 1: Matching por PS
# ──────────────────────────────────────────────────────────────────────────────

# Matching 1:1 con nearest neighbor (NN) con caliper de 0.2 SD
m_out <- matchit(formula_ps, data = dat,
 method = "nearest",
 distance = "logit",
 ratio = 1, # 1 control por tratado
 caliper = 0.2, # tolerancia en SD del logit del PS
 replace = FALSE)
summary(m_out, standardize = TRUE) # balance pre/post
# Extracción del dataset emparejado
dat_matched <- match.data(m_out)
# Love plot
library(cobalt)
love.plot(m_out, threshold = 0.10, abs = TRUE,
 var.order = "unadjusted",
 colors = c(unname(col_libro["rojo"]), unname(col_libro["verde"])),
 title = "Balance de covariables antes y después del emparejamiento")
# Análisis sobre los matched (estima ATT)
library(survival)
fit_cox_matched <- coxph(Surv(tiempo, evento) ~ tratado,
 data = dat_matched,
 cluster = subclass, # cluster por par
 robust = TRUE)
summary(fit_cox_matched)
# El HR estimado es el ATT

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estrategia 2: IPTW
# ──────────────────────────────────────────────────────────────────────────────

# Estimación del peso IPTW para ATE
dat$w_ate <- with(dat, ifelse(tratado == 1, 1/ps, 1/(1 - ps)))
# Peso IPTW para ATT (controles ponderados, tratados sin peso)
dat$w_att <- with(dat, ifelse(tratado == 1, 1, ps / (1 - ps)))
# Diagnóstico de pesos: distribución y truncado
summary(dat$w_ate)
quantile(dat$w_ate, c(0.01, 0.99))
# Truncado al 1.º y 99.º percentil si hay extremos
dat$w_ate_t <- pmin(pmax(dat$w_ate,
 quantile(dat$w_ate, 0.01)),
 quantile(dat$w_ate, 0.99))
# Alternativa: WeightIt con stabilización
library(WeightIt)
w_obj <- weightit(formula_ps, data = dat,
 method = "ps",
 estimand = "ATE",
 stabilize = TRUE)
dat$w_stab <- w_obj$weights
# Balance tras ponderación
bal.tab(w_obj, threshold = 0.10)
# Análisis ponderado para ATE
library(survey)
design_ate <- svydesign(ids = ~ id, weights = ~ w_stab, data = dat)
fit_cox_ipt <- svycoxph(Surv(tiempo, evento) ~ tratado, design = design_ate)
summary(fit_cox_ipt)
# El HR es el ATE bajo IPTW

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estrategia 3: Estratificación por quintiles
# ──────────────────────────────────────────────────────────────────────────────

# Quintiles del PS
dat$q_ps <- cut(dat$ps,
 breaks = quantile(dat$ps, probs = seq(0, 1, 0.2)),
 include.lowest = TRUE, labels = 1:5)
table(dat$q_ps, dat$tratado)
# Análisis dentro de cada estrato (Cox)
fit_cox_strat <- coxph(Surv(tiempo, evento) ~ tratado + strata(q_ps),
 data = dat)
summary(fit_cox_strat)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estrategia 4: Ajuste por regresión con PS
# ──────────────────────────────────────────────────────────────────────────────

# Inclusión del PS (mejor un spline) como covariable
library(splines)
fit_cox_reg <- coxph(Surv(tiempo, evento) ~ tratado + ns(ps, df = 4),
 data = dat)
summary(fit_cox_reg)
# El HR de 'tratado' es el efecto ajustado por PS
# El spline modela la relación PS-Y de forma flexible

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación del balance — diagnóstico universal
# ──────────────────────────────────────────────────────────────────────────────

# Balance antes vs después con cobalt
library(cobalt)
# Para matching
bal.tab(m_out, threshold = 0.10, m.threshold = 0.10)
# Para IPTW
bal.tab(w_obj, threshold = 0.10)
# Love plot combinado (varios métodos)
love.plot(formula_ps,
 data      = dat,
 estimand  = "ATT",
 weights   = list("Matching" = get.w(m_out), "IPTW" = get.w(w_obj)),
 threshold = 0.10,
 abs       = TRUE,
 var.order = "unadjusted",
 sample.names = c("Sin ajuste", "Matching", "IPTW"),
 colors       = unname(col_libro[c("gris", "navy", "ochre")]),
 title        = "Comparación del balance: sin ajuste vs. matching vs. IPTW",
 stars        = "std")   # distingue dif. estandarizadas (continuas) de brutas (binarias)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 14 · G-methods
# ║  Método: gfoRmula/g-formula · g-computación, IPTW marginal structural models, g-estimación
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Implementación en R
# ──────────────────────────────────────────────────────────────────────────────

library(gfoRmula)
library(data.table)
library(dplyr)
# Datos longitudinales en formato long (un registro por sujeto-tiempo)
# Variables: id, tiempo (0, 1, 2, ...), A (tratamiento), L (covariable),
# Y (desenlace al final)
# Definir el conjunto de covariables tiempo-dependientes
covar_set <- c("L", "edad_basal", "sexo")
# Ajustar la g-formula
# Dataset longitudinal simulado: 500 sujetos, 5 periodos (formato long)
set.seed(2024)
n_id <- 500; n_t <- 5
dat_gf <- data.table(
 id         = rep(1:n_id, each = n_t),
 tiempo     = rep(0:(n_t - 1), times = n_id),
 edad_basal = rep(rnorm(n_id, 55, 10), each = n_t),
 sexo       = rep(rbinom(n_id, 1, 0.5), each = n_t)
)
dat_gf[, L := rnorm(.N, 0.3 * tiempo + 0.1 * edad_basal, 1)]
dat_gf[, A := rbinom(.N, 1, plogis(-1 + 0.5 * L + 0.2 * tiempo))]
dat_gf[, Y := ifelse(tiempo == n_t - 1,
                     rnorm(.N, 2 * A + 0.5 * L + 0.01 * edad_basal, 1),
                     NA_real_)]
# Ajustar la g-formula
fit_gf <- gformula_continuous_eof(
 obs_data     = dat_gf,
 id           = "id",
 time_name    = "tiempo",
 basecovs     = c("edad_basal", "sexo"),  # covariables fijas en el tiempo
 covnames     = c("L", "A"),
 covtypes     = c("normal", "binary"),
 histvars     = list(c("L", "A")),
 histories    = c(lagged),           # infiere lags del prefijo lag1_ en covparams
 outcome_name = "Y",
 ymodel       = Y ~ A + L + edad_basal + sexo,
 covparams = list(covmodels = c(
 L ~ lag1_L + lag1_A + edad_basal,
 A ~ lag1_L + lag1_A + L + edad_basal
 )),
 intvars = list("A", "A"),
 interventions = list(list(c(static, rep(0, 5))), # nunca tratado
 list(c(static, rep(1, 5)))), # siempre tratado
 int_descript = c("Nunca tratado", "Siempre tratado"),
 ref_int = 1,
 nsamples = 200, # Monte Carlo
 nsimul = 5000,
 seed = 2024
)
print(fit_gf$result)
# Ê[Y(0)] vs Ê[Y(1)] y contraste con IC 95%

# · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·

library(ipw) # cálculo de IPTW longitudinales
library(survey) # análisis ponderado
library(survival)
library(dplyr)
# Reutilizar el dataset longitudinal simulado de la g-formula
dat_msm <- copy(dat_gf)          # data.table: copy() evita modificar dat_gf
setDT(dat_msm)
setorder(dat_msm, id, tiempo)
dat_msm[, lag1_A := shift(A, 1, type = "lag"), by = id]
dat_msm[, lag1_L := shift(L, 1, type = "lag"), by = id]
dat_msm <- dat_msm[!is.na(lag1_A)]  # eliminar tiempo 0 (sin lag previo)
# Estimación manual de pesos IPTW estabilizados (ipwtm sustituido por glm directo)
# Modelo numerador: P(A_t | V)
fit_num <- glm(A ~ edad_basal + sexo + factor(tiempo),
              data = dat_msm, family = binomial())
# Modelo denominador: P(A_t | L_t, A_{t-1}, V)
fit_den <- glm(A ~ L + lag1_A + edad_basal + sexo + factor(tiempo),
              data = dat_msm, family = binomial())
eps <- 1e-6
p_num <- pmin(pmax(predict(fit_num, type = "response"), eps), 1 - eps)
p_den <- pmin(pmax(predict(fit_den, type = "response"), eps), 1 - eps)
dat_msm <- dat_msm %>%
  mutate(prob_num = ifelse(A == 1, p_num, 1 - p_num),
         prob_den = ifelse(A == 1, p_den, 1 - p_den),
         w_inst   = prob_num / prob_den) %>%
  arrange(id, tiempo) %>% group_by(id) %>%
  mutate(sw = cumprod(w_inst)) %>% ungroup()
# Diagnóstico de pesos
summary(dat_msm$sw)
quantile(dat_msm$sw, c(0.01, 0.25, 0.5, 0.75, 0.99))
# Idealmente: media ≈ 1, máximo < 10
# Truncado
dat_msm$sw_t <- pmin(pmax(dat_msm$sw, quantile(dat_msm$sw, 0.01)),
 quantile(dat_msm$sw, 0.99))
# MSM con GLM ponderado sobre el último periodo
dat_msm_eof <- subset(dat_msm, tiempo == max(dat_msm$tiempo))
setDT(dat_msm_eof)
dat_msm_eof[, evento := as.integer(Y > mean(Y, na.rm = TRUE))]
design_msm <- svydesign(ids = ~ id, weights = ~ sw_t, data = dat_msm_eof)
fit_msm <- svyglm(evento ~ A + edad_basal + sexo,
  design = design_msm,
  family = quasibinomial())
summary(fit_msm)
# OR causal de A bajo el MSM
exp(coef(fit_msm)["A"])
exp(confint(fit_msm)["A", ])

# · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·

# gesttools fue archivado en CRAN (2025); se usa DTRreg como alternativa mantenida
library(DTRreg)
# DTRreg: g-estimation de SNMM en formato wide (una fila por sujeto)
# Añadir desenlace continuo Y al final del seguimiento en dat_msm
dat_msm <- dat_msm %>% group_by(id) %>%
  mutate(Y = ifelse(tiempo == max(tiempo),
    -0.5 * mean(A) + 0.3 * last(L) + 0.02 * first(edad_basal) + rnorm(1, 0, 1),
    NA_real_)) %>% ungroup()
dat_ge <- dat_msm %>% filter(!is.na(Y)) %>%
  dplyr::select(id, Y, A, L, edad_basal, sexo) %>% as.data.frame()
fit_dtr <- DTRreg(
  outcome   = dat_ge$Y,
  blip.mod  = list(~ L + edad_basal),       # efecto del tto depende de L
  treat.mod = list(A ~ L + edad_basal),     # modelo del tratamiento
  tf.mod    = list(~ L + edad_basal + sexo),
  data      = dat_ge,
  method    = "gest",                        # g-estimation de SNMM
  var.estim = "none"
)
summary(fit_dtr)
# psi: efecto causal del tratamiento ajustado por L y edad_basal
