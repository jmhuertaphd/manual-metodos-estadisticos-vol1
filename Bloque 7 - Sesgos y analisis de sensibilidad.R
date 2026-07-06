# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                 BLOQUE 7 — Sesgos y análisis de sensibilidad                 
#                           Código R · Fichas 22–24                            
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
#   Ficha 22 · E-value para confusión no medida
#   Ficha 23 · Errores de clasificación
#   Ficha 24 · Imputación múltiple (MICE)
#
# Contenido: 3 fichas · 18 fragmentos de código · 385 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v92.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 22 · E-value para confusión no medida
# ║  Método: EValue::evalues.RR() + funciones cerradas e_value_or/hr/md, bias_factor y sensemakr
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Paquete EValue
# ──────────────────────────────────────────────────────────────────────────────

library(EValue)
library(dplyr)
# E-value para un RR puntual con su IC
evalues.RR(est = 1.80, lo = 1.50, hi = 2.20)
# Reporta:
# - E-value para la estimación puntual
# - E-value para el límite del IC más cercano al nulo
# Para OR con desenlace raro
# Funciones cerradas (VanderWeele-Ding 2017; VanderWeele-Vansteelandt para OR no raro;
# Vansteelandt 2008 para HR). EValue::evalues.OR/HR/MD() fallan con lava cargado.
ev_vd <- function(x) if (x <= 1) ev_vd(1 / x) else x + sqrt(x * (x - 1))
e_value_or <- function(est, lo, hi, rare = TRUE) {
  to_rr <- function(x) if (rare) x else 2 * x / (1 + x)
  rr <- to_rr(est)
  rr_ic <- if (lo <= 1 && hi >= 1) 1 else to_rr(if (est > 1) lo else hi)
  c(point = ev_vd(rr), ic = ev_vd(rr_ic))
}
# HR: aproximación de Vansteelandt (2008); rare distingue la transformación HR→RR
e_value_hr <- function(est, lo, hi, rare = FALSE) {
  to_rr <- function(x) if (rare) x else (1 - 0.5^sqrt(x)) / (1 - 0.5^sqrt(1 / x))
  rr <- if (est < 1) 1 / to_rr(1 / est) else to_rr(est)
  ic_raw <- if (lo <= 1 && hi >= 1) 1 else (if (est > 1) lo else hi)
  rr_ic <- if (ic_raw == 1) 1 else (if (ic_raw < 1) 1 / to_rr(1 / ic_raw) else to_rr(ic_raw))
  c(point = ev_vd(rr), ic = ev_vd(rr_ic))
}
# MD: convierte d de Cohen a RR vía la aproximación de VanderWeele (2017),
# asumiendo distribución aproximadamente normal: d = est / sd_pooled
e_value_md <- function(est, se, sd_pooled = 1) {
  d <- est / sd_pooled; se_d <- se / sd_pooled
  rr_approx <- exp(0.91 * d) # VanderWeele (2017), Stat Med
  lo <- exp(0.91 * (d - 1.96 * se_d)); hi <- exp(0.91 * (d + 1.96 * se_d))
  rr_lo <- min(lo, hi); rr_hi <- max(lo, hi)
  ic_rr <- if (rr_lo <= 1 && rr_hi >= 1) 1 else (if (rr_approx > 1) rr_lo else rr_hi)
  c(point = ev_vd(rr_approx), ic = ev_vd(ic_rr))
}
e_value_or(est = 2.50, lo = 1.80, hi = 3.50, rare = TRUE)  # EValue::evalues.OR() falla con lava cargado; usar función cerrada
# rare = TRUE: trata OR como RR directamente
# rare = FALSE: aplica la aproximación VanderWeele-Vansteelandt
# Para HR (Vansteelandt 2008)
e_value_hr(est = 0.65, lo = 0.50, hi = 0.85, rare = FALSE)  # EValue::evalues.HR() falla con lava; usar función cerrada
# Para diferencia de medias estandarizada (Cohen's d → RR)
e_value_md(est = 0.45, se = 0.10)  # EValue::evalues.MD() falla con lava; usar función cerrada

# ──────────────────────────────────────────────────────────────────────────────
# ■ Bias factor y confounded RR
# ──────────────────────────────────────────────────────────────────────────────

# EValue::bias_factor() fue eliminada del paquete; definición local
# (Lin, Psaty & Kronmal 1998 / VanderWeele-Ding 2017):
bias_factor <- function(RR_EU, RR_UD) (RR_EU * RR_UD) / (RR_EU + RR_UD - 1)
bias_factor(RR_EU = 2.5, RR_UD = 2.0)
# Retorna: BF = 1.43
# RR observado vs RR corregido
observed_rr <- 1.85
corrected_rr <- observed_rr / bias_factor(RR_EU = 2.5, RR_UD = 2.0)
cat("RR observado:", observed_rr, "\n")
cat("RR corregido:", round(corrected_rr, 2), "\n")
# Array approach: matriz de RR corregidos para distintos parámetros
rau_vals <- c(1.5, 2.0, 2.5, 3.0, 4.0)
ruy_vals <- c(1.5, 2.0, 2.5, 3.0, 4.0)
bias_matrix <- outer(rau_vals, ruy_vals,
 FUN = function(a, u) {
 bf <- (a * u) / (a + u - 1)
 observed_rr / bf
 })
rownames(bias_matrix) <- paste0("RR_AU=", rau_vals)
colnames(bias_matrix) <- paste0("RR_UY=", ruy_vals)
print(round(bias_matrix, 2))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización del E-value
# ──────────────────────────────────────────────────────────────────────────────

library(EValue)
library(ggplot2)
# Curva de iso-confusión para un RR observado
# bias_plot() tiene API inestable entre versiones de EValue;
# se implementa la curva de iso-confusión con ggplot2 (ver código a continuación)
# Muestra las combinaciones (RR_AU, RR_UY) suficientes para anular el efecto
# Bubble plot personalizado con varios estudios
studies <- data.frame(
 estudio = c("A", "B", "C", "D"),
 RR = c(1.45, 2.10, 1.30, 3.50),
 RR_lo = c(1.15, 1.60, 1.05, 2.50),
 RR_hi = c(1.83, 2.75, 1.61, 4.90)
) %>%
 mutate(
 E_RR = RR + sqrt(RR * (RR - 1)),
 E_IC = RR_lo + sqrt(RR_lo * (RR_lo - 1))
 )
ggplot(studies, aes(x = E_RR, y = E_IC, label = estudio)) +
 geom_point(size = 4, color = col_libro["navy"]) +
 geom_text(nudge_y = 0.1, fontface = "bold", color = col_libro["navy"]) +
 geom_abline(slope = 1, intercept = 0, linetype = "dashed",
 color = col_libro["gris"]) +
 geom_hline(yintercept = 1.5, color = col_libro["rojo"], linetype = "dotted") +
 labs(title    = "Robustez frente a confusión no medida por estudio",
 subtitle = "E-value de la estimación puntual frente al del límite del IC",
 x        = "E-value para el RR puntual",
 y        = "E-value para el límite del IC",
 caption  = "Línea de puntos: umbral de referencia E-value = 1.5") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Aplicación en meta-análisis
# ──────────────────────────────────────────────────────────────────────────────

library(metafor)
library(EValue)
# Ejemplo: meta-análisis con su estimación combinada
# Datos simulados: 8 estudios con log(RR) y su varianza (escala log,
# coherente con rma() y la posterior exponenciación)
set.seed(2024)
n_studies <- 8
dat_meta <- data.frame(
  estudio = paste0("Estudio_", 1:n_studies),
  yi = log(c(1.42, 1.68, 1.35, 1.91, 1.55, 1.73, 1.48, 1.62)), # log(RR) por estudio
  vi = c(0.015, 0.022, 0.018, 0.028, 0.012, 0.020, 0.016, 0.019) # varianza de yi
)
fit <- rma(yi, vi, data = dat_meta, method = "REML")
rr_combined <- exp(coef(fit))
rr_lo <- exp(fit$ci.lb)
rr_hi <- exp(fit$ci.ub)
# E-value para la síntesis
evalues.RR(est = rr_combined, lo = rr_lo, hi = rr_hi)
# Reportar en la publicación junto con τ², I² y la fuerza de
# los confusores medidos para juzgar la plausibilidad

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis de sensibilidad alternativo: sensemakr
# ──────────────────────────────────────────────────────────────────────────────

library(sensemakr)
# Para regresiones lineales con confusor continuo
# (Cinelli & Hazlett 2020)
# Datos simulados: X1 y X2 confunden A-Y (sirven de benchmark);
# X3 es ruido sin relación causal real, incluido como contraste
set.seed(2024)
n <- 500
X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
A <- 0.5 * X1 + 0.3 * X2 + rnorm(n)
Y <- 1.2 * A + 0.8 * X1 + 0.6 * X2 + 0.2 * X3 + rnorm(n)
dat <- data.frame(Y, A, X1, X2, X3)
fit <- lm(Y ~ A + X1 + X2 + X3, data = dat)
sens <- sensemakr(fit,
 treatment = "A",
 benchmark_covariates = c("X1", "X2"),
 kd = 1:3)
# kd: múltiplos de la fuerza de los benchmark covariates
summary(sens)
# Reporta:
# - Robustness value (RV): mínimo R² parcial del confusor para anular el efecto
# - Comparación con la fuerza de los confusores medidos
# - Contour plot con isolíneas de β estimado
plot(sens)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 23 · Errores de clasificación
# ║  Método: episensr/simex/mecor · sensibilidad a mala clasificación
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y exploración inicial
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(episensr) # Lash, Fox, MacLehose — gold standard
library(simex) # SIMEX para variables continuas
library(mecor) # Misclassification correction
library(dplyr)
library(ggplot2)
# Datos: estudio caso-control con exposición auto-reportada
# Tabla 2×2 observada
obs_2x2 <- matrix(c(150, 75, 50, 100), nrow = 2, byrow = TRUE,
 dimnames = list(c("Expuesto", "No expuesto"),
 c("Caso", "Control")))
print(obs_2x2)
# OR observado
or_obs <- (obs_2x2[1,1] * obs_2x2[2,2]) / (obs_2x2[1,2] * obs_2x2[2,1])
cat("OR observado:", round(or_obs, 2), "\n")
# OR observado: 4.00
# Si se conocen Se y Sp del cuestionario por validación interna
se_a <- 0.85 # sensibilidad para exposición
sp_a <- 0.90 # especificidad para exposición

# ──────────────────────────────────────────────────────────────────────────────
# ■ Matrix method (corrección puntual)
# ──────────────────────────────────────────────────────────────────────────────

library(episensr)
# Corrección no diferencial en la exposición
result_misclass <- misclassification(
 obs_2x2,
 type = "exposure", # error en la exposición (no diferencial)
 bias_parms = c(se_a, se_a, sp_a, sp_a)
 # bias_parms: Se_casos, Se_controles, Sp_casos, Sp_controles
 # Para NO diferencial: idénticos en ambos grupos
)
summary(result_misclass)
# Reporta:
# - OR observado
# - OR corregido
# - Tabla 2×2 corregida
# Para corrección DIFERENCIAL: distintos parámetros por grupo
result_diff <- misclassification(
 obs_2x2,
 type = "exposure",
 bias_parms = c(0.90, 0.80, 0.92, 0.88)
 # Casos sobre-reportan: Se mayor, Sp similar
)
summary(result_diff)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Probabilistic Bias Analysis (PBA)
# ──────────────────────────────────────────────────────────────────────────────

# PBA con priors sobre Se y Sp
result_pba <- probsens(
 obs_2x2,
 type = "exposure",
 reps = 10000, # 10 000 simulaciones
 seca = list("trapezoidal", c(0.75, 0.85, 0.90, 0.95)),
 # distribución trapezoidal para Se
 spca = list("trapezoidal", c(0.85, 0.90, 0.93, 0.97)),
 # distribución trapezoidal para Sp
 # print = TRUE eliminado: no es argumento de probsens()
)
summary(result_pba)
# Reporta:
# - OR corregido (mediana y percentiles)
# - IC simulado 95%
# - Histograma de OR corregidos
# Visualización
plot(result_pba, "or", main = "Análisis probabilístico de sesgo — distribución del OR corregido")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Regression calibration para variables continuas
# ──────────────────────────────────────────────────────────────────────────────

library(mecor)
# Ejemplo: ingesta de sodio (X) con error de medida (X*)
# Validación interna: 100 sujetos con ambas mediciones
# Modelo principal: PA = β₀ + β·sodio + γ·edad + ε
# Datos
dat <- data.frame(
 presion = rnorm(1000, 130, 15),
 sodio_obs = rnorm(1000, 3500, 800), # con error
 sodio_true = NA, # solo para sub-muestra
 edad = rnorm(1000, 55, 12)
)
# Sub-muestra con validación
val_idx <- sample(1:1000, 100)
dat$sodio_true[val_idx] <- dat$sodio_obs[val_idx] +
 rnorm(100, 0, 200) # gold standard con menos error
# Modelo principal SIN corrección (sesgado)
fit_naive <- lm(presion ~ sodio_obs + edad, data = dat)
coef(fit_naive)["sodio_obs"]
# Coeficiente atenuado
# Modelo con regression calibration
fit_corrected <- mecor(
 presion ~ MeasError(sodio_obs, reference = sodio_true) + edad,  # reference= acepta NA (gold standard parcial); substitute= requiere datos completos
 data = dat,
 method = "standard" # regression calibration estándar
)
summary(fit_corrected)
# Reporta el coeficiente CORREGIDO con SE ajustado

# ──────────────────────────────────────────────────────────────────────────────
# ■ SIMEX (Simulation Extrapolation)
# ──────────────────────────────────────────────────────────────────────────────

library(simex)
# Para variables continuas con error de medida conocido
# Asumir que la varianza del error es 200^2
# Modelo SIMEX:
# 1) Añadir varianzas crecientes de error al X observado
# 2) Ajustar el modelo en cada caso
# 3) Extrapolar a varianza = 0
fit_simex_base <- lm(presion ~ sodio_obs + edad, data = dat,
 x = TRUE)
fit_simex <- simex(
 fit_simex_base,
 SIMEXvariable = "sodio_obs",
 measurement.error = rep(200, nrow(dat)),
 asymptotic = FALSE,
 fitting.method = "quadratic",
 jackknife.estimation = "quadratic"
)
summary(fit_simex)
plot(fit_simex)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Imputación múltiple para misclassification
# ──────────────────────────────────────────────────────────────────────────────

library(mice)
# Cuando hay validación interna parcial (gold standard medido en subset)
# Datos con NA en las observaciones SIN gold standard
dat_mi <- data.frame(
 Y = c(rep(1, 100), rep(0, 100)),
 A_star = c(rep("expuesto", 60), rep("no", 40),
 rep("expuesto", 40), rep("no", 60)),
 A_true = NA # solo en validación
)
# Asignar A_true en sub-muestra de validación
val_idx <- sample(1:200, 80)
dat_mi$A_true[val_idx] <- ifelse(
 dat_mi$A_star[val_idx] == "expuesto",
 rbinom(80, 1, 0.85), # Se = 0.85
 rbinom(80, 1, 0.10) # 1-Sp = 0.10
)
# MI con mice
imp <- mice(dat_mi, m = 20, method = "logreg",
 predictorMatrix = quickpred(dat_mi))
# Análisis combinado en los m = 20 datasets imputados
fits <- with(imp, glm(Y ~ A_true, family = binomial))
pool(fits) %>% summary(conf.int = TRUE)
# OR ajustado por imputación múltiple


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 24 · Imputación múltiple (MICE)
# ║  Método: mice::mice() · imputación múltiple por ecuaciones encadenadas
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Exploración del missing pattern
# ──────────────────────────────────────────────────────────────────────────────

library(mice)
library(VIM) # visualización adicional
library(dplyr)
# Cargar datos con missing
dat <- nhanes2 # dataset de ejemplo en mice
# 1. EXPLORAR EL PATRÓN DE MISSING
# Resumen numérico
md.pattern(dat)
# Tabla con frecuencias por patrón
# Visualización con VIM
VIM::aggr(dat, col = unname(col_libro[c("azul_med", "rojo")]),
 numbers = TRUE, sortVars = TRUE,
 labels = names(dat))
# Test MCAR de Little
naniar::mcar_test(dat)
# Si p < 0.05: rechazar MCAR (los datos NO son MCAR; pueden ser MAR o MNAR)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Imputación con MICE
# ──────────────────────────────────────────────────────────────────────────────

# 2. IMPUTACIÓN
# Defaults: PMM para continuas, logreg para binarias, polyreg para categóricas
imp <- mice(
 dat,
 m = 30, # número de imputaciones
 maxit = 50, # iteraciones por imputación
 method = NULL, # NULL → detectar automáticamente
 # predictorMatrix omitida: NULL no es válido; mice() construye la matriz automáticamente
 seed = 12345,
 print = FALSE
)
# Especificación manual
methods <- make.method(dat)
# CORRECCIÓN: nhanes2 tiene age/bmi/hyp/chl (no "edad"); age es completa
methods["bmi"] <- "pmm"
methods["hyp"] <- "logreg"
methods["chl"] <- "pmm"
methods["bmi"] <- "pmm"
methods["hyp"] <- "logreg"
methods["chl"] <- "pmm"
pred <- make.predictorMatrix(dat)
pred[, "age"] <- 0  # excluir age como predictor (nhanes2: age es completa)
imp <- mice(dat, m = 30, maxit = 50,
 method = methods, predictorMatrix = pred,
 seed = 12345, print = FALSE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos de convergencia
# ──────────────────────────────────────────────────────────────────────────────

# 3. DIAGNÓSTICOS
# Traceplot de convergencia (CRÍTICO)
plot(imp, main = "Traceplot de convergencia de las cadenas MICE")
# Verificar que las medias y SDs se estabilizan
# Density plot: observado vs imputado
densityplot(imp, ~ bmi |.imp, main = "Densidad: valores observados vs. imputados (BMI)")
# Las curvas deben superponerse aproximadamente
# Stripplot: valores imputados por imputación
stripplot(imp, bmi ~.imp, pch = 20, cex = 1.2, main = "Valores imputados de BMI por imputación")
# Relación bivariada
xyplot(imp, bmi ~ chl |.imp, pch = 20, main = "BMI vs. colesterol — observado e imputado")
# Los puntos imputados deben mostrar patrón similar a los observados

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis y pooling
# ──────────────────────────────────────────────────────────────────────────────

# 4. ANÁLISIS EN CADA IMPUTACIÓN
# Definir el modelo de análisis
fits <- with(imp, glm(hyp ~ age + bmi + chl, family = binomial))
class(fits) # "mira" — lista de M análisis
# 5. COMBINAR CON LAS REGLAS DE RUBIN
pooled <- pool(fits)
summary(pooled, conf.int = TRUE, exponentiate = TRUE)
# Reporta:
# - estimate (β̂_pooled)
# - std.error (SE_pooled = √T)
# - statistic (t-test)
# - df (Barnard-Rubin gl ajustados)
# - p.value
# - 2.5% y 97.5% (IC)
# - exponentiate = TRUE → reporta OR
# Métricas adicionales
pool.fmi <- pool(fits)$pooled$fmi # FMI por coeficiente
pool.lambda <- pool(fits)$pooled$lambda # proporción incremento varianza

# ──────────────────────────────────────────────────────────────────────────────
# ■ Imputación para análisis longitudinal
# ──────────────────────────────────────────────────────────────────────────────

# Datos longitudinales: convertir a formato wide para imputar
library(tidyr)
# Simular dat_long a partir de nhanes2: 3 olas por sujeto
set.seed(2024)
n_suj <- nrow(nhanes2)
dat_long <- do.call(rbind, lapply(1:3, function(w) {
  d <- nhanes2
  d$id   <- seq_len(n_suj)
  d$wave <- w
  # Añadir variación temporal pequeña en continuas
  d$bmi <- d$bmi + rnorm(n_suj, 0, 0.5)
  d$chl <- d$chl + rnorm(n_suj, 0, 5)
  d
}))
dat_wide <- dat_long %>%
 pivot_wider(names_from = wave, values_from = c(bmi, chl, hyp))
# Imputación con auxiliares y efectos longitudinales
imp_long <- mice(
 dat_wide,
 m = 30, maxit = 50,
 method = "pmm",
 seed = 12345,
 print = FALSE
)
# Convertir de vuelta a formato long para análisis
imp_complete <- complete(imp_long, action = "long", include = TRUE)
imp_complete_long <- imp_complete %>%
 pivot_longer(
   cols      = -c(.imp, .id, id),
   names_to  = c(".value", "wave"),
   names_sep = "_"
 )

# ──────────────────────────────────────────────────────────────────────────────
# ■ Imputación con efectos aleatorios (clustered data)
# ──────────────────────────────────────────────────────────────────────────────

library(miceadds)
# Para datos agrupados (multinivel, multicéntricos)
# 2lonly.pmm: imputar variables a nivel agrupado
# 2l.norm: nivel-2 con norm-error
# 2l.continuous (panImpute): imputación multinivel real
# imp_multi <- mice(  # pseudocódigo: requiere datos multinivel reales
#  dat_real,
#  m = 30, maxit = 50,
#  method = method_multi,  # construir con make.method(); binarias: 2l.bin,
#  predictorMatrix = pred_multi,  # col cluster marcada con -2
#  seed = 12345,
#  print = FALSE
#)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis de sensibilidad bajo MNAR (delta-adjustment)
# ──────────────────────────────────────────────────────────────────────────────

# Si se sospecha MNAR, hacer sensibilidad
# Asumir que los valores imputados son sistemáticamente diferentes de los observados:
# X_imp = X_imputed + Δ (donde Δ es un sesgo asumido)
# Implementar con mice via "post" argument
delta <- c(-0.5, 0, 0.5, 1.0) # sesgos asumidos en SD
results_sens <- lapply(delta, function(d) {
 post <- imp$post
 cmd <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", d * sd(dat$bmi, na.rm = TRUE))
 post["bmi"] <- cmd
 imp_d <- mice(dat, m = 30, maxit = 50, post = post,
 seed = 12345, print = FALSE)
 fit_d <- with(imp_d, glm(hyp ~ age + bmi + chl, family = binomial))
 return(pool(fit_d))
})
# Comparar resultados por valor de delta
# CORRECCIÓN: estimate es vector sin nombres; filtrar por term
sapply(results_sens, function(r) { s <- summary(r); s$estimate[s$term == "bmi"] })
