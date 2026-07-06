# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                           BLOQUE 6 — Meta-análisis                           
#                           Código R · Fichas 19–21                            
# ==============================================================================
#
# Fichas incluidas en este bloque:
#   Ficha 19 · Meta-análisis de efectos fijos y aleatorios
#   Ficha 20 · Meta-análisis dosis-respuesta
#   Ficha 21 · Meta-regresión y análisis de subgrupos
#
# Contenido: 3 fichas · 24 fragmentos de código · 453 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v90.docx.
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



# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 19 · Meta-análisis de efectos fijos y aleatorios
# ║  Método: meta::metabin()/metafor::rma() · síntesis cuantitativa; heterogeneidad y sesgo
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Paquete meta de Schwarzer
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(meta) # Schwarzer — sintaxis intuitiva, output completo
library(metafor) # Viechtbauer — más flexible, gold standard
library(dmetar) # diagnósticos avanzados
library(dplyr)
library(ggplot2)
# Datos: ejemplo de 10 estudios con eventos en grupo tratamiento y control
dat <- data.frame(
 estudio = c("Smith 2015", "Johnson 2016", "García 2017", "Müller 2018",
 "Tanaka 2019", "Pereira 2019", "Anderson 2020", "Rossi 2021",
 "Lindqvist 2022", "Kim 2023"),
 eventos_t = c(25, 18, 32, 42, 15, 22, 38, 12, 20, 16),
 total_t = c(125, 90, 155, 210, 75, 112, 190, 60, 100, 82),
 eventos_c = c(35, 22, 47, 56, 18, 35, 49, 11, 24, 22),
 total_c = c(125, 90, 155, 210, 75, 113, 190, 60, 100, 83),
 year = c(2015, 2016, 2017, 2018, 2019, 2019, 2020, 2021, 2022, 2023)
)
# Meta-análisis binario con OR (efectos aleatorios por defecto)
ma <- metabin(
 event.e = eventos_t, # eventos en experimental
 n.e = total_t,
 event.c = eventos_c, # eventos en control
 n.c = total_c,
 studlab = estudio,
 data = dat,
 sm = "OR", # "OR", "RR", "RD"
 method = "MH", # Mantel-Haenszel para FE
 method.tau = "REML", # REML para τ² (recomendado)
 method.random.ci = "HK", # Hartung-Knapp para IC (recomendado)
 common = TRUE, # mostrar también FE
 random = TRUE # mostrar RE
)
print(ma)
# Salida incluye:
# - Efectos individuales con IC
# - Síntesis FE y RE
# - τ², I², Q
# - IC predicción

# ──────────────────────────────────────────────────────────────────────────────
# ■ Forest plot canónico
# ──────────────────────────────────────────────────────────────────────────────

# Forest plot básico
forest(ma)
# Forest plot personalizado
forest(ma,
 sortvar = TE, # ordenar por efecto
 leftcols = c("studlab"),
 rightcols = c("effect", "ci", "w.random"),
 label.left = "Favorece tratamiento",
 label.right = "Favorece control",
 col.diamond = col_libro["verde"],
 col.square = col_libro["navy"],
 common = TRUE,
 random = TRUE,
 print.tau2 = TRUE,
 print.I2 = TRUE,
 text.predict = "IC predicción 95%",
 prediction = TRUE)
# Guardar como PDF
pdf("forest_plot.pdf", width = 11, height = 7)
forest(ma, prediction = TRUE)
dev.off()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Funnel plot y tests de sesgo
# ──────────────────────────────────────────────────────────────────────────────

# Funnel plot
# Funnel plot con contornos — ver bloque metafor más abajo (funnel(fit_re))
# Egger's test
metabias(ma, method.bias = "linreg")
# Begg's test (correlación de rango)
metabias(ma, method.bias = "rank")
# Trim-and-fill (Duval-Tweedie)
tf <- trimfill(ma)
summary(tf)
funnel(tf)
title("Funnel plot — corrección trim-and-fill (Duval-Tweedie)", col.main = col_libro["navy"], font.main = 2)
# Compara la estimación original con la corregida por estudios «faltantes»

# ──────────────────────────────────────────────────────────────────────────────
# ■ Meta-análisis para efectos continuos
# ──────────────────────────────────────────────────────────────────────────────

# Datos: media y SD por grupo
dat_cont <- data.frame(
 estudio = c("A", "B", "C", "D", "E"),
 mean_t = c(85, 82, 90, 78, 86),
 sd_t = c(12, 14, 11, 13, 10),
 n_t = c(50, 65, 80, 55, 70),
 mean_c = c(90, 88, 95, 84, 92),
 sd_c = c(13, 15, 12, 14, 11),
 n_c = c(50, 65, 80, 55, 70)
)
# Diferencia de medias estandarizada (Cohen's d / Hedges' g)
ma_cont <- metacont(
 n.e = n_t, mean.e = mean_t, sd.e = sd_t,
 n.c = n_c, mean.c = mean_c, sd.c = sd_c,
 studlab = estudio,
 data = dat_cont,
 sm = "SMD", # "MD" para diferencia bruta
 method.smd = "Hedges", # corrección Hedges para muestras pequeñas
 method.tau = "REML",
 method.random.ci = "HK"
)
print(ma_cont)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Meta-análisis de HR (tiempo a evento)
# ──────────────────────────────────────────────────────────────────────────────

# Para HR: se trabaja con log(HR) y SE(log(HR))
dat_hr <- data.frame(
 estudio = c("A", "B", "C", "D"),
 log_hr = log(c(0.75, 0.82, 0.68, 0.79)),
 se_log_hr = c(0.12, 0.10, 0.15, 0.08)
)
ma_hr <- metagen(
 TE = log_hr,
 seTE = se_log_hr,
 studlab = estudio,
 data = dat_hr,
 sm = "HR",
 method.tau = "REML",
 method.random.ci = "HK",
 exponentiate = TRUE # devolver HR (no log-HR)
)
print(ma_hr)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Paquete metafor (más flexible)
# ──────────────────────────────────────────────────────────────────────────────

library(metafor)
# Cálculo de tamaño del efecto y SE
dat_es <- escalc(measure = "OR",
 ai = eventos_t, n1i = total_t,
 ci = eventos_c, n2i = total_c,
 data = dat)
# Modelo RE con REML
fit_re <- rma(yi, vi, data = dat_es,
 method = "REML",
 test = "knha") # Hartung-Knapp
print(fit_re)
# Reporta τ², I², H², R², QE, p-value
# Forest plot (col/border del diamante resumen con la paleta del libro)
forest(fit_re,
 slab = dat$estudio,
 atransf = exp, # exponentiar para mostrar OR
 at = log(c(0.5, 1, 2)),
 col    = col_libro["navy"],   # diamante resumen
 border = col_libro["navy"],
 xlab = "Odds Ratio (escala logarítmica)",
 header = c("Estudio", "OR [IC 95%]"))
title("Forest plot — modelo de efectos aleatorios", col.main = col_libro["navy"], font.main = 2)
# Funnel y test
funnel(fit_re,
  level   = c(90, 95, 99),
  shade   = c("#CBD5E0", "#A0AEC0", "#718096"),
  refline = 0, atransf = exp,
  at      = log(c(0.25, 0.5, 1, 2, 4)),
  xlab    = "Odds Ratio", legend = TRUE)
title("Funnel plot — contornos de significación (90/95/99%)", col.main = col_libro["navy"], font.main = 2)
regtest(fit_re) # Egger
ranktest(fit_re) # Begg
# Intervalo de predicción
predict(fit_re, transf = exp)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis de influencia
# ──────────────────────────────────────────────────────────────────────────────

# Influencia de cada estudio (leave-one-out)
inf <- metainf(ma, pooled = "random")
forest(inf)
title("Análisis de influencia — leave-one-out", col.main = col_libro["navy"], font.main = 2)
# Estudios potencialmente influyentes
inf2 <- influence(fit_re)
plot(inf2)
# Reporta: residual estandarizado, hat, Cook's distance, dfbetas, weight

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis acumulativo
# ──────────────────────────────────────────────────────────────────────────────

# Meta-análisis acumulativo por año
cum <- metacum(ma, pooled = "random", sortvar = dat$year)
forest(cum)
title("Meta-análisis acumulativo por año de publicación", col.main = col_libro["navy"], font.main = 2)
# Muestra cómo la estimación cambia conforme se añaden estudios

# ──────────────────────────────────────────────────────────────────────────────
# ■ Sensibilidad a métodos
# ──────────────────────────────────────────────────────────────────────────────

# Comparar distintos estimadores de τ²
methods <- c("DL", "REML", "ML", "EB", "HE", "SJ")
results <- lapply(methods, function(m) {
 rma(yi, vi, data = dat_es, method = m)
})
# Tabla comparativa
sapply(results, function(x) c(beta = x$b, tau2 = x$tau2, ci.lb = x$ci.lb, ci.ub = x$ci.ub))


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 20 · Meta-análisis dosis-respuesta
# ║  Método: dosresmeta::dosresmeta() · RR continuo por dosis; splines restringidos
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y formato
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(dosresmeta) # Crippa & Orsini — gold standard
library(metafor) # Viechtbauer — utilidades de meta-análisis
library(rms) # Harrell — splines restringidos
library(mvmeta) # Gasparrini — meta-análisis multivariante
library(ggplot2)
library(dplyr)
# Datos típicos: una fila por categoría-estudio
# Columnas: id (estudio), dose (nivel), cases (eventos), n (total), logrr, se
# Ejemplo: alcohol y cáncer de mama
# Datos sintéticos representativos (sustituye al CSV externo):
# Dataset dosis-respuesta simulado (sustituye al CSV externo)
set.seed(2024)
n_studies <- 8; doses <- c(0, 5, 15, 30, 45)
dat_dr <- do.call(rbind, lapply(1:n_studies, function(i) {
  n_per     <- round(runif(length(doses), 500, 2000))
  base_rate <- runif(1, 0.05, 0.15)
  true_rr   <- exp(0.015 * doses)
  cases     <- rbinom(length(doses), n_per,
                      base_rate * true_rr / mean(true_rr))
  se_i      <- runif(length(doses), 0.06, 0.14)
  logrr_i   <- log(true_rr) + rnorm(length(doses), 0, se_i)
  logrr_i[1] <- 0; se_i[1] <- 0   # categoría referencia
  data.frame(id = i, dose = doses, cases = cases,
             n = n_per, logrr = logrr_i, se = se_i)
}))
head(dat_dr)
# Estructura esperada:
# id dose cases n logrr se
# 1 0 120 1500 0.000 0.000 <- categoría referencia
# 1 5 145 1450 0.106 0.080
# 1 15 185 1300 0.328 0.090
# 1 30 240 1100 0.531 0.110
# 2 0 55 700 0.000 0.000
# 2 8 72 680 0.182 0.075
# ...
# Importante: para cada estudio, el RR de referencia es 1.00 (log = 0)
# La dose se asigna al punto medio de la categoría

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelo lineal con dosresmeta
# ──────────────────────────────────────────────────────────────────────────────

# Modelo lineal (RR proporcional a la dosis)
fit_lin <- dosresmeta(
 formula = logrr ~ dose, # modelo lineal en log
 type = "ir", # incidence rate (puede ser "cc" o "ci")
 id = id,
 se = se,
 cases = cases,
 n = n,
 data = dat_dr,
 method = "reml" # efectos aleatorios REML
)
summary(fit_lin)
# Reporta:
# - β̂ con SE robusto e IC 95%
# - τ² (heterogeneidad entre-estudios)
# - I² estandarizado
# - p-valor de heterogeneidad
# Efecto por incremento (ej: por 10 g/día)
exp(coef(fit_lin) * 10)
exp(confint(fit_lin) * 10)
# Interpretación: RR aumenta exp(β·10) veces por cada 10 g/día

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelo con splines restringidos
# ──────────────────────────────────────────────────────────────────────────────

# Splines restringidos con 4 nudos (recomendado)
# Nudos en percentiles 5, 35, 65, 95 de la distribución conjunta de dosis
# 3 nudos sobre dosis positivas (4 valores únicos: 5,15,30,45)
knots <- quantile(dat_dr$dose[dat_dr$dose > 0], c(0.10, 0.50, 0.90))
fit_spl <- dosresmeta(
 formula = logrr ~ rcs(dose, knots),
 type = "ir",
 id = id,
 se = se,
 cases = cases,
 n = n,
 data = dat_dr,
 method = "reml"
)
summary(fit_spl)
# Test de no-linealidad: H₀ = todos los términos no lineales = 0
# Equivalente a comparar fit_spl con fit_lin
waldtest(Sigma = vcov(fit_spl), b = coef(fit_spl), Terms = 2:length(coef(fit_spl)))
# Predicciones a lo largo de la dosis
new_dose <- seq(0, 50, by = 1)
pred <- predict(fit_spl, newdata = data.frame(dose = new_dose),
 xref = 0, # referencia: 0 g/día
 exp = TRUE) # devolver RR (no log)
# Tabla de RR predicho con IC
head(pred)
# Columnas: dose, pred, ci.lb, ci.ub

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización con ggplot2
# ──────────────────────────────────────────────────────────────────────────────

# Curva dosis-respuesta canónica
# Construir data frame de predicción con columna dose explícita
new_doses <- seq(0, max(dat_dr$dose), length.out = 100)
pred_mat  <- predict(fit_spl, newdata = data.frame(dose = new_doses),
                     xref = 0, expo = TRUE)
pred <- data.frame(dose  = new_doses,
                   pred  = pred_mat$pred,
                   ci.lb = pred_mat$ci.lb,
                   ci.ub = pred_mat$ci.ub)
ggplot(pred, aes(x = dose, y = pred)) +
 geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub),
 fill = col_libro["azul_med"], alpha = 0.20) +
 geom_line(color = col_libro["navy"], linewidth = 1.3) +
 geom_hline(yintercept = 1, linetype = "dashed", color = col_libro["gris"]) +
 # Puntos de los estudios individuales
 geom_point(data = dat_dr %>% filter(dose > 0),
 aes(y = exp(logrr)),
 color = col_libro["ochre"], size = 2.5, alpha = 0.7) +
 scale_y_log10(breaks = c(0.7, 1, 1.5, 2, 3),
 limits = c(0.7, 3)) +
 labs(title    = "Curva dosis-respuesta del consumo de alcohol",
 subtitle = "Spline cúbico restringido con IC 95%; puntos = estudios individuales",
 x        = "Consumo de alcohol (g/día)",
 y        = "RR (escala logarítmica)",
 caption  = "Línea discontinua: RR = 1 (ausencia de efecto)") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Selección del mejor modelo
# ──────────────────────────────────────────────────────────────────────────────

# Comparar modelos (esquema; "..." son los argumentos comunes ya usados
# en fit_lin/fit_spl: type="ir", id=id, se=se, cases=cases, n=n, data=dat_dr,
# method="reml"; NO ejecutar estas líneas literalmente, "..." no es código R)
fit_lin <- dosresmeta(logrr ~ dose, ...)
fit_quad <- dosresmeta(logrr ~ dose + I(dose^2), ...)
fit_spl3 <- dosresmeta(logrr ~ rcs(dose, 3), ...)
fit_spl4 <- dosresmeta(logrr ~ rcs(dose, 4), ...)
fit_spl5 <- dosresmeta(logrr ~ rcs(dose, 5), ...)
# Comparación por AIC
AIC(fit_lin, fit_quad, fit_spl3, fit_spl4, fit_spl5)
# Test no-linealidad lineal vs spline
waldtest(Sigma = vcov(fit_spl4), b = coef(fit_spl4), Terms = 2:length(coef(fit_spl4)))
# H₀: el modelo lineal es suficiente

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelo one-stage (Crippa-Orsini 2019)
# ──────────────────────────────────────────────────────────────────────────────

# Para datos IPD parciales o muchos estudios
fit_1s <- dosresmeta(
 formula = logrr ~ rcs(dose, c(5, 15, 25, 40)),
 type = "ir",
 id = id,
 se = se,
 cases = cases,
 n = n,
 data = dat_dr,
 method = "ml",
 proc = "1stage" # una etapa
)
summary(fit_1s)
# La principal ventaja: permite efectos aleatorios sobre TODA la curva
# (cada estudio tiene su propia curva, no solo su propio intercepto/intensidad)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Predicción para dosis específicas
# ──────────────────────────────────────────────────────────────────────────────

# Tabla de RR para dosis clave
dose_key <- c(5, 10, 20, 30, 50)
dose_key <- c(0, dose_key)  # incluir referencia (dose=0) requerida por predict.dosresmeta
pred_key <- predict(fit_spl, newdata = data.frame(dose = dose_key),
 xref = 0, exp = TRUE)
pred_key <- as.data.frame(pred_key); pred_key$dose <- dose_key
pred_key <- pred_key[pred_key$dose > 0, ]; print(pred_key)
# Reporta RR a cada dosis con IC 95%
# Encontrar la dosis a la que RR alcanza un umbral (ej: 1.5)
threshold_dose <- approx(pred$pred, pred$dose, xout = 1.5)$y
cat("RR = 1.5 a:", round(threshold_dose, 1), "g/día\n")
# Encontrar nadir (mínimo) de la curva (para J-shape)
nadir_idx <- which.min(pred$pred)
cat("Nadir en:", pred$dose[nadir_idx],
 "con RR =", round(pred$pred[nadir_idx], 2), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 21 · Meta-regresión y análisis de subgrupos
# ║  Método: metafor::rma(mods=) · moderadores continuos y categóricos; bubble plot
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y meta-análisis base
# ──────────────────────────────────────────────────────────────────────────────

library(meta)
library(metafor)
library(dplyr)
library(ggplot2)
# Datos: 18 estudios sobre tratamiento antihipertensivo y eventos CV
dat <- data.frame(
 estudio = paste0("Estudio_", 1:18),
 yi = c(-0.18, -0.25, -0.32, -0.40, -0.15, -0.22, -0.30, -0.45,
 -0.10, -0.20, -0.50, -0.38, -0.28, -0.42, -0.05, -0.35,
 -0.48, -0.55),
 vi = c(0.012, 0.008, 0.015, 0.009, 0.020, 0.011, 0.014, 0.007,
 0.022, 0.013, 0.006, 0.010, 0.016, 0.008, 0.025, 0.012,
 0.005, 0.009),
 edad_media = c(58, 62, 70, 75, 45, 55, 68, 72,
 42, 55, 78, 70, 60, 73, 38, 65, 76, 80),
 pct_mujeres = c(45, 52, 48, 58, 55, 48, 50, 60,
 42, 48, 55, 52, 50, 56, 60, 47, 58, 65),
 duracion = c(2.0, 2.5, 3.0, 4.0, 1.5, 2.0, 2.8, 3.5,
 1.0, 2.0, 4.5, 3.5, 2.5, 4.0, 1.5, 3.0, 4.0, 5.0),
 tipo_control = factor(c("placebo", "placebo", "activo", "activo", "placebo",
 "placebo", "activo", "activo", "placebo", "placebo",
 "activo", "activo", "placebo", "activo", "placebo",
 "activo", "activo", "activo"))
)
# Meta-análisis base sin moderadores
fit_base <- rma(yi, vi, data = dat, method = "REML")
summary(fit_base)
# τ², I², Q, etc.

# ──────────────────────────────────────────────────────────────────────────────
# ■ Meta-regresión con un moderador continuo
# ──────────────────────────────────────────────────────────────────────────────

# Meta-regresión: ¿el efecto depende de la edad media?
fit_age <- rma(yi, vi, mods = ~ edad_media, data = dat,
 method = "REML", test = "knha")
summary(fit_age)
# Reporta:
# - β para edad_media con SE, t y p-valor
# - τ²_residual (heterogeneidad que QUEDA)
# - R² estandarizado: % de τ² explicada
# Calcular R² explícitamente
tau2_base <- fit_base$tau2
tau2_age <- fit_age$tau2
R2_age <- (tau2_base - tau2_age) / tau2_base
cat("R² =", round(R2_age * 100, 1), "%\n")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Bubble plot canónico
# ──────────────────────────────────────────────────────────────────────────────

# Bubble plot con metafor
regplot(fit_age,
 xlab = "Edad media del estudio (años)",
 ylab = "log(OR)",
 col = col_libro["azul_med"],
 bg = col_libro["navy"],
 main = "Meta-regresión: log(OR) frente a edad media (bubble plot)",
 psize = 1.5, # tamaño base de las burbujas
 las = 1)
# Versión con ggplot2 (más personalizable)
dat$weight <- 1 / dat$vi
ggplot(dat, aes(x = edad_media, y = yi, size = weight)) +
 geom_point(alpha = 0.55, color = col_libro["navy"]) +
 geom_abline(intercept = coef(fit_age)["intrcpt"],
 slope = coef(fit_age)["edad_media"],
 color = col_libro["ochre"], linewidth = 1.2) +
 geom_hline(yintercept = 0, linetype = "dashed", color = col_libro["gris"]) +
 scale_size(range = c(2, 12), name = "Peso (1/varianza)") +
 labs(title    = "Meta-regresión: efecto según la edad media del estudio",
 subtitle = "Tamaño del punto proporcional al peso; recta = ajuste meta-regresivo",
 x        = "Edad media del estudio (años)",
 y        = "log(OR)",
 caption  = "Línea discontinua: log(OR) = 0 (ausencia de efecto)") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Meta-regresión con múltiples moderadores
# ──────────────────────────────────────────────────────────────────────────────

# Modelo completo con varios moderadores
fit_full <- rma(yi, vi,
 mods = ~ edad_media + pct_mujeres + duracion + tipo_control,
 data = dat,
 method = "REML",
 test = "knha")
summary(fit_full)
# Reporta β para cada moderador, R² conjunto, test de moderación global
# Comparación con modelo base
anova(fit_base, fit_full)
# Test de la inclusión conjunta de moderadores
# Importante: evitar sobreajuste con pocos estudios
# Cochrane recomienda ≥ 10 estudios por moderador

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis por subgrupos (moderador categórico)
# ──────────────────────────────────────────────────────────────────────────────

# Análisis estratificado por tipo de control
fit_subgrp <- update(fit_base,
 mods = ~ tipo_control - 1) # sin intercepto: una estimación por grupo
summary(fit_subgrp)
# Alternativa con meta::metabin
ma <- metagen(TE = yi, seTE = sqrt(vi), studlab = estudio, data = dat,
 method.tau = "REML", method.random.ci = "HK")
update(ma, subgroup = tipo_control)
# El test Q_between determina si los subgrupos difieren significativamente

# ──────────────────────────────────────────────────────────────────────────────
# ■ Forest plot estratificado por subgrupos
# ──────────────────────────────────────────────────────────────────────────────

# Categorizar la edad media en tres grupos
dat$grupo_edad <- cut(dat$edad_media,
 breaks = c(0, 50, 65, 100),
 labels = c("Jóvenes (<50)", "Adultos (50-65)",
 "Mayores (>65)"))
ma_subgrp <- metagen(TE = yi, seTE = sqrt(vi),
 studlab = estudio,
 data = dat,
 method.tau = "REML", method.random.ci = "HK")
# Forest plot estratificado
forest(update(ma_subgrp, subgroup = grupo_edad),
 sortvar = TE,
 leftcols = c("studlab"),
 rightcols = c("effect", "ci"),
 label.left = "Favorece tratamiento",
 label.right = "Favorece control",
 common = FALSE, random = TRUE,
 print.subgroup.labels = TRUE,
 test.subgroup = TRUE)
title("Forest plot estratificado por grupo de edad", col.main = col_libro["navy"], font.main = 2)
# Reporta Q_between en el panel inferior

# ──────────────────────────────────────────────────────────────────────────────
# ■ Detección de la falacia ecológica
# ──────────────────────────────────────────────────────────────────────────────

# Si se dispone de IPD, se puede comparar β_between vs β_within
# Simulación de IPD: el efecto modificación A:X varía SISTEMÁTICAMENTE
# con una covariable de nivel-estudio (edad media), para ilustrar la falacia ecológica
set.seed(2024)
n_est <- 6; n_i <- 80
edad_estudio <- c(40, 50, 55, 60, 65, 70) # covariable de nivel-estudio
beta_AX_estudio <- 0.02 * (edad_estudio - 55) # interacción A:X varía con la edad
ipd_data <- do.call(rbind, lapply(1:n_est, function(s) {
  Xi <- rnorm(n_i); Ai <- rbinom(n_i, 1, 0.4)
  data.frame(estudio = s, edad_estudio_s = edad_estudio[s], A = Ai, X = Xi,
    Y = rbinom(n_i, 1, plogis(-1 + 0.5*Ai + 0.3*Xi + beta_AX_estudio[s]*Ai*Xi)))
}))
ipd_data$estudio <- factor(ipd_data$estudio)
# Modelo IPD con efecto modificación
library(lme4)
fit_ipd <- glmer(
 Y ~ A + X + A:X + (1 + A | estudio),
 family = binomial,
 data = ipd_data       # ipd_data ahora definida arriba
)
summary(fit_ipd)
# El coeficiente de A:X es el efecto modificación INDIVIDUAL (β_within)
# Comparar con β de la meta-regresión (β_between)
# Modelo two-stage para descomponer:
# 1) Dentro de cada estudio: efecto modificación específico
# 2) Meta-analizar los efectos modificación individuales
fit_AX_por_estudio <- lapply(1:n_est, function(s) {
  d_s <- subset(ipd_data, estudio == s)
  fit_s <- glm(Y ~ A * X, data = d_s, family = binomial)
  c(beta_AX = unname(coef(fit_s)["A:X"]),
    se_AX = unname(summary(fit_s)$coefficients["A:X", "Std. Error"]))
})
dat_AX <- as.data.frame(do.call(rbind, fit_AX_por_estudio))
dat_AX$edad_media <- edad_estudio
# Meta-regresión sobre la covariable de nivel-estudio (edad media):
# el coeficiente de edad_media es β_between (efecto APARENTE a nivel-estudio)
fit_metareg_AX <- rma(yi = beta_AX, sei = se_AX, mods = ~ edad_media,
 data = dat_AX, method = "REML", test = "knha")
summary(fit_metareg_AX)
# Comparación explícita β_within (glmer, individuo) vs β_between (meta-regresión):
round(fixef(fit_ipd)["A:X"], 4)         # β_within
round(coef(fit_metareg_AX)["edad_media"], 4)  # β_between
# Si difieren en magnitud o dirección, se evidencia la falacia ecológica:
# inferir el efecto individual a partir del efecto entre estudios (o viceversa) puede ser engañoso

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis acumulativo por moderador
# ──────────────────────────────────────────────────────────────────────────────

# Meta-análisis acumulativo por año de publicación
dat$year <- 1995 + cumsum(rep(1, 18)) # ejemplo
fit_cum <- rma.uni(yi, vi, data = dat, method = "REML")
cum_results <- cumul(fit_cum, order = dat$year)
forest(cum_results, transf = exp)
title("Meta-análisis acumulativo por año de publicación", col.main = col_libro["navy"], font.main = 2)
# Muestra cómo la estimación cambia conforme se añaden estudios cronológicamente
