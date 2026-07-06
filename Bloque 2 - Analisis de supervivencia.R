# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                     BLOQUE 2 — Análisis de supervivencia                     
#                            Código R · Fichas 4–7                             
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
#   Ficha 04 · Kaplan-Meier y log-rank
#   Ficha 05 · Regresión de Cox
#   Ficha 06 · Modelos paramétricos de supervivencia
#   Ficha 07 · Riesgos competitivos
#
# Contenido: 4 fichas · 31 fragmentos de código · 526 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v91.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 04 · Kaplan-Meier y log-rank
# ║  Método: survival::survfit()/survdiff() · curvas, comparación de grupos
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos de ejemplo y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival) # Surv(), survfit(), survdiff()
library(survminer) # ggsurvplot() — gráficos publicables
library(broom) # tidy()
library(dplyr)
library(ggplot2)
data(lung, package = "survival")
dat <- lung %>%
 mutate(
 death = ifelse(status == 2, 1, 0),
 sex = factor(sex, levels = c(1, 2), labels = c("hombre", "mujer")),
 ph_cat = factor(ph.ecog, levels = 0:3,
 labels = c("ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3"))
 ) %>%
 filter(!is.na(ph.ecog))
# Definir el objeto de supervivencia (UNA vez, reutilizable)
surv_obj <- with(dat, Surv(time = time, event = death))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estimador de Kaplan-Meier
# ──────────────────────────────────────────────────────────────────────────────

# KM global (sin estratificación)
fit_km <- survfit(surv_obj ~ 1, data = dat,
 conf.type = "log-log") # IC log-log preferido
print(fit_km)
summary(fit_km, times = c(180, 365, 730)) # Ŝ a 6 m, 1 a, 2 a
# Mediana de supervivencia con IC95%
quantile(fit_km, probs = c(0.25, 0.50, 0.75), conf.int = TRUE)
# Tabla resumen extendida
summary(fit_km)$table # n, eventos, mediana, IC, *rmean*, etc.
# Extracción ordenada con broom
broom::tidy(fit_km) %>% head(10)

# ──────────────────────────────────────────────────────────────────────────────
# ■ KM por grupos y test de log-rank
# ──────────────────────────────────────────────────────────────────────────────

# KM estratificado por sexo
fit_sex <- survfit(surv_obj ~ sex, data = dat, conf.type = "log-log")
# Mediana por grupo
summary(fit_sex)$table
# Test de log-rank (Mantel-Cox)
lr <- survdiff(surv_obj ~ sex, data = dat, rho = 0) # rho=0 → log-rank
print(lr)
# p-valor exacto
1 - pchisq(lr$chisq, df = length(lr$n) - 1)
# Log-rank con peso de Wilcoxon (Peto-Peto): rho = 1
survdiff(surv_obj ~ sex, data = dat, rho = 1)
# Log-rank ESTRATIFICADO (ajustado por una variable categórica)
survdiff(surv_obj ~ sex + strata(ph_cat), data = dat)

# ──────────────────────────────────────────────────────────────────────────────
# ■ RMST (Restricted Mean Survival Time)
# ──────────────────────────────────────────────────────────────────────────────

# RMST con paquete survRM2
library(survRM2)
# tau = horizonte temporal (típicamente el percentil 90 del seguimiento mínimo)
tau_max <- min(max(dat$time[dat$sex == "hombre"]),
 max(dat$time[dat$sex == "mujer"]))
rmst_fit <- rmst2(time = dat$time,
 status = dat$death,
 arm = as.numeric(dat$sex == "mujer"),
 tau = tau_max)
print(rmst_fit)
# Diferencia de RMST entre grupos con IC95%
rmst_fit$unadjusted.result

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación de supuestos
# ──────────────────────────────────────────────────────────────────────────────

# 1. Censura no informativa: comparar características de censurados vs no censurados
dat %>%
 group_by(censurado = ifelse(death == 0, "censurado", "evento")) %>%
 summarise(across(c(age, ph.karno, meal.cal, wt.loss),
 list(media = ~mean(.x, na.rm = TRUE),
 n_na = ~sum(is.na(.x)))))
# 2. Verificación visual de hazards proporcionales (HP) entre grupos
# Si HP se cumple → log(-log Ŝ(t)) vs log(t) deben ser paralelas
plot(fit_sex, fun = "cloglog",
 col = unname(col_libro[c("navy","rojo")]), lwd = 2,
 xlab = "log(tiempo)", ylab = "log(-log Ŝ(t))",
 main = "Diagnóstico visual de hazards proporcionales")
legend("topleft", levels(dat$sex), col = unname(col_libro[c("navy","rojo")]), lwd = 2)
# 3. Curvas KM con tabla "at risk"
survminer::ggsurvplot(
 fit_sex,
 data = dat,
 conf.int = TRUE,
 pval = TRUE, # test de log-rank automático
 pval.method = TRUE,
 risk.table = TRUE,
 risk.table.height = 0.25,
 break.time.by = 180,
 title = "Supervivencia de Kaplan-Meier por sexo",
 xlab = "Tiempo (días)",
 ylab = "Probabilidad de supervivencia",
 legend.title = "Sexo",
 palette = unname(col_libro[c("navy", "rojo")])
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparaciones múltiples y subgrupos
# ──────────────────────────────────────────────────────────────────────────────

# Test de log-rank entre múltiples grupos (>2): omnibus + por pares
survdiff(surv_obj ~ ph_cat, data = dat)
# Comparaciones por pares con corrección de Bonferroni
survminer::pairwise_survdiff(Surv(time, death) ~ ph_cat,
 data = dat,
 p.adjust.method = "bonferroni")
# Forest plot de subgrupos (frecuente en ensayos clínicos)
# (requiere modelo Cox por subgrupo; ver Ficha 5)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 05 · Regresión de Cox
# ║  Método: survival::coxph() · HR, riesgos proporcionales, diagnósticos
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Firma del método
# ──────────────────────────────────────────────────────────────────────────────

survival::coxph(Surv(time, event) ~ X)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos de ejemplo y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival) # coxph(), Surv(), cox.zph()
library(survminer) # ggsurvplot, ggcoxzph, ggcoxdiagnostics
library(broom) # tidy()
library(rms) # cph(), validate(), nomogram
library(timeROC) # AUC dependiente del tiempo
library(dplyr)
library(ggplot2)
data(lung, package = "survival")
dat <- lung %>%
 mutate(
 death = ifelse(status == 2, 1, 0),
 sex = factor(sex, levels = c(1, 2), labels = c("hombre", "mujer")),
 ph_cat = factor(ph.ecog, levels = 0:3,
 labels = c("ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3"))
 ) %>%
 filter(!is.na(ph.ecog), !is.na(wt.loss), !is.na(ph.karno))
# Conjunto de covariables (definido UNA vez)
covars <- c("age", "sex", "ph_cat", "wt.loss", "ph.karno")
formula_cox <- as.formula(
 paste("Surv(time, death) ~", paste(covars, collapse = " + "))
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste del modelo de Cox
# ──────────────────────────────────────────────────────────────────────────────

# Ajuste por verosimilitud parcial (Efron por defecto)
fit <- coxph(formula_cox, data = dat, ties = "efron")
summary(fit)
# Extracción ordenada de HR con IC95%
res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
 mutate(across(c(estimate, conf.low, conf.high), ~round(.x, 3)),
 p.value = format.pval(p.value, digits = 3, eps = 0.001))
print(res)
# Tests globales del modelo (Wald, LR, Score)
glance(fit)
# Test parcial: contribución de cada covariable
drop1(fit, test = "Chisq")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Predicción y supervivencia ajustada
# ──────────────────────────────────────────────────────────────────────────────

# Predicción del riesgo lineal beta·X
dat$lin_pred <- predict(fit, type = "lp")
# Predicción de supervivencia para perfiles concretos
new_profiles <- data.frame(
 age = c(55, 70, 70),
 sex = factor(c("mujer", "hombre", "hombre")),
 ph_cat = factor(c("ECOG 0", "ECOG 1", "ECOG 3")),
 wt.loss = c(0, 5, 15),
 ph.karno = c(90, 80, 60)
)
# Curvas Ŝ(t|X) para cada perfil
sf <- survfit(fit, newdata = new_profiles)
# Visualización con survminer
ggsurvplot(sf, data = new_profiles,
 conf.int = TRUE,
 legend.labs = c("Joven, mujer, buen estado",
 "Mayor, hombre, leve",
 "Mayor, hombre, ECOG 3"),
 legend.title = "Perfil de covariables",
 palette = unname(col_libro[c("navy", "ochre", "rojo")]),
 title = "Supervivencia ajustada por perfil de covariables (Cox)",
 xlab = "Tiempo (días)",
 ylab = "Probabilidad de supervivencia")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Cox estratificado y covariables variables en el tiempo
# ──────────────────────────────────────────────────────────────────────────────

# Cox estratificado por una covariable que viola HP
fit_strat <- coxph(Surv(time, death) ~ age + sex + wt.loss +
 ph.karno + strata(ph_cat),
 data     = dat)
summary(fit_strat)
# Cox con efecto dependiente del tiempo (interacción β·log(t))
# Forma canónica (Therneau): tt(sex) = coeficiente β·log(t); sin reformatear a (start,stop)
fit_tv <- coxph(Surv(time, death) ~ age + sex + ph_cat + wt.loss +
 ph.karno + tt(sex), data = dat,
 tt = function(x, t, ...) as.numeric(x == "mujer") * log(t))
summary(fit_tv)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación del supuesto HP — residuos de Schoenfeld
# ──────────────────────────────────────────────────────────────────────────────

# Test de Schoenfeld por covariable + global
test_zph <- cox.zph(fit, transform = "km")
print(test_zph)
# Salida:
# chisq df p
# age 0.7 1 0.40
# sex 2.1 1 0.15
# ph_cat 5.2 3 0.16
# wt.loss 1.4 1 0.24
# ph.karno 0.9 1 0.34
# GLOBAL 9.8 7 0.20
# Visualización de los residuos por covariable
ggcoxzph(test_zph)
# Si la línea es horizontal → HP se cumple
# Si la línea muestra pendiente → considerar:
# • estratificar por la covariable
# • introducir interacción con función del tiempo
# • cambiar a modelo paramétrico flexible (Ficha 6)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Otros diagnósticos
# ──────────────────────────────────────────────────────────────────────────────

# Linealidad en covariables continuas: residuos de martingala
ggcoxdiagnostics(fit, type = "martingale", linear.predictions = FALSE,
 ox.scale = "linear.predictions")
# Influencia: residuos de devianza
ggcoxdiagnostics(fit, type = "deviance", linear.predictions = FALSE)
# DFBETAS para detectar observaciones influyentes
ggcoxdiagnostics(fit, type = "dfbeta")
# Concordancia (C-index)
fit$concordance["concordance"]
# o equivalentemente:
concordance(fit)
# Validación interna con bootstrap (rms package)
fit_rms <- cph(formula_cox, data = dat, x = TRUE, y = TRUE, surv = TRUE)
validate(fit_rms, B = 200) # Optimismo-corregido C-index, R²D, ...

# ──────────────────────────────────────────────────────────────────────────────
# ■ Discriminación y calibración
# ──────────────────────────────────────────────────────────────────────────────

# AUC dependiente del tiempo (para predicción)
library(timeROC)
roc_t <- timeROC(T = dat$time, delta = dat$death,
 marker = dat$lin_pred,
 cause = 1, weighting = "marginal",
 times = c(180, 365, 730), iid = TRUE)
print(roc_t)
plot(roc_t, time = 365, main = "AUC dependiente del tiempo (IPCW) — t = 365 días")
# Calibración por bootstrap
val <- rms::calibrate(fit_rms, u = 365, B = 200)
plot(val, main = "Calibración por bootstrap — supervivencia a 365 días")


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 06 · Modelos paramétricos de supervivencia
# ║  Método: flexsurv/rstpm2 · Weibull, log-normal, splines
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos de ejemplo
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival) # survreg() — paramétricos clásicos
library(flexsurv) # flexsurvreg(), flexsurvspline() — flexibles
library(rstpm2) # stpm2() — implementación moderna de RP
library(broom) # tidy()
library(ggplot2)
library(dplyr)
data(lung, package = "survival")
dat <- lung %>%
 mutate(
 death = ifelse(status == 2, 1, 0),
 sex = factor(sex, levels = c(1, 2), labels = c("hombre", "mujer")),
 ph_cat = factor(ph.ecog, levels = 0:3,
 labels = c("ECOG 0", "ECOG 1", "ECOG 2", "ECOG 3"))
 ) %>%
 filter(!is.na(ph.ecog), !is.na(wt.loss), !is.na(ph.karno))
# Conjunto de covariables (definido UNA vez)
covars <- c("age", "sex", "ph_cat", "wt.loss", "ph.karno")
formula_surv <- as.formula(
 paste("Surv(time, death) ~", paste(covars, collapse = " + "))
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelos paramétricos clásicos con flexsurv
# ──────────────────────────────────────────────────────────────────────────────

# Ajuste de varias distribuciones para comparar
fit_exp <- flexsurvreg(formula_surv, data = dat, dist = "exp")
fit_wei <- flexsurvreg(formula_surv, data = dat, dist = "weibull")
fit_gomp <- flexsurvreg(formula_surv, data = dat, dist = "gompertz")
fit_lnorm <- flexsurvreg(formula_surv, data = dat, dist = "lnorm")
fit_llog <- flexsurvreg(formula_surv, data = dat, dist = "llogis")
fit_ggam <- flexsurvreg(formula_surv, data = dat, dist = "gengamma")
# Comparación por AIC
aic_tab <- data.frame(
 model = c("exp", "weibull", "gompertz", "lnorm", "llogis", "gengamma"),
 aic = sapply(list(fit_exp, fit_wei, fit_gomp, fit_lnorm, fit_llog, fit_ggam),
 AIC),
 npar = sapply(list(fit_exp, fit_wei, fit_gomp, fit_lnorm, fit_llog, fit_ggam),
 function(f) length(coef(f)))
) %>%
 arrange(aic) %>%
 mutate(delta_aic = aic - min(aic))
print(aic_tab)
# Resumen del mejor modelo
summary(fit_wei, type = "hazard") # h(t) por perfil
summary(fit_wei, type = "survival") # S(t) por perfil

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelos Royston-Parmar con flexsurv
# ──────────────────────────────────────────────────────────────────────────────

# RP con distintos K nudos sobre la escala log-cumulative-hazard (PH)
fit_rp1 <- flexsurvspline(formula_surv, data = dat, k = 1, scale = "hazard")
fit_rp2 <- flexsurvspline(formula_surv, data = dat, k = 2, scale = "hazard")
fit_rp3 <- flexsurvspline(formula_surv, data = dat, k = 3, scale = "hazard")
fit_rp4 <- flexsurvspline(formula_surv, data = dat, k = 4, scale = "hazard")
# Comparación de K nudos por AIC
aic_rp <- data.frame(
 k = 1:4,
 aic = sapply(list(fit_rp1, fit_rp2, fit_rp3, fit_rp4), AIC)
)
print(aic_rp)
# El mejor RP suele estar en K = 2-3
print(fit_rp3)
# HR ajustados (escala "hazard" da PH directo)
res_rp <- as.data.frame(fit_rp3$res.t) # coef en escala log
res_rp$HR <- exp(res_rp$est)
res_rp$HR_lo <- exp(res_rp$L95)
res_rp$HR_hi <- exp(res_rp$U95)
print(res_rp)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Implementación moderna con rstpm2 (equivalente a stpm2 de Stata)
# ──────────────────────────────────────────────────────────────────────────────

library(rstpm2)
# Modelo Royston-Parmar — sintaxis muy similar a stpm2 de Stata
fit_stpm <- stpm2(
 Surv(time, death) ~ age + sex + ph_cat + wt.loss + ph.karno,
 data = dat,
 df = 4 # df = K+1 nudos internos + nudos boundary
)
summary(fit_stpm)
# Predicción de S(t) y h(t) suavizados
new_perfil <- data.frame(
 age = 65, sex = factor("mujer"),
 ph_cat = factor("ECOG 1"), wt.loss = 5, ph.karno = 80
)
# Predicciones con IC
pred_S <- predict(fit_stpm, newdata = new_perfil, type = "surv",
 grid = TRUE, full = TRUE, se.fit = TRUE)
pred_h <- predict(fit_stpm, newdata = new_perfil, type = "hazard",
 grid = TRUE, full = TRUE, se.fit = TRUE)
# HR variable en el tiempo: stpm2 con interacción spline·covariable
dat$sexM <- as.integer(dat$sex == "mujer") # indicador 0/1 para el HR
fit_stpm_tvc <- stpm2(
 Surv(time, death) ~ age + sexM + ph_cat + wt.loss + ph.karno,
 data = dat, df = 4,
 tvc = list(sexM = 2) # spline temporal con 2 df para sexM (mujer vs hombre)
)
summary(fit_stpm_tvc)
# Predicción de HR(t) para sexo
hr_t <- predict(fit_stpm_tvc, newdata = transform(new_perfil, sexM = 0),
 type = "hr", var = "sexM", grid = TRUE, full = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelos AFT clásicos con survreg
# ──────────────────────────────────────────────────────────────────────────────

# survreg es el motor clásico de paramétricos AFT
fit_aft_wei <- survreg(formula_surv, data = dat, dist = "weibull")
fit_aft_lnorm <- survreg(formula_surv, data = dat, dist = "lognormal")
fit_aft_llog <- survreg(formula_surv, data = dat, dist = "loglogistic")
# Time Ratios (TR) con IC
res_aft <- summary(fit_aft_lnorm)$table
TR <- exp(res_aft[, "Value"])
TR_lo <- exp(res_aft[, "Value"] - 1.96 * res_aft[, "Std. Error"])
TR_hi <- exp(res_aft[, "Value"] + 1.96 * res_aft[, "Std. Error"])
data.frame(TR = round(TR, 3),
 IC_inf = round(TR_lo, 3),
 IC_sup = round(TR_hi, 3))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación con Kaplan-Meier
# ──────────────────────────────────────────────────────────────────────────────

# Comparación visual con KM por sexo
km_fit <- survfit(Surv(time, death) ~ sex, data = dat)
# Predicciones del modelo paramétrico
new_dat <- expand.grid(
 sex = factor(c("hombre", "mujer")),
 age = mean(dat$age),
 ph_cat = factor("ECOG 1", levels = levels(dat$ph_cat)),
 wt.loss = mean(dat$wt.loss),
 ph.karno = mean(dat$ph.karno)
)
# Plot
plot(km_fit, col = unname(col_libro[c("navy","rojo")]), lwd = 2,
 xlab = "Tiempo (días)", ylab = "S(t)",
 main = "KM (escalonado) vs Royston-Parmar (suavizado)")
lines(fit_rp3, newdata = new_dat[1,, drop = FALSE], col = col_libro["navy"], lty = 2)
lines(fit_rp3, newdata = new_dat[2,, drop = FALSE], col = col_libro["rojo"], lty = 2)
legend("topright", c("KM hombre", "KM mujer", "RP hombre", "RP mujer"),
 col = unname(col_libro[c("navy","rojo","navy","rojo")]),
 lwd = 2, lty = c(1, 1, 2, 2))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Residuos de Cox-Snell
# ──────────────────────────────────────────────────────────────────────────────

# Si el modelo es correcto, los residuos de Cox-Snell siguen exp(1)
# y -log(Ŝ_KM(rᵢ)) frente a rᵢ debe ser una recta de pendiente 1
cox_snell <- function(fit, dat) {
 # Residuo de Cox-Snell del sujeto i: H(tᵢ|xᵢ) = -log S(tᵢ|xᵢ),
 # riesgo acumulado evaluado en el tiempo observado de cada sujeto.
 vapply(seq_len(nrow(dat)), function(i) {
  s <- summary(fit, newdata = dat[i, , drop = FALSE],
              type = "cumhaz", t = dat$time[i], ci = FALSE)
  s[[1]]$est[1]
 }, numeric(1))
}
cs_res <- cox_snell(fit_rp3, dat)
km_cs <- survfit(Surv(cs_res, dat$death) ~ 1)
plot(km_cs$time, -log(km_cs$surv), pch = 16,
 col = adjustcolor(col_libro["navy"], alpha.f = 0.6),
 xlab = "Residuo de Cox-Snell", ylab = "-log Ŝ(rᵢ)",
 main = "QQ-plot Cox-Snell — si el modelo es correcto, sigue la diagonal")
abline(0, 1, col = col_libro["rojo"], lwd = 2)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Selección de número de nudos en RP
# ──────────────────────────────────────────────────────────────────────────────

# Búsqueda sistemática de K óptimo
k_grid <- 1:5
aic_grid <- numeric(length(k_grid))
for (k in k_grid) {
 fit_k <- flexsurvspline(formula_surv, data = dat, k = k, scale = "hazard")
 aic_grid[k] <- AIC(fit_k)
}
data.frame(k = k_grid, AIC = aic_grid,
 delta_AIC = aic_grid - min(aic_grid))
# Test LR entre modelos anidados (df adicional por nudo)
# Manualmente:
LR <- 2 * (logLik(fit_rp3) - logLik(fit_rp2))
pchisq(as.numeric(LR), df = 1, lower.tail = FALSE) # k=3 vs k=2

# ──────────────────────────────────────────────────────────────────────────────
# ■ Validación con bootstrap (rstpm2 / pec)
# ──────────────────────────────────────────────────────────────────────────────

# C-index dependiente del tiempo
library(pec)
library(prodlim)
cox_ref <- coxph(formula_surv, data = dat, x = TRUE, y = TRUE)
# pec/riskRegression no traen método para flexsurv: lo definimos
predictSurvProb.flexsurvreg <- function(object, newdata, times, ...) {
 s <- summary(object, newdata = newdata, type = "survival", t = times, ci = FALSE)
 t(vapply(s, function(d) d$est, numeric(length(times))))
}
cindex_t <- pec::cindex(
 list("Cox" = cox_ref, "RP-3" = fit_rp3),
 formula = Surv(time, death) ~ 1,
 data = dat,
 eval.times = c(180, 365, 730)
)
print(cindex_t)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 07 · Riesgos competitivos
# ║  Método: cmprsk/tidycmprsk · subdistribución de Fine-Gray
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos de ejemplo y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival)
library(cmprsk) # cuminc(), crr() — gold standard clásico
library(tidycmprsk) # tidycuminc(), tidycrr() — interfaz moderna
library(survminer) # ggcompetingrisks
library(mstate) # marco multistate más general
library(rstpm2) # paramétrico flexible para CR
library(broom)
library(dplyr)
library(ggplot2)
# Datos de ejemplo
data(bmtcrr, package = "casebase")
# bmtcrr (casebase): Status 0=censura, 1=recaída, 2=evento competidor (muerte sin recaída)
# ftime (meses), Sex (F/M), D (enfermedad ALL/AML), Phase, Source (injerto), Age
bmt <- transform(bmtcrr,
                 cause = Status, # numérico 0/1/2 (cmprsk::cuminc()/crr())
                 group = D,      # enfermedad ALL/AML (alternativa: Source = injerto)
                 age   = Age,
                 sex   = Sex)
dat <- bmt %>%
 mutate(
 # 0 = censura, 1 = recaída, 2 = muerte sin recaída
 cause = factor(Status, levels = 0:2,
 labels = c("censura", "recaida", "muerte_sr"))
 )
table(dat$cause)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estimador no paramétrico (Aalen-Johansen)
# ──────────────────────────────────────────────────────────────────────────────

# Con cmprsk::cuminc — sintaxis clásica
ci_global <- cmprsk::cuminc(ftime = bmt$ftime, fstatus = bmt$cause)
print(ci_global)
# Por grupo (variable group: enfermedad ALL/AML)
ci_grupo <- cmprsk::cuminc(ftime = bmt$ftime, fstatus = bmt$cause,
 group = bmt$group)
print(ci_grupo)
# Test de Gray automático en la salida
# Visualización con survminer
ggcompetingrisks(ci_grupo,
 multiple_panels = FALSE,
 xlab = "Tiempo (meses)",
 ylab = "Incidencia acumulada",
 title = "CIF por enfermedad (ALL/AML)") +
 scale_color_manual(values = unname(col_libro[c("navy","ochre","azul_med","rojo")])) +
 theme_libro()
# Equivalente moderno con tidycmprsk
tcuminc <- tidycmprsk::cuminc(Surv(ftime, cause) ~ group, data = dat)
tidycmprsk::tbl_cuminc(tcuminc, times = c(12, 24, 60),
 label_header = "**{time} meses**")
library(ggsurvfit)
ggcuminc(tcuminc, outcome = "recaida") +
 add_confidence_interval() +
 add_risktable() +
 ggplot2::labs(x = "Tiempo (meses)", y = "Incidencia acumulada de recaída")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelo cause-specific Cox
# ──────────────────────────────────────────────────────────────────────────────

# Recodificar para cause-specific: cada evento se trata como evento;
# los demás como censura
covars <- c("age", "sex", "group") # ajustar a las variables del dataset
formula_cs <- as.formula(
 paste("Surv(ftime, cause == 'recaida') ~", paste(covars, collapse = " + "))
)
# Cox cause-specific para recaída
fit_cs_relapse <- coxph(formula_cs, data = dat)
summary(fit_cs_relapse)
# Y para el evento competidor (muerte sin recaída)
formula_cs2 <- as.formula(
 paste("Surv(ftime, cause == 'muerte_sr') ~", paste(covars, collapse = " + "))
)
fit_cs_death <- coxph(formula_cs2, data = dat)
summary(fit_cs_death)
# Extracción ordenada
res_cs <- broom::tidy(fit_cs_relapse, exponentiate = TRUE, conf.int = TRUE)
print(res_cs)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelo de Fine-Gray (sub-distribution)
# ──────────────────────────────────────────────────────────────────────────────

# crr() exige diseño matricial — no acepta fórmula directa con factores
# Construir la matriz de diseño manualmente
covars_num <- model.matrix(~ age + sex + group, data = dat)[, -1]
fit_fg <- cmprsk::crr(ftime = bmt$ftime,
 fstatus = bmt$cause, # numérico 0/1/2
 cov1 = covars_num,
 failcode = 1, # recaída (evento de interés)
 cencode = 0)  # censura
summary(fit_fg)
# Más cómodo con tidycmprsk (interfaz moderna sobre crr)
fit_fg2 <- tidycmprsk::crr(
 Surv(ftime, cause) ~ age + sex + group,
 data = dat,
 failcode = "recaida"
)
print(fit_fg2)
broom::tidy(fit_fg2, exponentiate = TRUE, conf.int = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelo paramétrico flexible para riesgos competitivos
# ──────────────────────────────────────────────────────────────────────────────

# rstpm2::stpm2 con cause-specific approach (más natural)
# Ajustar un modelo flexible para cada hazard cause-specific
fit_stpm_relapse <- stpm2(
 Surv(ftime, cause == "recaida") ~ age + sex + group,
 data = dat, df = 4
)
fit_stpm_death <- stpm2(
 Surv(ftime, cause == "muerte_sr") ~ age + sex + group,
 data = dat, df = 4
)
# La CIF se obtiene por integración numérica de los hazards
# (rstpm2 facilita esto con predict)
# Para CIF directa con splines flexibles:
# library(flexsurv)
# fit_fg_spline <- flexsurvspline(
# Surv(ftime, cause == "recaida") ~ age + sex + group,
# data = dat, k = 3,
# anc = list(gamma1 = ~ group) # spline-coef varía por grupo
# )

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación visual: CIF observada vs predicha
# ──────────────────────────────────────────────────────────────────────────────

# Comparación gráfica de la CIF estimada por el modelo
# frente a la no paramétrica de Aalen-Johansen
new_dat <- expand.grid(
 age = mean(dat$age),
 sex = unique(dat$sex)[1],
 group = unique(dat$group)
)
# Predicción de CIF según Fine-Gray (cmprsk::predict)
pred_fg <- predict(fit_fg, cov1 = model.matrix(~ age + sex + group, new_dat)[, -1])
plot(pred_fg, lty = 1:nrow(new_dat), col = unname(col_libro[c("navy","rojo")]),
 xlab = "Tiempo (meses)", ylab = "F̂(t) — Recaída",
 main = "CIF predichas — Fine-Gray")
# Superponer estimaciones AJ por grupo (no paramétricas)
ci_aj <- cmprsk::cuminc(ftime = bmt$ftime, fstatus = bmt$cause, group = bmt$group)
plot(ci_aj, col = unname(col_libro[c("navy","rojo")]), lty = 2,
 curvlab = c(""), xlab = "", ylab = "")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación de proporcionalidad en Fine-Gray
# ──────────────────────────────────────────────────────────────────────────────

# El test de Schoenfeld estándar no se aplica directamente a Fine-Gray
# Alternativa: introducir interacción tiempo·covariable
# (cmprsk no la soporta nativamente; usar rstpm2 o splitting)
# Vía 1: comparación con modelo cause-specific
# Si HR_cs ≈ sHR para covariables que no afectan al competidor → coherencia interna
# Si difieren mucho → reportar ambos análisis
# Vía 2: validación cruzada temporal — partir el seguimiento en periodos
dat_split <- survSplit(
 Surv(ftime, cause == "recaida") ~.,
 data = dat,
 cut = quantile(dat$ftime[dat$cause == "recaida"], c(0.33, 0.66))
)
# Y ajustar el modelo permitiendo coeficientes distintos por intervalo

# ──────────────────────────────────────────────────────────────────────────────
# ■ C-index para riesgos competitivos
# ──────────────────────────────────────────────────────────────────────────────

# C-index de Wolbers (2009) — adaptado para competing risks
library(pec)
# Necesita un modelo que prediga F_k(t)
# FG con riskRegression::FGR (tiene predictEventProb; pec::cindex lo requiere)
library(riskRegression)
fgr <- FGR(Hist(ftime, cause) ~ age + sex + group, data = bmt, cause = 1)
# C-index dependiente del tiempo sobre F_k(t)
cindex_cr <- pec::cindex(
 list("Fine-Gray" = fgr),
 formula = Hist(ftime, cause) ~ 1,
 cause = 1,
 data = bmt,
 eval.times = c(12, 24, 60)
)
print(cindex_cr)
# --- Horizonte fijo: AUC(t) y Brier (recomendado) ---
# El C-index no es "propio" para riesgos a t fijo (Blanche et al., 2019):
# un modelo mal especificado puede superar al verdadero. La AUC dependiente
# del tiempo sí lo es y el Brier/IPA añade calibración.
score_fg <- Score(
 list("Fine-Gray" = fgr),
 formula = Hist(ftime, cause) ~ 1,
 data = bmt, cause = 1,
 times = c(12, 24, 60),
 metrics = c("AUC", "Brier"),
 summary = "ipa", null.model = TRUE, conf.int = TRUE
)
score_fg
