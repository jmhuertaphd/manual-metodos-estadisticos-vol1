# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                   BLOQUE 1 — Modelos de regresión clásicos                   
#                            Código R · Fichas 1–3                             
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
#   Ficha 01 · Regresión logística
#   Ficha 02 · Regresión Poisson y binomial negativa
#   Ficha 03 · Regresión lineal
#
# Contenido: 3 fichas · 22 fragmentos de código · 298 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v90.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 01 · Regresión logística
# ║  Método: glm(family=binomial) · OR, ajuste, diagnósticos
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Firma del método
# ──────────────────────────────────────────────────────────────────────────────

glm(..., family = binomial)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Sintaxis básica con datos reales
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival) # datos clínicos de ejemplo
library(broom) # tidy() para extraer coeficientes
library(performance) # diagnósticos integrales
library(car) # VIF, tests de hipótesis
library(ggplot2)
library(dplyr)
# Datos de ejemplo: PBC (cirrosis biliar primaria)
data(pbc, package = "survival")
dat <- pbc %>%
 filter(!is.na(trt)) %>%
 mutate(
 death = ifelse(status == 2, 1, 0),
 sex = factor(sex, levels = c("m","f")),
 edema_cat = factor(edema, levels = c(0, 0.5, 1),
 labels = c("ninguno","leve","grave")),
 age_decade = age / 10
 )
# Conjunto de covariables (definido UNA vez, reutilizable)
covars <- c("age_decade", "sex", "edema_cat", "bili", "albumin")
# Fórmula explícita (sin construcciones eval/parse)
formula_log <- as.formula(paste("death ~", paste(covars, collapse = " + ")))
# Ajuste GLM
fit <- glm(formula_log, data = dat, family = binomial(link = "logit"))
# Resumen estadístico
summary(fit)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Extracción de OR e IC 95%
# ──────────────────────────────────────────────────────────────────────────────

# Coeficientes con OR e IC perfil de verosimilitud
res <- tidy(fit,
 exponentiate = TRUE,
 conf.int = TRUE,
 conf.method = "profile") %>%
 mutate(across(c(estimate, conf.low, conf.high), ~round(.x, 3)),
 p.value = format.pval(p.value, digits = 3, eps = 0.001))
print(res)
# Test global del modelo (LR test contra modelo nulo)
fit_null <- glm(death ~ 1, data = dat, family = binomial)
anova(fit_null, fit, test = "LRT")
# Test parcial: contribución de edema_cat
drop1(fit, test = "LRT")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Regresión logística condicional (caso-control con matching)
# ──────────────────────────────────────────────────────────────────────────────

library(survival)
# Para datos emparejados se usa clogit() — basado en Cox stratified
# Cada caso emparejado con N controles según variable de matching ('set')
clogit(case ~ exposicion + confusores + strata(set), data = dat_match)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos integrales
# ──────────────────────────────────────────────────────────────────────────────

# Batería completa de diagnósticos
check_model(fit) # gráficos múltiples (binned residuals, influencia, colinealidad)
model_performance(fit) # AIC, BIC, R², RMSE, log_loss, score
check_collinearity(fit) # VIF generalizado (GVIF)
check_outliers(fit) # Cook's D, leverage, residuos estandarizados

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos detallados
# ──────────────────────────────────────────────────────────────────────────────

# --- a) Linealidad del logit (variables continuas) ---
library(splines)
fit_spline <- glm(death ~ ns(age_decade, 4) + sex + edema_cat + bili + albumin,
 data = dat, family = binomial)
anova(fit, fit_spline, test = "LRT") # p-valor pequeño => no linealidad
# Visualización: logit empírico por deciles
dat %>%
 mutate(decil_age = ntile(age_decade, 10)) %>%
 group_by(decil_age) %>%
 summarise(media_age = mean(age_decade),
 logit_emp = log(mean(death) / (1 - mean(death)))) %>%
 ggplot(aes(media_age, logit_emp)) +
 geom_point(size = 3, colour = col_libro["navy"]) +
 geom_smooth(method = "loess", se = FALSE, colour = col_libro["ochre"], linewidth = 1.2) +
 labs(title    = "Linealidad del logit frente a la edad",
 subtitle = "Logit empírico por decil de edad con suavizado LOESS",
 x        = "Edad (décadas)",
 y        = "Logit empírico [log(p/(1-p))]",
 caption  = "Desviaciones de la tendencia lineal sugieren no linealidad") +
 theme_libro()
# --- b) Multicolinealidad ---
car::vif(fit) # GVIF^(1/(2*Df)) > 2.5 indica problema
# --- c) Bondad de ajuste Hosmer-Lemeshow ---
library(ResourceSelection)
hoslem.test(fit$y, fitted(fit), g = 10)
# Alternativas: rms::val.prob() o calibration plot manual con ggplot2.
# --- d) Discriminación: ROC y AUC ---
library(pROC)
roc_obj <- roc(dat$death, fitted(fit), quiet = TRUE)
auc(roc_obj)
ci.auc(roc_obj)
# --- e) Calibración ---
library(rms)
val.prob(fitted(fit), dat$death, m = 30)
# --- f) Observaciones influyentes ---
infl <- influence.measures(fit)
summary(infl)
which(apply(infl$is.inf, 1, any)) # filas marcadas como influyentes
# --- g) Residuos: deviance y Pearson ---
res_dev <- residuals(fit, type = "deviance")
res_pea <- residuals(fit, type = "pearson")
plot(fitted(fit), res_dev,
 xlab = "Predicción (probabilidad)", ylab = "Residuo de devianza",
 main = "Residuos de devianza frente a valores ajustados",
 pch = 19, col = adjustcolor(col_libro["navy"], alpha.f = 0.5),
 col.main = col_libro["navy"], font.main = 2)
abline(h = c(-2, 0, 2), lty = 2, col = col_libro["rojo"])


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 02 · Regresión Poisson y binomial negativa
# ║  Método: glm(family=poisson)/glm.nb() · IRR, sobredispersión
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos de ejemplo y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(MASS) # Insurance, glm.nb()
library(broom) # tidy()
library(performance) # diagnósticos
library(AER) # dispersiontest()
library(DHARMa) # residuos cuantil-aleatorizados
library(sandwich) # varianza robusta
library(lmtest) # coeftest() con SE robustos
library(ggplot2)
library(dplyr)
data(Insurance, package = "MASS")
dat <- Insurance %>%
 mutate(
 Group = factor(Group, ordered = FALSE), # cilindrada
 Age = factor(Age, ordered = FALSE),
 District = factor(District)
 )
# Conjunto de covariables (definido UNA vez)
covars <- c("District", "Group", "Age")
# Fórmula con offset (log de exposición = pólizas-año)
formula_pois <- as.formula(
 paste("Claims ~", paste(covars, collapse = " + "), "+ offset(log(Holders))")
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste del modelo de Poisson
# ──────────────────────────────────────────────────────────────────────────────

# Modelo de Poisson con offset
fit_pois <- glm(formula_pois, data = dat, family = poisson(link = "log"))
summary(fit_pois)
# IRR con IC95% perfil de verosimilitud
res_pois <- tidy(fit_pois,
 exponentiate = TRUE,
 conf.int = TRUE,
 conf.method = "profile") %>%
 filter(term != "(Intercept)") %>%
 mutate(across(c(estimate, conf.low, conf.high), ~round(.x, 3)))
print(res_pois)
# Test global LR
fit_null <- glm(Claims ~ offset(log(Holders)), data = dat, family = poisson)
anova(fit_null, fit_pois, test = "LRT")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comprobación de sobredispersión
# ──────────────────────────────────────────────────────────────────────────────

# 1. Estadístico de dispersión phi
phi_hat <- sum(residuals(fit_pois, type = "pearson")^2) /
 df.residual(fit_pois)
phi_hat # > 1.2 sugiere sobredispersión
# 2. Test formal de Cameron-Trivedi
AER::dispersiontest(fit_pois, trafo = 1) # H1: Var = mu + alpha*mu
AER::dispersiontest(fit_pois, trafo = 2) # H1: Var = mu + alpha*mu^2 (NB2)
# 3. Comparación AIC Poisson vs NB
AIC(fit_pois)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste de la binomial negativa
# ──────────────────────────────────────────────────────────────────────────────

# Modelo NB2 (parametrización por defecto)
fit_nb <- glm.nb(formula_pois, data = dat)
summary(fit_nb)
# Parámetro de dispersión theta (kappa) y su SE
fit_nb$theta
fit_nb$SE.theta
# IRR con IC95%
res_nb <- tidy(fit_nb,
 exponentiate = TRUE,
 conf.int = TRUE,
 conf.method = "profile") %>%
 filter(term != "(Intercept)") %>%
 mutate(across(c(estimate, conf.low, conf.high), ~round(.x, 3)))
print(res_nb)
# Test LR Poisson vs NB (no estándar: alpha = 1/theta = 0 en frontera)
LR <- 2 * (logLik(fit_nb) - logLik(fit_pois))
pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE) / 2
# lmtest::lrtest(fit_pois, fit_nb) da el p-valor SIN el /2.
data.frame(LR = round(as.numeric(LR), 2), p_value = signif(pval, 3))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Quasi-Poisson y varianza robusta
# ──────────────────────────────────────────────────────────────────────────────

# Alternativa 1: Quasi-Poisson (mismos beta, SE inflados por phi)
fit_qp <- glm(formula_pois, data = dat, family = quasipoisson(link = "log"))
summary(fit_qp)
# Alternativa 2: SE robustos sandwich (Poisson modificada de Zou para RR)
coef_robust <- lmtest::coeftest(fit_pois, vcov = sandwich::vcovHC(fit_pois, type = "HC0"))
print(coef_robust)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos integrales
# ──────────────────────────────────────────────────────────────────────────────

# Batería completa con performance
performance::check_model(fit_nb)
performance::check_overdispersion(fit_pois)
performance::check_zeroinflation(fit_pois)
performance::model_performance(fit_nb)
# Comparación entre modelos
performance::compare_performance(fit_pois, fit_qp, fit_nb,
 metrics = c("AIC", "BIC", "RMSE"))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Residuos cuantil-aleatorizados (DHARMa)
# ──────────────────────────────────────────────────────────────────────────────

# Simulación de residuos escalados
sim_res <- DHARMa::simulateResiduals(fit_nb, n = 1000)
# Tests automáticos
DHARMa::testResiduals(sim_res) # KS, dispersión, outliers
DHARMa::testZeroInflation(sim_res) # exceso de ceros
DHARMa::testDispersion(sim_res) # sobredispersión residual
DHARMa::testOutliers(sim_res)
# Gráficos diagnósticos
plot(sim_res) # QQ + residuos vs predicción
DHARMa::plotResiduals(sim_res, dat$Group) # residuos por grupo

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos detallados
# ──────────────────────────────────────────────────────────────────────────────

# --- Multicolinealidad ---
car::vif(fit_nb)
# --- Rootogram (Kleiber & Zeileis, 2016) ---
library(topmodels)
rootogram(fit_pois, main = "Poisson")
rootogram(fit_nb, main = "Binomial Negativa")
# --- Observaciones influyentes ---
infl <- influence.measures(fit_nb)
which(apply(infl$is.inf, 1, any))
# --- Predicciones vs observados ---
dat$pred_nb <- predict(fit_nb, type = "response")
ggplot(dat, aes(pred_nb, Claims)) +
 geom_point(alpha = 0.6, colour = col_libro["navy"]) +
 geom_abline(slope = 1, linetype = "dashed", colour = col_libro["rojo"], linewidth = 1) +
 labs(title    = "Calibración del modelo binomial negativo",
 subtitle = "Recuentos predichos frente a observados",
 x        = "Recuento predicho (binomial negativa)",
 y        = "Recuento observado",
 caption  = "Línea discontinua: calibración perfecta (predicho = observado)") +
 theme_libro()


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 03 · Regresión lineal
# ║  Método: lm() · coeficientes, diagnósticos, multicolinealidad
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Firma del método
# ──────────────────────────────────────────────────────────────────────────────

lm(formula, data)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos de ejemplo y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(MASS) # Boston
library(broom) # tidy() para extraer coeficientes
library(performance) # diagnósticos integrales
library(DHARMa)  # required by performance::check_model()
library(car) # VIF, tests de hipótesis, ncvTest()
library(sandwich) # vcovHC para SE robustos
library(lmtest) # coeftest, bptest
library(ggplot2)
library(dplyr)
data(Boston, package = "MASS")
dat <- Boston %>%
 mutate(
 chas = factor(chas, levels = c(0, 1), labels = c("no", "sí")),
 rad = factor(rad) # accesibilidad a autopistas (categórica)
 )
# Conjunto de covariables (definido UNA vez)
covars <- c("crim", "rm", "age", "dis", "tax", "ptratio", "lstat", "chas")
formula_lin <- as.formula(
 paste("medv ~", paste(covars, collapse = " + "))
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste del modelo lineal
# ──────────────────────────────────────────────────────────────────────────────

# Ajuste OLS
fit <- lm(formula_lin, data = dat)
summary(fit)
# Extracción de coeficientes con IC95%
res <- tidy(fit, conf.int = TRUE) %>%
 mutate(across(c(estimate, std.error, conf.low, conf.high), ~round(.x, 4)),
 p.value = format.pval(p.value, digits = 3, eps = 0.001))
print(res)
# R² y otros estadísticos globales
glance(fit)
# Test F global y test parcial
fit_null <- lm(medv ~ 1, data = dat)
anova(fit_null, fit) # comparación contra modelo nulo
drop1(fit, test = "F") # contribución de cada covariable

# ──────────────────────────────────────────────────────────────────────────────
# ■ Coeficientes estandarizados
# ──────────────────────────────────────────────────────────────────────────────

# Estandarización: coef estandarizados beta*
# (cambio en SD de Y por SD de X) — útil para comparar covariables
dat_std <- dat %>%
 mutate(across(all_of(c("crim", "rm", "age", "dis", "tax", "ptratio", "lstat", "medv")),
 ~as.numeric(scale(.x))))
fit_std <- lm(formula_lin, data = dat_std)
tidy(fit_std, conf.int = TRUE)
# Alternativa con paquete dedicado
library(parameters)
parameters::standardise_parameters(fit, method = "refit")

# ──────────────────────────────────────────────────────────────────────────────
# ■ SE robustos (heterocedasticidad)
# ──────────────────────────────────────────────────────────────────────────────

# HC3: recomendado en muestras pequeñas (Long & Ervin, 2000)
coef_robust <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC3"))
print(coef_robust)
# IC robustos al 95%
lmtest::coefci(fit, vcov = sandwich::vcovHC(fit, type = "HC3"), level = 0.95)
# SE clúster (datos agrupados; ej. agrupado por 'rad')
library(lmtest)
clustered_se <- sandwich::vcovCL(fit, cluster = ~ rad)
coeftest(fit, vcov = clustered_se)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos integrales
# ──────────────────────────────────────────────────────────────────────────────

# Batería completa con performance
performance::check_model(fit)
performance::model_performance(fit)
performance::check_normality(fit)
performance::check_heteroscedasticity(fit)
performance::check_collinearity(fit)
performance::check_outliers(fit)
# plot.lm(): los 4 gráficos diagnósticos clásicos
par(mfrow = c(2, 2))
plot(fit)
par(mfrow = c(1, 1))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos detallados
# ──────────────────────────────────────────────────────────────────────────────

# --- a) Linealidad ---
# Component-residual plots (CR plots): detectan no linealidad por covariable
car::crPlots(fit, smooth = list(span = 0.5))
# --- b) Normalidad de residuos ---
# QQ-plot con bandas de confianza
car::qqPlot(fit, simulate = TRUE, reps = 1000)
# Test formal (sensibles con n grande; uso visual preferido)
shapiro.test(residuals(fit))
# --- c) Homocedasticidad ---
# Test de Breusch-Pagan
lmtest::bptest(fit)
# Test de Cook-Weisberg / Score
car::ncvTest(fit)
# --- d) Independencia ---
# Test de Durbin-Watson (autocorrelación de residuos)
car::durbinWatsonTest(fit)
# --- e) Multicolinealidad ---
car::vif(fit) # VIF; > 5–10 problemático
car::vif(fit, type = "predictor") # GVIF para factores con >1 g.l.
# --- f) Observaciones influyentes ---
infl <- influence.measures(fit)
summary(infl)
which(apply(infl$is.inf, 1, any)) # filas marcadas como influyentes
# Cook's D, DFBETAS, DFFITS, hat values
plot(cooks.distance(fit), type = "h",
 ylab = "Distancia de Cook", xlab = "Índice de observación",
 main = "Observaciones influyentes (distancia de Cook)",
 col = col_libro["navy"], col.main = col_libro["navy"], font.main = 2)
abline(h = 4 / nrow(dat), col = col_libro["rojo"], lty = 2)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Predicción y validación
# ──────────────────────────────────────────────────────────────────────────────

# Predicciones con IC para la media y para nuevas observaciones
new_dat <- data.frame(
 crim = 0.5, rm = 6.5, age = 50, dis = 5,
 tax = 300, ptratio = 18, lstat = 10, chas = "no"
)
predict(fit, newdata = new_dat, interval = "confidence", level = 0.95)
predict(fit, newdata = new_dat, interval = "prediction", level = 0.95)
# Validación cruzada k-fold
library(caret)
ctrl <- trainControl(method = "cv", number = 10)
caret::train(formula_lin, data = dat, method = "lm", trControl = ctrl)
