# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                         BLOQUE 5 — Datos especiales                          
#                           Código R · Fichas 15–18                            
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
#   Ficha 15 · Modelos mixtos (efectos aleatorios)
#   Ficha 16 · Ecuaciones de estimación generalizadas (GEE)
#   Ficha 17 · Regresión ZIP y ZINB
#   Ficha 18 · Análisis de mediación
#
# Contenido: 4 fichas · 25 fragmentos de código · 531 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v90.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 15 · Modelos mixtos (efectos aleatorios)
# ║  Método: lme4::lmer()/glmer() · efectos fijos y aleatorios; BLUPs, ICC, selección
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y exploración inicial
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(lme4) # lmer() y glmer() — gold standard
library(lmerTest) # tests con df ajustados (Satterthwaite)
library(nlme) # lme() — alternativa con estructuras de correlación
library(performance) # check_model(), icc(), r2()
library(broom.mixed) # tidy() para modelos mixtos
library(ggplot2)
library(dplyr)
# Datos longitudinales simulados: 300 sujetos x 6 visitas
set.seed(2024)
n_subj <- 300
n_visits <- 6
dat <- expand.grid(id = 1:n_subj, visita = 0:(n_visits - 1)) %>%
 arrange(id, visita)
# Efectos basales (entre sujetos)
dat <- dat %>%
 group_by(id) %>%
 mutate(
 edad_basal = rnorm(1, 55, 10),
 sexo = sample(c("hombre", "mujer"), 1),
 tratado = sample(c(0, 1), 1, prob = c(0.5, 0.5)),
 u0 = rnorm(1, 0, 6), # variabilidad intercepto
 u1 = rnorm(1, 0, 1.2) # variabilidad pendiente
 ) %>%
 ungroup() %>%
 mutate(
 # Generar Y con efecto fijo y aleatorio
 Y = (52 + u0) + (-2.5 + u1) * visita +
 (-1.8) * tratado + (-1.0) * tratado * visita +
 0.6 * (edad_basal - 55)/10 +
 rnorm(n(), 0, 1.5)
 )
# Exploración: spaghetti plot
ggplot(dat, aes(x = visita, y = Y, group = id)) +
 geom_line(alpha = 0.2, colour = col_libro["gris"]) +
 geom_smooth(aes(group = NULL), method = "lm", se = FALSE,
 color = col_libro["rojo"], linewidth = 1.5) +
 labs(title    = "Trayectorias individuales a lo largo del seguimiento",
 subtitle = "Cada línea gris es un sujeto; la línea roja es la tendencia media (MCO)",
 x        = "Visita",
 y        = "Respuesta (Y)",
 caption  = "La variabilidad entre sujetos motiva el uso de efectos aleatorios") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste con lme4::lmer()
# ──────────────────────────────────────────────────────────────────────────────

# Conjunto único de covariables fijas
covars_fijas <- c("visita", "tratado", "visita:tratado",
 "edad_basal", "sexo")
# Modelo con random intercept + random slope para 'visita'
fit_lmm <- lmerTest::lmer(
 Y ~ visita * tratado + edad_basal + sexo + (1 + visita | id),
 data = dat,
 REML = TRUE
)
summary(fit_lmm)
# Salida incluye:
# - Fixed effects: β̂ con SE, t y p-valor (Satterthwaite por lmerTest)
# - Random effects: σ²_intercepto, σ²_pendiente, correlación
# - Residual: σ²_ε
# Extracción ordenada de efectos fijos
broom.mixed::tidy(fit_lmm, effects = "fixed", conf.int = TRUE)
# Componentes de varianza
broom.mixed::tidy(fit_lmm, effects = "ran_pars")
# ICC (proporción de varianza explicada por el cluster)
performance::icc(fit_lmm)
# Reporta ICC ajustado (recomendado) y no ajustado
# R² (Nakagawa-Schielzeth: marginal y condicional)
performance::r2(fit_lmm)
# R²_marginal: solo efectos fijos
# R²_condicional: fijos + aleatorios

# ──────────────────────────────────────────────────────────────────────────────
# ■ BLUPs y predicciones individualizadas
# ──────────────────────────────────────────────────────────────────────────────

# Extracción de BLUPs (efectos aleatorios por sujeto)
blups <- ranef(fit_lmm)$id
head(blups)
# Columnas: (Intercept) y visita (los random effects)
# Predicciones individuales
dat$Y_pred_pop <- predict(fit_lmm, re.form = NA) # solo fijos
dat$Y_pred_ind <- predict(fit_lmm) # fijos + BLUPs
# Visualización
ggplot(dat, aes(x = visita, group = id)) +
 geom_line(aes(y = Y, color = "Observado"), alpha = 0.2) +
 geom_line(aes(y = Y_pred_ind, color = "BLUP individual"), alpha = 0.5) +
 geom_line(aes(y = Y_pred_pop, color = "Población"), linewidth = 1.2) +
 scale_color_manual(values = c("Observado" = unname(col_libro["gris"]),
 "BLUP individual" = unname(col_libro["navy"]),
 "Población" = unname(col_libro["rojo"]))) +
 labs(title    = "Predicciones del modelo mixto: BLUP frente a media poblacional",
 subtitle = "Los BLUP encogen las trayectorias individuales hacia la media",
 x        = "Visita",
 y        = "Respuesta (Y)",
 color    = "Trayectoria",
 caption  = "BLUP: mejores predictores lineales insesgados (empirical Bayes)") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación de modelos y selección
# ──────────────────────────────────────────────────────────────────────────────

# Comparar estructuras de efectos aleatorios (usar ML, no REML)
# Re-ajustar con lmerTest::lmer() → clase lmerModLmerTest, compatible con PBmodcomp
fit_ri    <- lmerTest::lmer(Y ~ visita + tratado + (1 | id), data = dat, REML = FALSE)
fit_ri_rs <- lmerTest::lmer(Y ~ visita + tratado + (1 + visita | id), data = dat, REML = FALSE)
anova(fit_ri, fit_ri_rs)
# LR test: ¿mejora añadir random slope?
# CUIDADO: el test es 0.5·χ²₁ + 0.5·χ²₂ por la frontera
# Bootstrap paramétrico (más fiable)
library(pbkrtest)
PBmodcomp(fit_ri_rs, fit_ri, nsim = 500)
# Comparación por AIC/BIC
AIC(fit_ri, fit_ri_rs)
BIC(fit_ri, fit_ri_rs)

# ──────────────────────────────────────────────────────────────────────────────
# ■ GLMM: extensión a desenlaces no normales
# ──────────────────────────────────────────────────────────────────────────────

# Para desenlaces binarios (logística mixta)
library(glmmTMB)
# Ejemplo: binario por visita
dat$Y_bin <- as.numeric(dat$Y > 45)
fit_glmm_logit <- glmer(
 Y_bin ~ visita * tratado + edad_basal + (1 + visita | id),
 data = dat,
 family = binomial(link = "logit"),
 control = glmerControl(optimizer = "bobyqa")
)
# Extracción y OR
broom.mixed::tidy(fit_glmm_logit, effects = "fixed",
 exponentiate = TRUE, conf.int = TRUE)
# Para recuentos (Poisson o NB mixta)
fit_glmm_pois <- glmmTMB(
 Y_bin ~ visita * tratado + edad_basal + (1 + visita | id),
 data = dat,
 family = poisson
)
# Para sobredispersión, usar negative binomial
fit_glmm_nb <- glmmTMB(
 Y_bin ~ visita * tratado + edad_basal + (1 + visita | id),
 data = dat,
 family = nbinom2
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos esenciales
# ──────────────────────────────────────────────────────────────────────────────

# Diagnóstico completo con performance
library(performance)
check_model(fit_lmm)
# Genera 6 paneles: residuos, normalidad, homocedasticidad,
# multicolinealidad, influencia, normalidad de efectos aleatorios
# Verificaciones manuales clave:
# NOTA sobre el dispositivo gráfico: plot(fit_lmm) usa lattice/grid, mientras que
# qqnorm() usa el sistema base. Ambos NO se limpian mutuamente en el mismo device,
# por lo que deben ejecutarse por separado (cada uno en su propia página/ventana)
# para evitar que las figuras se solapen. Los dos Q-Q base se agrupan aparte.

# 1) Residuos vs valores ajustados (homocedasticidad) — lattice, página propia
print(plot(fit_lmm, main = "Residuos vs. valores ajustados (modelo mixto)"))

# 2-3) Q-Q de residuos y de BLUPs — base graphics, en una figura de 2 paneles
op_qq <- par(mfrow = c(1, 2))
qqnorm(resid(fit_lmm),
 main = "Normalidad de los residuos (Q-Q normal)",
 xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales de los residuos")
qqline(resid(fit_lmm), col = col_libro["rojo"], lwd = 2)
qqnorm(ranef(fit_lmm)$id[, "(Intercept)"],
 main = "Normalidad de los efectos aleatorios (Q-Q normal)",
 xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales de los BLUPs")
qqline(ranef(fit_lmm)$id[, "(Intercept)"], col = col_libro["rojo"], lwd = 2)
par(op_qq)
# 4) Caterpillar plot de los BLUPs
library(lattice)
# dotplot.ranef.mer usa 'main' como lógico (nombre del factor), no como título:
# main="..." da error. Se extrae el objeto trellis del factor 'id' y se titula
# con update(), conservando relation="free" (escala propia por panel) y fijando
# la escala del eje abajo (alternating = FALSE).
bp_ranef <- dotplot(ranef(fit_lmm, condVar = TRUE))$id
update(bp_ranef,
 main = "Caterpillar plot de los efectos aleatorios (BLUPs ± IC 95%)",
 xlab = "Desviación respecto a la media poblacional",
 scales = list(alternating = FALSE, x = list(relation = "free")))
# 5) Convergencia (importante en GLMM)
# Si fit_lmm tiene warning de convergencia:
# start solo acepta theta (parámetros de varianza); fixef no es válido
ss <- getME(fit_lmm, "theta")
fit_lmm_2 <- update(fit_lmm, start = ss,
 control = lmerControl(optCtrl = list(maxfun = 2e5)))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Estructuras de correlación residual
# ──────────────────────────────────────────────────────────────────────────────

# Ejemplo con AR(1) usando nlme
library(nlme)
fit_ar1 <- lme(
 Y ~ visita * tratado + edad_basal + sexo,
 random = ~ 1 + visita | id,
 correlation = corAR1(form = ~ visita | id),
 data = dat
)
summary(fit_ar1)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 16 · Ecuaciones de estimación generalizadas (GEE)
# ║  Método: geepack::geeglm() · GEE para datos correlacionados; estructuras de correlación
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y ajuste con geepack::geeglm
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(geepack) # geeglm() — gold standard
library(gee) # gee() — más antiguo, sintaxis distinta
library(dplyr)
library(broom) # tidy()
library(ggplot2)
# Datos longitudinales: 300 sujetos × 6 visitas (mismo ejemplo Ficha 15)
set.seed(2024)
n_subj <- 300
n_visits <- 6
dat <- expand.grid(id = 1:n_subj, visita = 0:(n_visits - 1)) %>%
 arrange(id, visita)
dat <- dat %>%
 group_by(id) %>%
 mutate(
 edad_basal = rnorm(1, 55, 10),
 sexo = sample(c("hombre", "mujer"), 1),
 tratado = sample(c(0, 1), 1, prob = c(0.5, 0.5)),
 u0 = rnorm(1, 0, 6),
 u1 = rnorm(1, 0, 1.2)
 ) %>%
 ungroup() %>%
 mutate(
 Y = (52 + u0) + (-2.5 + u1) * visita +
 (-1.8) * tratado + (-1.0) * tratado * visita +
 0.6 * (edad_basal - 55) / 10 +
 rnorm(n(), 0, 1.5),
 Y_bin = as.numeric(Y > 45)
 )
# REQUISITO: los datos deben estar ordenados por id (clave para identificar clusters)
dat <- dat %>% arrange(id, visita)

# ──────────────────────────────────────────────────────────────────────────────
# ■ GEE con desenlace continuo
# ──────────────────────────────────────────────────────────────────────────────

# Conjunto único de covariables
covars <- c("visita", "tratado", "visita:tratado",
 "edad_basal", "sexo")
# Modelo GEE con correlación exchangeable
fit_gee_exch <- geeglm(
 Y ~ visita * tratado + edad_basal + sexo,
 id = id, # variable que identifica el cluster
 family = gaussian,
 corstr = "exchangeable",
 data = dat
)
summary(fit_gee_exch)
# Salida incluye:
# - Estimaciones de β con SE robusto y p-valor
# - Estimación del parámetro de correlación (α)
# - Phi (escala)
# Tabla limpia
broom::tidy(fit_gee_exch, conf.int = TRUE)
# Distintas estructuras para comparar (ANÁLISIS DE SENSIBILIDAD)
fit_gee_ind <- update(fit_gee_exch, corstr = "independence")
fit_gee_ar1 <- update(fit_gee_exch, corstr = "ar1")
fit_gee_unstr <- update(fit_gee_exch, corstr = "unstructured")
# Comparación: los β̂ deben ser CASI IDÉNTICOS (robustez del sandwich)
comparison <- rbind(
 cbind(modelo = "Independence", tidy(fit_gee_ind)),
 cbind(modelo = "Exchangeable", tidy(fit_gee_exch)),
 cbind(modelo = "AR(1)", tidy(fit_gee_ar1)),
 cbind(modelo = "Unstructured", tidy(fit_gee_unstr))
) %>% filter(term == "tratado")
print(comparison, digits = 3)

# ──────────────────────────────────────────────────────────────────────────────
# ■ GEE con desenlace binario
# ──────────────────────────────────────────────────────────────────────────────

# Modelo logístico marginal con GEE
fit_gee_logit <- geeglm(
 Y_bin ~ visita * tratado + edad_basal + sexo,
 id = id,
 family = binomial(link = "logit"),
 corstr = "exchangeable",
 data = dat
)
summary(fit_gee_logit)
# OR con IC robustos
broom::tidy(fit_gee_logit, exponentiate = TRUE, conf.int = TRUE)
# La interpretación es MARGINAL: «el OR poblacional pasa de OR_0 a OR_0·exp(β)»

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación GEE vs GLMM
# ──────────────────────────────────────────────────────────────────────────────

# Modelo equivalente con GLMM (Ficha 15)
library(lme4)
fit_glmm_logit <- glmer(
 Y_bin ~ visita * tratado + edad_basal + sexo + (1 | id),
 family = binomial,
 data = dat
)
# Comparación: GEE (marginal) vs GLMM (condicional)
res_gee <- tidy(fit_gee_logit, exponentiate = TRUE, conf.int = TRUE)
res_glmm <- tidy(fit_glmm_logit, effects = "fixed",
 exponentiate = TRUE, conf.int = TRUE)
# Los β̂ marginales (GEE) suelen ser MENORES en magnitud que los condicionales (GLMM)
# Esto NO es un error: refleja la distinta pregunta de investigación
print(dplyr::select(res_gee,  term, estimate, conf.low, conf.high))
print(dplyr::select(res_glmm, term, estimate, conf.low, conf.high))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Selección de estructura con QIC
# ──────────────────────────────────────────────────────────────────────────────

# QIC y CIC para comparar estructuras
QIC(fit_gee_ind)
QIC(fit_gee_exch)
QIC(fit_gee_ar1)
QIC(fit_gee_unstr)
# El QIC más bajo sugiere mejor ajuste (similar a AIC)
# CICu (penaliza solo efectos fijos) para seleccionar la media
# Conviene reportar:
# - β̂ con varias estructuras (mostrar robustez)
# - QIC para justificar la estructura principal
# - SE naive vs sandwich (debe haber poco cambio si Rᵢ es correcta)

# ──────────────────────────────────────────────────────────────────────────────
# ■ GEE con offset (tasas de Poisson)
# ──────────────────────────────────────────────────────────────────────────────

# Para análisis de tasas con persona-tiempo (Ficha 2 y 10)
# Simulamos persona_tiempo (años de seguimiento por observación)
set.seed(42)
dat$persona_tiempo <- runif(nrow(dat), 0.5, 2.0)
fit_gee_poi <- geeglm(
 Y ~ tratado + edad_basal + sexo + offset(log(persona_tiempo)),
 id = id,
 family = poisson(link = "log"),
 corstr = "exchangeable",
 data = dat
)
broom::tidy(fit_gee_poi, exponentiate = TRUE, conf.int = TRUE)
# exp(β) = IRR marginal


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 17 · Regresión ZIP y ZINB
# ║  Método: pscl::zeroinfl()/hurdle() · conteos con exceso de ceros; rootogram
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y exploración inicial
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(pscl) # zeroinfl(), hurdle()
library(MASS) # glm.nb()
library(glmmTMB) # zero-inflated con efectos aleatorios
library(topmodels) # rootogram(), pit
library(performance) # check_overdispersion, check_zeroinflation
library(dplyr)
library(ggplot2)
library(broom)
# Datos simulados: hospitalizaciones anuales en una cohorte
set.seed(2024)
N <- 2000
dat <- data.frame(
 id = 1:N,
 edad = rnorm(N, 55, 12),
 sexo = factor(sample(c("hombre", "mujer"), N, replace = TRUE)),
 tratamiento = rbinom(N, 1, 0.4),
 comorbilidad = rbinom(N, 1, 0.30)
)
# Generar el desenlace con mezcla:
# - 40% son "no susceptibles" (siempre 0)
# - 60% siguen Poisson con λ que depende de covariables
no_susceptible <- rbinom(N, 1, 0.40)
lambda <- exp(0.5 + 0.03 * (dat$edad - 55) +
 0.40 * dat$comorbilidad - 0.30 * dat$tratamiento)
dat$Y <- ifelse(no_susceptible == 1, 0,
 rpois(N, lambda))
# Exploración del exceso de ceros
table(dat$Y == 0)
mean(dat$Y == 0) # proporción de ceros observada
# Esperada bajo Poisson(λ̂)
lambda_hat <- mean(dat$Y)
exp(-lambda_hat) # esperada de Poisson clásica
# Histograma observado
ggplot(dat, aes(x = Y)) +
 geom_histogram(binwidth = 1, fill = col_libro["navy"],
 color = "white", alpha = 0.85) +
 labs(title    = "Distribución del número de hospitalizaciones",
 subtitle = "Exceso de ceros y sobredispersión frente a la Poisson",
 x        = "Hospitalizaciones anuales",
 y        = "Frecuencia",
 caption  = "El pico en cero motiva los modelos inflados de ceros (ZIP/ZINB)") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste de los modelos competidores
# ──────────────────────────────────────────────────────────────────────────────

# Conjunto único de covariables
covars <- c("edad", "sexo", "comorbilidad", "tratamiento")
form_count <- as.formula(paste("Y ~", paste(covars, collapse = " + ")))
# 1) Poisson clásico (baseline)
fit_pois <- glm(form_count, data = dat, family = poisson)
summary(fit_pois)
# 2) Binomial negativa
fit_nb <- glm.nb(form_count, data = dat)
summary(fit_nb)
# 3) ZIP (mismas covariables en ambos componentes)
# Sintaxis: y ~ count_covars | zero_covars
fit_zip <- zeroinfl(Y ~ edad + sexo + comorbilidad + tratamiento |
 edad + sexo + comorbilidad + tratamiento,
 data = dat, dist = "poisson")
summary(fit_zip)
# 4) ZINB
fit_zinb <- zeroinfl(Y ~ edad + sexo + comorbilidad + tratamiento |
 edad + sexo + comorbilidad + tratamiento,
 data = dat, dist = "negbin")
summary(fit_zinb)
# 5) Hurdle Poisson
fit_hurdle_p <- hurdle(Y ~ edad + sexo + comorbilidad + tratamiento |
 edad + sexo + comorbilidad + tratamiento,
 data = dat, dist = "poisson", zero.dist = "binomial")
# 6) Hurdle NB
fit_hurdle_nb <- hurdle(Y ~ edad + sexo + comorbilidad + tratamiento |
 edad + sexo + comorbilidad + tratamiento,
 data = dat, dist = "negbin", zero.dist = "binomial")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación de modelos
# ──────────────────────────────────────────────────────────────────────────────

# Tabla comparativa por AIC y BIC
AIC(fit_pois, fit_nb, fit_zip, fit_zinb, fit_hurdle_p, fit_hurdle_nb)
BIC(fit_pois, fit_nb, fit_zip, fit_zinb, fit_hurdle_p, fit_hurdle_nb)
# Test de Vuong (con cautela)
pscl::vuong(fit_pois, fit_zip)
pscl::vuong(fit_zip, fit_zinb) # NO anidados
# LR test para modelos ANIDADOS
lmtest::lrtest(fit_pois, fit_nb) # Poisson vs NB
lmtest::lrtest(fit_zip, fit_zinb) # ZIP vs ZINB (NB vs Poisson en el recuento)
# Detección de exceso de ceros con performance
performance::check_zeroinflation(fit_nb)
# Reporta:
# - Ceros observados vs predichos
# - Si los predichos < observados → considerar ZI
# Sobredispersión
performance::check_overdispersion(fit_pois)
# Si dispersion ratio >> 1 → usar NB en lugar de Poisson

# ──────────────────────────────────────────────────────────────────────────────
# ■ Interpretación de los coeficientes del modelo ZIP/ZINB
# ──────────────────────────────────────────────────────────────────────────────

# Los modelos ZI tienen DOS sets de coeficientes
summary(fit_zinb)
# Salida:
# Count model coefficients (negbin with log link): 
# β para el componente de recuento (exp(β) = IRR entre susceptibles)
# Zero-inflation model coefficients (binomial with logit link):
# γ para el componente de inflación (exp(γ) = OR de SER no susceptible)
# Extraer coeficientes con IC
# broom::tidy() no tiene método para zeroinfl;
# usar parameters::model_parameters() (dep. de performance, ya cargado)
library(parameters)
model_parameters(fit_zinb, component = "conditional",
                 exponentiate = TRUE, ci = 0.95)  # componente recuento
model_parameters(fit_zinb, component = "zero_inflated",
                 exponentiate = TRUE, ci = 0.95)  # componente inflación
# Predicciones
# E[Y | X] = (1 - π) * λ
pred_mean <- predict(fit_zinb, type = "response")
# Probabilidad de cero
pred_zero <- predict(fit_zinb, type = "zero") # solo la componente π
# Probabilidades de cada recuento
pred_prob <- predict(fit_zinb, type = "prob")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnóstico con rootogram
# ──────────────────────────────────────────────────────────────────────────────

# Rootogram: la herramienta diagnóstica canónica
library(topmodels)
par(mfrow = c(2, 2))
rootogram(fit_pois, main = "Poisson clásico")
rootogram(fit_nb, main = "Binomial negativa")
rootogram(fit_zip, main = "ZIP")
rootogram(fit_zinb, main = "ZINB")
# Buen ajuste: barras terminan cerca del eje 0
# Mal ajuste en ceros: barras NEGATIVAS en y=0 (Poisson clásico)
# Buen ajuste en ceros pero mal en la cola: NB sin ZI
# Buen ajuste global: ZIP o ZINB

# ──────────────────────────────────────────────────────────────────────────────
# ■ ZINB con efectos aleatorios (datos correlacionados)
# ──────────────────────────────────────────────────────────────────────────────

# Para datos correlacionados (medidas repetidas, multicéntricos)
library(glmmTMB)
# ZIP con random intercept por sujeto
fit_zinb_re <- glmmTMB(
 Y ~ edad + sexo + comorbilidad + tratamiento + (1 | id),
 zi = ~ edad + sexo + comorbilidad + tratamiento,
 family = nbinom2,
 data = dat
)
summary(fit_zinb_re)
# Ventajas: permite combinar zero-inflated con jerarquía (Ficha 15)
# La componente zi acepta cualquier modelo del usuario


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 18 · Análisis de mediación
# ║  Método: mediation/CMAverse/regmedint · NDE, NIE, CDE; E-value para confusión M-Y
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y exploración inicial
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(mediation) # Imai, Keele, Tingley — paquete clásico
library(CMAverse) # VanderWeele — moderno y completo
library(regmedint) # Yoshida & Mathur — análisis paramétrico
library(dplyr)
library(broom)
# Datos de ejemplo: estrés (A) → insomnio (M) → depresión (Y)
set.seed(2024)
N <- 2000
dat <- data.frame(
 id = 1:N,
 edad = rnorm(N, 50, 12),
 sexo = factor(sample(c("hombre", "mujer"), N, replace = TRUE)),
 educ = factor(sample(c("baja", "media", "alta"), N, replace = TRUE))
)
# A: estrés crónico (binaria)
dat$A <- rbinom(N, 1, plogis(-1 + 0.02 * (dat$edad - 50) +
 ifelse(dat$sexo == "mujer", 0.3, 0)))
# M: insomnio (binaria, afectado por A)
dat$M <- rbinom(N, 1, plogis(-1.5 + 1.2 * dat$A + 0.02 * (dat$edad - 50)))
# Y: depresión (binaria)
# Efecto directo de A + efecto de M + efectos basales
dat$Y <- rbinom(N, 1, plogis(-2.5 + 0.5 * dat$A + 0.8 * dat$M +
 0.02 * (dat$edad - 50)))
# Conjunto de covariables
covars <- c("edad", "sexo", "educ")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Método clásico: paquete mediation (Imai-Keele-Tingley)
# ──────────────────────────────────────────────────────────────────────────────

# Paso 1: modelo del mediador
fit_M <- glm(M ~ A + edad + sexo + educ, data = dat, family = binomial)
# Paso 2: modelo del desenlace
fit_Y <- glm(Y ~ A + M + edad + sexo + educ, data = dat, family = binomial)
# Análisis de mediación con bootstrap
med_result <- mediate(
 fit_M, fit_Y,
 treat = "A", # exposición
 mediator = "M", # mediador
 boot = TRUE, # IC por bootstrap (recomendado)
 sims = 1000, # nº de réplicas bootstrap
 conf.level = 0.95
)
summary(med_result)
# Salida:
# - ACME (Average Causal Mediation Effect) = NIE
# - ADE (Average Direct Effect) = NDE
# - Total Effect = TCE
# - Prop. Mediated = PM
# forest plot interno; se amplía el margen superior para que el título
# no se solape con el marco del gráfico (plot.mediate lo dibuja pegado al box)
op <- par(mar = c(5, 6, 5, 2) + 0.1)
plot(med_result, main = "")
title(main = "Efectos de mediación: ACME (indirecto), ADE (directo) y total",
 line = 3, col.main = col_libro["navy"], font.main = 2, cex.main = 1.0)
par(op)
# Test de sensibilidad a confusión M-Y
# medsens() solo soporta lm+lm o lm+glm(gaussian); con dos glm(binomial) falla.
# Alternativas: (a) cmsens() de CMAverse; (b) E-value (Ficha 22)
# sens <- medsens(med_result, rho.by = 0.1)  # solo para lm+lm
# plot(sens)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Método moderno: paquete CMAverse (VanderWeele)
# ──────────────────────────────────────────────────────────────────────────────

library(CMAverse)
# Análisis completo con descomposición 4-way
cma_result <- cmest(
 data = dat,
 model = "rb", # Regression-based (Baron-Kenny tipo, pero causal)
 outcome = "Y",
 exposure = "A",
 mediator = "M",
 basec = covars, # confusores basales
 EMint = TRUE, # incluir interacción A×M
 yreg = "logistic", # tipo de modelo de Y
 mreg = list("logistic"), # tipo de modelo de M
 astar = 0, # valor de referencia de A
 a = 1, # valor de comparación de A
 mval = list(0), # nivel de M para CDE
 estimation = "imputation", # alternativa: "paramfunc"
 inference = "bootstrap",
 nboot = 1000
)
summary(cma_result)
# Reporta:
# - Te (Total Effect)
# - Pnde (Pure Natural Direct Effect)
# - Tnie (Total Natural Indirect Effect)
# - Cde (Controlled Direct Effect)
# - INTref, INTmed, PIE (descomposición 4-way)
# - pm (Proportion Mediated)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis con efectos no lineales: paquete regmedint
# ──────────────────────────────────────────────────────────────────────────────

library(regmedint)
# regmedint exige que c_cond sea NUMÉRICO: crear dummies para factores
dat$sexo_mujer  <- as.integer(dat$sexo == "mujer")
dat$educ_media  <- as.integer(dat$educ == "media")
dat$educ_alta   <- as.integer(dat$educ == "alta")
# Análisis paramétrico con interacción A×M
# c_cond: valores de referencia numéricos (media/moda de la muestra)
reg_med <- regmedint(
 data = dat,
 yvar = "Y",
 avar = "A",
 mvar = "M",
 cvar = c("edad", "sexo_mujer", "educ_media", "educ_alta"),
 c_cond = c(edad = 50, sexo_mujer = 0, educ_media = 1, educ_alta = 0),
 a0 = 0, a1 = 1, m_cde = 0,
 interaction = TRUE,
 yreg = "logistic",
 mreg = "logistic"
)
summary(reg_med)
# Reporta NDE, NIE, CDE con SE delta (más rápido que bootstrap)

# ──────────────────────────────────────────────────────────────────────────────
# ■ G-methods para confusor M-Y afectado por A
# ──────────────────────────────────────────────────────────────────────────────

# Cuando hay confusor L tal que A → L → M y L → Y
# (violación del supuesto 4), usar g-methods (Ficha 14)
# Añadir confusor L (afectado por A, afecta a M e Y) para ilustrar postc
set.seed(99)
dat$L <- rbinom(nrow(dat), 1, plogis(-0.5 + 0.8 * dat$A + 0.01 * (dat$edad - 50)))
library(CMAverse)
cma_gmethod <- cmest(
 data = dat,
 model = "gformula", # g-formula
 outcome = "Y",
 exposure = "A",
 mediator = "M",
 postc = "L", # confusor M-Y afectado por A
 basec = covars,
 EMint = TRUE, # obligatorio en cmest (interacción A×M)
 yreg = "logistic",
 mreg = list("logistic"),
 postcreg = list("logistic"), # modelo de L
 astar = 0, a = 1, mval = list(0),
 yval = "1",
 estimation = "imputation",
 inference = "bootstrap",
 nboot = 1000
)
summary(cma_gmethod)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis de sensibilidad a confusión no observada
# ──────────────────────────────────────────────────────────────────────────────

# E-value para mediación (VanderWeele & Ding 2017)
# Permite cuantificar la robustez del NIE a confusión M-Y no medida
# E-value para mediación (VanderWeele & Ding 2017)
# evalues.OR() falla con lava cargado; función cerrada con rare = FALSE
ev_vd <- function(x) if (x <= 1) ev_vd(1 / x) else x + sqrt(x * (x - 1))
e_value_or <- function(est, lo, hi, rare = FALSE) {
  to_rr <- function(x) if (rare) x else 2 * x / (1 + x)
  rr <- to_rr(est)
  rr_ic <- if (lo <= 1 && hi >= 1) 1 else to_rr(if (est > 1) lo else hi)
  c(point = ev_vd(rr), ic = ev_vd(rr_ic))
}
# Para el NIE (en escala OR)
NIE_or <- 1.20 # ejemplo
NIE_lo <- 1.08
NIE_hi <- 1.34
e_value_or(est = NIE_or, lo = NIE_lo, hi = NIE_hi, rare = FALSE)
# Reporta E-value puntual y para el límite del IC más cercano al nulo
