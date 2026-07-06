# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                      BLOQUE 3 — Diseños epidemiológicos                      
#                            Código R · Fichas 8–10                            
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
#   Ficha 08 · Estudios de casos y controles
#   Ficha 09 · Estudios de caso-cohorte
#   Ficha 10 · Estudios de cohorte
#
# Contenido: 3 fichas · 22 fragmentos de código · 462 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v90.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 08 · Estudios de casos y controles
# ║  Método: glm()/clogit() · OR crudo y ajustado; emparejamiento; E-value
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis no emparejado: tabla 2×2 y OR crudo
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(epitools) # oddsratio(), riskratio(), epitable()
library(survival) # clogit() para caso-control emparejado
library(broom) # tidy()
library(dplyr)
library(ggplot2)
# Datos de ejemplo (tabla 2x2 con valores conocidos)
# Estudio caso-control: tabaquismo y cáncer de pulmón
tab <- matrix(c(80, 120, # Casos: expuestos / no expuestos
 20, 180), # Controles: expuestos / no expuestos
 nrow = 2, byrow = FALSE,
 dimnames = list(Exposicion = c("Expuesto", "No expuesto"),
 Caso = c("Caso", "Control")))
print(tab)
# OR crudo con IC 95% (varios métodos)
epitools::oddsratio(tab, method = "wald") # Asintótico
epitools::oddsratio(tab, method = "midp") # Mid-p (recomendado, robusto a celdas pequeñas)
epitools::oddsratio(tab, method = "fisher") # Exacto (Fisher)
# Test de chi² o Fisher
chisq.test(tab)
fisher.test(tab)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis con regresión logística (no emparejado)
# ──────────────────────────────────────────────────────────────────────────────

# Datos individuales (ejemplo simulado)
set.seed(2024)
n_cas <- 250; n_ctrl <- 500
dat <- data.frame(
 D = c(rep(1, n_cas), rep(0, n_ctrl)),
 edad = c(rnorm(n_cas, 60, 10), rnorm(n_ctrl, 55, 10)),
 sexo = factor(sample(c("hombre", "mujer"), n_cas + n_ctrl, replace = TRUE)),
 tabaco = factor(sample(c("nunca", "ex", "actual"),
 n_cas + n_ctrl, replace = TRUE,
 prob = c(0.3, 0.3, 0.4))),
 alcohol = factor(sample(c("bajo", "moderado", "alto"),
 n_cas + n_ctrl, replace = TRUE)),
 af_familiares = rbinom(n_cas + n_ctrl, 1, 0.2)
)
# Ajustar modelo logístico
covars <- c("edad", "sexo", "tabaco", "alcohol", "af_familiares")
formula_lr <- as.formula(paste("D ~", paste(covars, collapse = " + ")))
fit_lr <- glm(formula_lr, data = dat, family = binomial)
summary(fit_lr)
# Tabla de OR ajustados con IC 95%
res_or <- broom::tidy(fit_lr, exponentiate = TRUE, conf.int = TRUE) %>%
 filter(term != "(Intercept)") %>%
 mutate(across(c(estimate, conf.low, conf.high), ~round(.x, 3)),
 p.value = format.pval(p.value, digits = 3, eps = 0.001))
print(res_or)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis emparejado: regresión logística condicional
# ──────────────────────────────────────────────────────────────────────────────

# Datos emparejados: cada caso con 2 controles (estrato)
set.seed(2024)
n_pairs <- 200
strata_id <- rep(1:n_pairs, each = 3)
dat_match <- data.frame(
 estrato = strata_id,
 D = rep(c(1, 0, 0), n_pairs), # 1 caso : 2 controles
 edad = rep(round(rnorm(n_pairs, 60, 10)), each = 3), # Matching exacto
 sexo = rep(sample(c("H", "M"), n_pairs, replace = TRUE), each = 3),
 tabaco = factor(sample(c("nunca", "actual"), 3 * n_pairs,
 replace = TRUE)),
 alcohol = rnorm(3 * n_pairs, 20, 10)
)
# IMPORTANTE: edad y sexo son variables de matching;
# NO se incluyen como covariables (colinealidad con strata)
# Modelo condicional con clogit() (paquete survival)
fit_clr <- clogit(D ~ tabaco + alcohol + strata(estrato), data = dat_match)
summary(fit_clr)
# OR ajustados con IC 95%
broom::tidy(fit_clr, exponentiate = TRUE, conf.int = TRUE)
# Test de Wald, LR y Score (igual que en Cox)
# salida estándar de summary() incluye los tres

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis estratificado (Mantel-Haenszel)
# ──────────────────────────────────────────────────────────────────────────────

# OR común estratificado (alternativa al ajuste por regresión)
library(epitools)
# Crear tabla 2x2 por estrato (ej. por edad categorizada)
dat_estr <- dat %>%
 mutate(grupo_edad = cut(edad, breaks = c(0, 50, 65, 100),
 labels = c("<50", "50-65", ">65"))) %>%
 mutate(tabaco_bin = ifelse(tabaco %in% c("ex", "actual"), 1, 0))
# Tabla estratificada
tab_estr <- xtabs(~ tabaco_bin + D + grupo_edad, data = dat_estr)
print(tab_estr)
# Test y OR de Mantel-Haenszel
mantelhaen.test(tab_estr)
# Test de homogeneidad de OR (Breslow-Day)
# Disponible en paquete vcdExtra
# vcdExtra::woolf_test(tab_estr)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación del emparejamiento
# ──────────────────────────────────────────────────────────────────────────────

# 1. Distribución de variables de matching antes/después
library(tableone)
# Comparación caso vs control en variables matched
CreateTableOne(vars = c("edad", "sexo"),
 strata = "D", data = dat_match,
 test = TRUE, smd = TRUE)
# SMD < 0.10 indica buen balance
# 2. Distribución de variables NO matched (deben diferir entre grupos)
CreateTableOne(vars = c("tabaco", "alcohol"),
 strata = "D", data = dat_match,
 test = TRUE)
# 3. Verificación gráfica: distribución de edad por estrato
ggplot(dat_match, aes(x = factor(D), y = edad, group = estrato)) +
 geom_line(alpha = 0.25, colour = col_libro["gris"]) +
 geom_point(alpha = 0.5, colour = col_libro["navy"]) +
 scale_x_discrete(labels = c("0" = "Control", "1" = "Caso")) +
 labs(title    = "Verificación del emparejamiento por edad",
 subtitle = "Cada línea conecta los miembros de un estrato caso-control",
 x        = "Condición",
 y        = "Edad (años)",
 caption  = "Un buen emparejamiento produce líneas aproximadamente horizontales") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos del modelo logístico
# ──────────────────────────────────────────────────────────────────────────────

# Diagnósticos estándar (igual que en Ficha 1)
library(performance)
performance::check_model(fit_lr)
performance::performance_hosmer(fit_lr) # Hosmer-Lemeshow
performance::check_collinearity(fit_lr) # VIF
# Curva ROC y AUC
library(pROC)
pred_prob <- predict(fit_lr, type = "response")
roc_obj <- pROC::roc(dat$D, pred_prob)
print(auc(roc_obj))
plot(roc_obj, print.auc = TRUE,
 col = col_libro["navy"], lwd = 2,
 main = "Curva ROC del modelo logístico",
 print.auc.col = col_libro["ochre"])
# Validación interna por bootstrap (rms)
library(rms)
fit_rms <- lrm(formula_lr, data = dat, x = TRUE, y = TRUE)
validate(fit_rms, B = 200)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis de sensibilidad
# ──────────────────────────────────────────────────────────────────────────────

# 1. E-value para confusión no medida (Ficha 23)
library(EValue)
# E-value de forma cerrada (VanderWeele & Ding, 2017).
# EValue::evalues.OR() falla con lava cargado; se usa función cerrada.
# Para OR con desenlace raro (rare = TRUE) el OR se trata como RR.
ev_vd <- function(x) if (x <= 1) ev_vd(1 / x) else x + sqrt(x * (x - 1))
e_value_or <- function(est, lo, hi, rare = TRUE) {
  to_rr <- function(x) if (rare) x else 2 * x / (1 + x)
  rr <- to_rr(est)
  rr_ic <- if (lo <= 1 && hi >= 1) 1 else to_rr(if (est > 1) lo else hi)
  c(point = ev_vd(rr), ic = ev_vd(rr_ic))
}
e_value_or(est = 3.45, lo = 2.20, hi = 5.41, rare = TRUE)
# rare = TRUE: el evento se considera raro y OR ≈ RR
# 2. Sensibilidad al sesgo de selección (selection bias)
# El paquete sensemakr permite cuantificarlo
# 3. Análisis estratificado por subgrupos para evaluar interacción
# Modelo con interacción
fit_int <- glm(D ~ tabaco * sexo + edad + alcohol,
 data = dat, family = binomial)
anova(fit_lr, fit_int, test = "Chisq") # Test LR de interacción


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 09 · Estudios de caso-cohorte
# ║  Método: survival::cch() · diseño de caso-cohorte de Prentice/Borgan
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival) # cch() — gold standard
library(dplyr)
library(broom)
library(ggplot2)
# Estructura mínima de datos esperada (ejemplo simulado)
set.seed(2024)
N <- 5000 # cohorte completa
m <- 750 # subcohorte (15% de N)
p_event <- 0.06 # 6% incidencia
dat <- data.frame(
 id = 1:N,
 time = rweibull(N, 1.4, 30),
 evento = rbinom(N, 1, p_event),
 edad = round(rnorm(N, 55, 10)),
 sexo = factor(sample(c("hombre", "mujer"), N, replace = TRUE)),
 centro = factor(sample(c("Centro D", "Asturias", "Granada", "Navarra"),
 N, replace = TRUE)),
 APS = rnorm(N, 0, 1), # Alcohol Proteomic Score (Z-score)
 tabaco = factor(sample(c("nunca", "ex", "actual"), N, replace = TRUE,
 prob = c(0.4, 0.3, 0.3))),
 IMC = rnorm(N, 26, 4)
)
# Selección de subcohorte aleatoria
sub_idx <- sample(N, m)
dat$subcohort <- as.numeric(seq_len(N) %in% sub_idx)
# En CC solo se miden APS, tabaco e IMC en sujetos de subcohorte O casos
dat_cc <- dat %>%
 mutate(en_estudio = (subcohort == 1) | (evento == 1)) %>%
 filter(en_estudio)
# Resumen de la estructura
table(dat_cc$evento, dat_cc$subcohort,
 dnn = c("Caso", "En subcohorte"))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis con cch() — método Self-Prentice
# ──────────────────────────────────────────────────────────────────────────────

# Definición del conjunto de covariables (UNA sola vez)
covars <- c("APS", "edad", "sexo", "tabaco", "IMC")
formula_cc <- as.formula(
 paste("Surv(time, evento) ~", paste(covars, collapse = " + "))
)
# Ajuste con cch() — Self-Prentice (recomendado por defecto)
fit_cc <- cch(
 formula = formula_cc,
 data = dat_cc,
 subcoh = ~ subcohort,
 id = ~ id,
 cohort.size = N, # tamaño de la cohorte completa
 method = "SelfPrentice"
)
summary(fit_cc)
# Extracción ordenada
res_cc <- data.frame(
 variable = names(coef(fit_cc)),
 HR = exp(coef(fit_cc)),
 SE = sqrt(diag(vcov(fit_cc))),
 HR_lo = exp(coef(fit_cc) - 1.96 * sqrt(diag(vcov(fit_cc)))),
 HR_hi = exp(coef(fit_cc) + 1.96 * sqrt(diag(vcov(fit_cc)))),
 z = coef(fit_cc) / sqrt(diag(vcov(fit_cc))),
 p_value = 2 * (1 - pnorm(abs(coef(fit_cc) / sqrt(diag(vcov(fit_cc))))))
)
print(res_cc, digits = 3)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Métodos alternativos de estimación
# ──────────────────────────────────────────────────────────────────────────────

# Comparación entre los cuatro métodos disponibles
fit_prent <- cch(formula_cc, data = dat_cc, subcoh = ~ subcohort,
 id = ~ id, cohort.size = N, method = "Prentice")
fit_sprent <- cch(formula_cc, data = dat_cc, subcoh = ~ subcohort,
 id = ~ id, cohort.size = N, method = "SelfPrentice")
fit_liny <- cch(formula_cc, data = dat_cc, subcoh = ~ subcohort,
 id = ~ id, cohort.size = N, method = "LinYing")
# I.Borgan: estimador para diseño ESTRATIFICADO; exige 'stratum' y un
# 'cohort.size' por estrato (NO estratificar por sexo, que está en el modelo).
csize_centro <- as.integer(table(dat$centro))    # cohorte completa por estrato
names(csize_centro) <- levels(dat$centro)         # nombres = niveles de centro
fit_borgan <- cch(formula_cc, data = dat_cc, subcoh = ~ subcohort,
 id = ~ id, stratum = ~ centro, cohort.size = csize_centro,
 method = "I.Borgan")
# Comparación de coeficientes y SE
methods <- list(Prentice = fit_prent, SelfPrentice = fit_sprent,
 LinYing = fit_liny, I.Borgan = fit_borgan)
sapply(methods, function(f) coef(f)["APS"])
sapply(methods, function(f) sqrt(diag(vcov(f)))["APS"])

# ──────────────────────────────────────────────────────────────────────────────
# ■ Subcohorte estratificada por centro y sexo
# ──────────────────────────────────────────────────────────────────────────────

# Cuando la subcohorte se selecciona estratificadamente
# (p. ej. estratificada por centro y sexo en una cohorte multicéntrica)
# Es necesario especificar los pesos de muestreo correctos
# Suponiendo que se han calculado los pesos de muestreo por estrato
dat_cc <- dat_cc %>%
 group_by(centro, sexo) %>%
 mutate(
 n_estrato = n(),
 n_sub_estrato = sum(subcohort == 1),
 peso_estrato = ifelse(subcohort == 1,
 n_estrato / n_sub_estrato,
 1)
 ) %>%
 ungroup()
# Cox ponderado equivalente (Barlow, 1994) implementado vía coxph
fit_barlow <- coxph(
 formula_cc,
 data = dat_cc,
 weights = dat_cc$peso_estrato,
 cluster = id,
 robust = TRUE
)
summary(fit_barlow)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis para múltiples desenlaces
# ──────────────────────────────────────────────────────────────────────────────

# Ventaja clave del CC: la MISMA subcohorte sirve para varios desenlaces
# (sin re-muestrear casos ni controles)
# Suponiendo varios desenlaces en el dataset:
# evento_1 = recurrencia, evento_2 = mortalidad cardiovascular, etc.
# (Datos ilustrativos) Simulamos tres desenlaces en la cohorte completa:
set.seed(11)
dat$evento_1 <- rbinom(nrow(dat), 1, 0.05)  # p. ej. recurrencia
dat$evento_2 <- rbinom(nrow(dat), 1, 0.04)  # p. ej. mortalidad cardiovascular
dat$evento_3 <- rbinom(nrow(dat), 1, 0.03)
desenlaces <- c("evento_1", "evento_2", "evento_3")
resultados <- list()
for (out in desenlaces) {
 # La MISMA subcohorte + los casos de ESE desenlace
 dat_out <- dat[dat$subcohort == 1 | dat[[out]] == 1, ]
 formula_i <- as.formula(
 paste("Surv(time,", out, ") ~", paste(covars, collapse = " + "))
 )
 fit_i <- cch(formula_i, data = dat_out, subcoh = ~ subcohort,
 id = ~ id, cohort.size = N, method = "SelfPrentice")
 co <- coef(fit_i); se <- sqrt(diag(vcov(fit_i)))   # cch no tiene método broom::tidy
 resultados[[out]] <- data.frame(
 variable = names(co), HR = exp(co),
 HR_lo = exp(co - 1.96 * se), HR_hi = exp(co + 1.96 * se)
 )
}
# Combinar resultados
do.call(rbind, lapply(names(resultados), function(out)
 cbind(desenlace = out, resultados[[out]])))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación de la subcohorte
# ──────────────────────────────────────────────────────────────────────────────

# 1. La subcohorte debe ser representativa de la cohorte basal
# Comparación de variables basales entre subcohorte y cohorte completa
library(tableone)
CreateTableOne(
 vars = c("edad", "sexo", "centro", "tabaco", "IMC"),
 strata = "subcohort",
 data = dat,
 test = TRUE, smd = TRUE
)
# SMD < 0.10 indica que la subcohorte es representativa
# 2. Tasa de eventos en subcohorte vs no subcohorte
dat %>%
 group_by(subcohort) %>%
 summarise(
 n = n(),
 eventos = sum(evento),
 tasa = sum(evento) / n
 )

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación del supuesto HP
# ──────────────────────────────────────────────────────────────────────────────

# Test de Schoenfeld adaptado al modelo CC
# (cox.zph no funciona directamente sobre cch; usar coxph con pesos)
fit_check <- coxph(
 formula_cc,
 data = dat_cc,
 weights = ifelse(dat_cc$subcohort == 1, 1/(m/N), 1),
 cluster = id,
 robust = TRUE
)
zph_test <- cox.zph(fit_check)
print(zph_test)
plot(zph_test, var = "APS", col = col_libro["rojo"],
 main = "Residuos de Schoenfeld — supuesto de riesgos proporcionales")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación con Cox sobre cohorte completa
# ──────────────────────────────────────────────────────────────────────────────

# Test de robustez: ajustar Cox sobre cohorte completa para variables
# ya medidas en todos (edad, sexo, tabaco básico)
# Si los HR son consistentes con los del CC → confianza en el diseño
covars_all <- c("edad", "sexo", "tabaco") # variables medidas en cohorte completa
formula_all <- as.formula(
 paste("Surv(time, evento) ~", paste(covars_all, collapse = " + "))
)
fit_full <- coxph(formula_all, data = dat) # cohorte completa
fit_cc_red <- cch(formula_all, data = dat_cc,
 subcoh = ~ subcohort, id = ~ id,
 cohort.size = N, method = "SelfPrentice")
# Comparación de HR
data.frame(
 variable = names(coef(fit_full)),
 HR_full = round(exp(coef(fit_full)), 3),
 HR_CC = round(exp(coef(fit_cc_red)), 3),
 ratio = round(exp(coef(fit_cc_red)) / exp(coef(fit_full)), 3)
)


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 10 · Estudios de cohorte
# ║  Método: glm()/coxph() · riesgo relativo, incidencia acumulada, tasas
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Datos y preparación
# ──────────────────────────────────────────────────────────────────────────────

# Carga de paquetes
library(survival) # Surv(), coxph()
library(epitools) # riskratio(), rateratio()
library(broom) # tidy()
library(dplyr)
library(ggplot2)
# Ejemplo simulado de cohorte
set.seed(2024)
N <- 5000
dat <- data.frame(
 id = 1:N,
 edad = rnorm(N, 55, 10),
 sexo = factor(sample(c("hombre", "mujer"), N, replace = TRUE)),
 expuesto = rbinom(N, 1, 0.30),
 fum_basal = factor(sample(c("nunca", "ex", "actual"),
 N, replace = TRUE, prob = c(0.4, 0.3, 0.3))),
 IMC = rnorm(N, 26, 4)
)
# Generar tiempos de seguimiento con hazard que depende de la exposición
dat$hazard_basal <- 0.02 * exp(0.6 * dat$expuesto + 0.03 * (dat$edad - 55))
dat$tiempo_ev <- rexp(N, dat$hazard_basal)
dat$tiempo_cens <- runif(N, 5, 15) # censura administrativa
dat$tiempo <- pmin(dat$tiempo_ev, dat$tiempo_cens)
dat$evento <- as.numeric(dat$tiempo_ev <= dat$tiempo_cens)
# Resumen descriptivo
dat %>%
 group_by(expuesto) %>%
 summarise(
 n = n(),
 eventos = sum(evento),
 pt_total = sum(tiempo),
 incidencia = mean(evento),
 tasa = sum(evento) / sum(tiempo) * 1000 # por 1000 personas-año
 )

# ──────────────────────────────────────────────────────────────────────────────
# ■ Tablas 2×2 y medidas crudas
# ──────────────────────────────────────────────────────────────────────────────

# Tabla 2x2 de incidencia (riesgo)
tab <- with(dat, table(expuesto, evento))
print(tab)
# RR e IC 95% con epitools
epitools::riskratio(tab, rev = "neither")
# RR e IC 95% manualmente (método Wald)
a <- tab[2, 2]; b <- tab[2, 1]; c <- tab[1, 2]; d <- tab[1, 1]
RR <- (a / (a + b)) / (c / (c + d))
se_logRR <- sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
ic_RR <- exp(log(RR) + c(-1, 1) * 1.96 * se_logRR)
data.frame(RR = round(RR, 3),
 IC_inf = round(ic_RR[1], 3),
 IC_sup = round(ic_RR[2], 3))
# Tasa e IRR con persona-tiempo
tasa_dat <- dat %>%
 group_by(expuesto) %>%
 summarise(eventos = sum(evento), pt = sum(tiempo))
tasa_dat
# IRR e IC 95% con epitools
epitools::rateratio(c(tasa_dat$eventos[2], tasa_dat$eventos[1]),
 c(tasa_dat$pt[2], tasa_dat$pt[1]))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Modelos ajustados
# ──────────────────────────────────────────────────────────────────────────────

# Conjunto de covariables (definido UNA vez)
covars <- c("expuesto", "edad", "sexo", "fum_basal", "IMC")
# 1) Logística para incidencia acumulada (OR ≈ RR si raro)
fit_lr <- glm(
 as.formula(paste("evento ~", paste(covars, collapse = " + "))),
 data = dat, family = binomial
)
summary(fit_lr)
broom::tidy(fit_lr, exponentiate = TRUE, conf.int = TRUE)
# 2) Poisson para tasas (con offset)
fit_poi <- glm(
 as.formula(paste("evento ~", paste(covars, collapse = " + "),
 "+ offset(log(tiempo))")),
 data = dat, family = poisson
)
summary(fit_poi)
broom::tidy(fit_poi, exponentiate = TRUE, conf.int = TRUE)
# 3) Cox para HR
fit_cox <- coxph(
 as.formula(paste("Surv(tiempo, evento) ~", paste(covars, collapse = " + "))),
 data = dat
)
summary(fit_cox)
broom::tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comparación entre las tres medidas
# ──────────────────────────────────────────────────────────────────────────────

# Extraer estimaciones de las tres aproximaciones
extract_or <- function(fit) {
 res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
 res[res$term == "expuesto", c("estimate", "conf.low", "conf.high")]
}
comparison <- rbind(
 cbind(modelo = "Logística (OR)", extract_or(fit_lr)),
 cbind(modelo = "Poisson (IRR)", extract_or(fit_poi)),
 cbind(modelo = "Cox (HR)", extract_or(fit_cox))
)
print(comparison, digits = 3)
# Cuando el evento es raro y HP se cumple, los tres son similares

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación del balance basal
# ──────────────────────────────────────────────────────────────────────────────

# Tabla 1: características basales por grupo de exposición
library(tableone)
CreateTableOne(
 vars = c("edad", "sexo", "fum_basal", "IMC"),
 strata = "expuesto",
 data = dat,
 test = TRUE, smd = TRUE
)
# SMD > 0.10 sugiere desequilibrio basal → potencial confusión

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis de pérdidas de seguimiento
# ──────────────────────────────────────────────────────────────────────────────

# 1. Tasa de pérdida por grupo
dat %>%
 group_by(expuesto) %>%
 summarise(
 n = n(),
 cens_admin = sum(evento == 0 & tiempo == tiempo_cens),
 cens_inform = sum(evento == 0 & tiempo < 14), # ejemplo
 eventos = sum(evento)
 )
# 2. Comparación de características entre completos y perdidos
dat$perdido <- as.numeric(dat$evento == 0 & dat$tiempo < 14)
CreateTableOne(
 vars = c("edad", "sexo", "fum_basal", "IMC", "expuesto"),
 strata = "perdido",
 data = dat,
 test = TRUE, smd = TRUE
)
# 3. Test de Schoenfeld para HP del Cox (Ficha 5)
cox.zph(fit_cox)
# 4. Análisis de sensibilidad: mejor/peor escenario
# Mejor escenario: todos los perdidos son no eventos
# Peor escenario: todos los perdidos son eventos
dat_mejor <- dat
dat_mejor$evento_alt <- dat$evento
dat_peor <- dat
dat_peor$evento_alt <- ifelse(dat$evento == 0 & dat$tiempo < 14, 1, dat$evento)
cox_mejor <- coxph(Surv(tiempo, evento) ~ expuesto, data = dat_mejor)
cox_peor <- coxph(Surv(tiempo, evento_alt) ~ expuesto, data = dat_peor)
data.frame(escenario = c("Original", "Mejor (pérdidas = no evento)", "Peor"),
 HR = c(exp(coef(fit_cox)["expuesto"]),
 exp(coef(cox_mejor)["expuesto"]),
 exp(coef(cox_peor)["expuesto"])))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Análisis estratificado por confusor
# ──────────────────────────────────────────────────────────────────────────────

# Análisis estratificado por edad (categorizada)
dat$grupo_edad <- cut(dat$edad, breaks = c(0, 50, 65, 100),
 labels = c("<50", "50-65", ">65"))
# OR/IRR estratificado
library(epitools)
# Tabla 2x2 por estrato
tab_estr <- xtabs(~ expuesto + evento + grupo_edad, data = dat)
print(tab_estr)
# Test de Mantel-Haenszel
mantelhaen.test(tab_estr)
# Comparar con modelo ajustado
fit_cox_ajust <- coxph(
 Surv(tiempo, evento) ~ expuesto + strata(grupo_edad),
 data = dat
)
summary(fit_cox_ajust)
