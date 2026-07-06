# ==============================================================================
#               MANUAL PRÁCTICO DE BIOESTADÍSTICA Y EPIDEMIOLOGÍA
#                 J. M. Huerta Castaño  ·  S. M. Colorado Yohar
# ------------------------------------------------------------------------------
#                  BLOQUE 8 — Multivariantes y clasificación                   
#                           Código R · Fichas 25–32                            
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
#   Ficha 25 · Análisis de componentes principales (PCA)
#   Ficha 26 · Análisis factorial exploratorio (EFA)
#   Ficha 27 · Análisis de conglomerados (cluster)
#   Ficha 28 · Análisis discriminante (LDA/QDA)
#   Ficha 29 · Análisis de correspondencias múltiples (MCA)
#   Ficha 30 · Análisis de clases latentes (LCA)
#   Ficha 31 · Random Forest
#   Ficha 32 · Regresión penalizada
#
# Contenido: 8 fichas · 56 fragmentos de código · 964 líneas.
# Codificación: UTF-8. Cada fragmento reproduce el código del libro;
# los comentarios en español son del manuscrito. Los títulos de apartado
# (marcados con «■») coinciden con los del texto del libro.
# Base: Libro_Metodos_Estadisticos_v93.docx.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 25 · Análisis de componentes principales (PCA)
# ║  Método: prcomp()/FactoMineR::PCA()
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ PCA básico con prcomp
# ──────────────────────────────────────────────────────────────────────────────

# ==============================================================================
# DATOS SINTÉTICOS — generados una sola vez antes de las fichas
# set.seed(2024) garantiza reproducibilidad completa del bloque.
# ==============================================================================
# ── Dataset A: FFQ sintético (Fichas 25–27) ───────────────────────────────────────
set.seed(2024)
N <- 500
F1 <- rnorm(N); F2 <- rnorm(N); F3 <- rnorm(N)
eps <- function() rnorm(N, 0, 0.8)
data_ffq <- data.frame(
  verduras=0.75*F1+eps(), fruta=0.70*F1+0.20*F3+eps(),
  pescado=0.65*F1+eps(), legumbres=0.60*F1+eps(),
  frutos_secos=0.55*F1+0.25*F3+eps(),
  carne_roja=0.72*F2+eps(), procesados=0.78*F2+eps(),
  dulces=0.65*F2+eps(), cereales=0.55*F2+eps(),
  lacteos=0.70*F3+eps())
data_ffq <- as.data.frame(lapply(data_ffq,
  function(x) pmax(0, x + abs(min(x)) + 0.5)))
data_ffq$edad <- round(rnorm(N,52,12))
data_ffq$sexo <- rbinom(N,1,0.5)
data_ffq$tabaquismo <- rbinom(N,1,0.25)
data_ffq$imc <- rnorm(N,27,4)
data_ffq$tiempo <- rexp(N,0.05)
data_ffq$evento <- rbinom(N,1,0.15)
food_vars <- c("verduras","fruta","pescado","carne_roja","procesados",
               "dulces","lacteos","cereales","legumbres","frutos_secos")
# ── Dataset B: psych::bfi para EFA (Ficha 26) ─────────────────────────────
library(psych)
bfi_clean <- na.omit(psych::bfi[, 1:25])
colnames(bfi_clean) <- paste0("item", 1:25)
data_cuestionario <- bfi_clean; items <- bfi_clean
set.seed(42)
idx_efa <- sample(nrow(bfi_clean), floor(nrow(bfi_clean)*0.5))
data_validacion <- bfi_clean[-idx_efa, ]
# ── Dataset C: Pacientes clínicos (Fichas 27–28) ───────────────────────────
set.seed(2024)
n_pac <- 300; g_true <- sample(1:3, n_pac, replace=TRUE)
data_pacientes <- data.frame(
  imc=rnorm(n_pac,c(22,28,35)[g_true],3),
  glucosa=rnorm(n_pac,c(90,110,140)[g_true],15),
  pcr=rlnorm(n_pac,c(.5,1.2,2)[g_true],.5),
  pa_sistolica=rnorm(n_pac,c(115,130,155)[g_true],12),
  colesterol=rnorm(n_pac,c(175,210,240)[g_true],25),
  trigliceridos=rlnorm(n_pac,c(4.5,5,5.5)[g_true],.4),
  hba1c=rnorm(n_pac,c(5.2,6.1,7.8)[g_true],.8),
  edad=rnorm(n_pac,c(45,58,65)[g_true],10),
  biomarcador1=rnorm(n_pac,c(1.2,2.5,4)[g_true],.8),
  biomarcador2=rnorm(n_pac,c(3,2,1)[g_true],.7),
  grupo=factor(g_true,labels=c("A","B","C")),
  sexo=rbinom(n_pac,1,.5),
  tiempo=rexp(n_pac,.03),
  evento=as.factor(rbinom(n_pac,1,.2)))
levels(data_pacientes$evento) <- c("no","si")
data_mixta <- data_pacientes
set.seed(42)
train_idx <- sample(nrow(data_pacientes), floor(nrow(data_pacientes)*0.7))
data_train <- data_pacientes[train_idx,]
data_test  <- data_pacientes[-train_idx,]
levels(data_train$evento) <- c("no","si")
levels(data_test$evento)  <- c("no","si")
predictors <- c("biomarcador1","biomarcador2","edad","imc")
# ── Dataset D: Variables categóricas para MCA (Ficha 29) ───────────────────
set.seed(2024); n_mca <- 400
data_estudio <- data.frame(
  tabaco=factor(sample(c("nunca","exfumador","activo"),n_mca,replace=TRUE,prob=c(.5,.3,.2))),
  diagnostico=factor(sample(c("sano","HTA","diabetes"),n_mca,replace=TRUE,prob=c(.5,.3,.2))),
  act_fisica=factor(sample(c("baja","moderada","alta"),n_mca,replace=TRUE)),
  grupo_edad=factor(sample(c("<40","40-60",">60"),n_mca,replace=TRUE)),
  nivel_educativo=factor(sample(c("bajo","medio","alto"),n_mca,replace=TRUE)),
  estado_civil=factor(sample(c("soltero","casado","otro"),n_mca,replace=TRUE)),
  mortalidad_5a=factor(rbinom(n_mca,1,.08),labels=c("vivo","fallecido")),
  edad_continua=rnorm(n_mca,52,14))
cat_vars <- c("tabaco","diagnostico","act_fisica",
              "grupo_edad","nivel_educativo","estado_civil")
# ── Dataset E: Cohorte con comorbilidades (Fichas 30–32) ──────────────────
set.seed(2024); n_lca <- 600
indicadores <- c("hipertension","diabetes","obesidad",
                  "epoc","depresion","artrosis","dislipidemia")
data_cohorte <- data.frame(
  hipertension=rbinom(n_lca,1,.35)+1L, diabetes=rbinom(n_lca,1,.15)+1L,
  obesidad=rbinom(n_lca,1,.25)+1L, epoc=rbinom(n_lca,1,.10)+1L,
  depresion=rbinom(n_lca,1,.20)+1L, artrosis=rbinom(n_lca,1,.30)+1L,
  dislipidemia=rbinom(n_lca,1,.40)+1L,
  edad=rnorm(n_lca,58,12), sexo=rbinom(n_lca,1,.5),
  educacion=sample(1:3,n_lca,replace=TRUE),
  pa_sistolica=rnorm(n_lca,128,18), tabaco=rbinom(n_lca,1,.25),
  colesterol=rnorm(n_lca,205,35), glucosa=rnorm(n_lca,102,25),
  act_fisica=rbinom(n_lca,1,.4), imc=rnorm(n_lca,27,4),
  tiempo=rexp(n_lca,.04),
  evento=as.factor(rbinom(n_lca,1,.18)))
levels(data_cohorte$evento) <- c("no","si")
predictores <- c("edad","pa_sistolica","tabaco","colesterol",
                  "imc","glucosa","sexo","act_fisica")
confusores_matrix <- model.matrix(~ edad + pa_sistolica + colesterol +
  imc + glucosa + sexo + act_fisica - 1, data=data_cohorte)
desenlace  <- as.integer(data_cohorte$evento == "si")
exposicion <- data_cohorte$tabaco
library(FactoMineR) # PCA con herramientas epidemiológicas
library(factoextra) # visualización
library(psych) # análisis paralelo, rotación
library(dplyr)
# Datos: cuestionario de frecuencia alimentaria (FFQ)
# N sujetos × P grupos de alimentos
# data_ffq: dataset sintético definido en el bloque de datos al inicio del bloque
# Variables de alimentos (excluir ID y covariables)
food_vars <- c("verduras", "fruta", "pescado", "carne_roja",
 "procesados", "dulces", "lacteos", "cereales",
 "legumbres", "frutos_secos")
X <- data_ffq[, food_vars]
# IMPORTANTE: estandarizar (scale = TRUE usa matriz de correlación)
pca <- prcomp(X, center = TRUE, scale. = TRUE)
# Resumen: desviación estándar, proporción de varianza, acumulada
summary(pca)
# Eigenvalues (= sdev^2)
eigenvalues <- pca$sdev^2
print(eigenvalues)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Selección del número de componentes
# ──────────────────────────────────────────────────────────────────────────────

# 1. CRITERIO DE KAISER (λ > 1)
sum(eigenvalues > 1)
# 2. SCREE PLOT
factoextra::fviz_eig(pca, addlabels = TRUE,
 barfill = unname(col_libro["azul_med"]),
 barcolor = unname(col_libro["navy"]),
 main = "Sedimentación (scree plot) del ACP",
 xlab = "Componente principal",
 ylab = "% de varianza explicada") +
 geom_hline(yintercept = 100/length(food_vars),
 linetype = "dashed", color = col_libro["rojo"]) +
 theme_libro()
# 3. VARIANZA ACUMULADA
cumsum(eigenvalues) / sum(eigenvalues) * 100
# 4. ANÁLISIS PARALELO DE HORN (el más riguroso)
psych::fa.parallel(X, fa = "pc", n.iter = 100,
 main = "Análisis paralelo de Horn")
# Retener componentes cuyo eigenvalue observado > eigenvalue aleatorio

# ──────────────────────────────────────────────────────────────────────────────
# ■ Rotación VARIMAX
# ──────────────────────────────────────────────────────────────────────────────

# Retener K componentes (ej: K = 3) y rotar VARIMAX
K <- 3
# Con psych::principal (incluye rotación)
pca_rotated <- psych::principal(
 X,
 nfactors = K,
 rotate = "varimax", # "none", "varimax", "promax", "oblimin"
 scores = TRUE
)
print(pca_rotated$loadings, cutoff = 0.4) # mostrar solo |loading| > 0.4
# Los scores rotados (para usar como exposición)
scores_rotated <- pca_rotated$scores
head(scores_rotated)
# Interpretar cada componente según sus loadings altos
# RC1: verduras, fruta, pescado → "patrón saludable"
# RC2: carne_roja, procesados, dulces → "patrón occidental"
# RC3: ...

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización: biplot y loadings
# ──────────────────────────────────────────────────────────────────────────────

# Biplot: scores + loadings simultáneos
factoextra::fviz_pca_biplot(
 pca,
 geom.ind = "point",
 col.ind = unname(col_libro["gris"]),
 col.var = unname(col_libro["navy"]),
 repel = TRUE,
 title = "Biplot del ACP — patrones dietéticos"
) +
 labs(x = "Primera componente principal (Dim 1)",
 y = "Segunda componente principal (Dim 2)") +
 theme_libro()
# Solo loadings (variables)
factoextra::fviz_pca_var(
 pca,
 col.var = "contrib", # color por contribución (gradiente semántico: bajo→alto)
 gradient.cols = unname(col_libro[c("gris", "ochre", "rojo")]),
 repel = TRUE,
 title = "Contribución de las variables a las componentes"
) +
 theme_libro()
# Heatmap de loadings
loadings_mat <- pca$rotation[, 1:K]
library(pheatmap)
pheatmap(loadings_mat,
 color = colorRampPalette(c(col_libro["navy"], "white", col_libro["ochre"]))(100),  # divergente: carga -/0/+
 display_numbers = TRUE,
 main = "Mapa de calor de las cargas de las componentes (ACP)",
 cluster_rows = FALSE, cluster_cols = FALSE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Uso de scores en análisis epidemiológico
# ──────────────────────────────────────────────────────────────────────────────

# Añadir los scores al dataframe principal
data_ffq$patron_saludable <- scores_rotated[, 1]
data_ffq$patron_occidental <- scores_rotated[, 2]
# Categorizar en cuantiles (terciles o cuartiles)
data_ffq$saludable_q <- cut(
 data_ffq$patron_saludable,
 breaks = quantile(data_ffq$patron_saludable, probs = 0:4/4),
 labels = c("Q1", "Q2", "Q3", "Q4"),
 include.lowest = TRUE
)
# Usar como exposición en modelo de Cox
library(survival)
fit <- coxph(
 Surv(tiempo, evento) ~ saludable_q + edad + sexo + tabaquismo + imc,
 data = data_ffq
)
summary(fit)
# HR por cuartil del patrón dietético
# O como variable continua (por incremento de 1 DT)
fit_cont <- coxph(
 Surv(tiempo, evento) ~ patron_saludable + edad + sexo + tabaquismo,
 data = data_ffq
)
summary(fit_cont)

# ──────────────────────────────────────────────────────────────────────────────
# ■ PCA con FactoMineR (más detalle epidemiológico)
# ──────────────────────────────────────────────────────────────────────────────

# PCA con variables suplementarias (no usadas en la construcción)
pca_fm <- PCA(
 data_ffq[, c(food_vars, "edad", "imc")],
 scale.unit = TRUE,
 quanti.sup = c(11, 12), # edad e imc como suplementarias
 ncp = 5, # número de dimensiones a retener
 graph = FALSE
)
# Contribuciones de las variables
pca_fm$var$contrib
# Calidad de representación (cos²)
pca_fm$var$cos2
# Coordenadas de individuos
pca_fm$ind$coord


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 26 · Análisis factorial exploratorio (EFA)
# ║  Método: psych::fa()
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Adecuación y exploración
# ──────────────────────────────────────────────────────────────────────────────

library(psych) # fa, KMO, cortest.bartlett, fa.parallel
library(GPArotation) # rotaciones
library(nFactors) # criterios para número de factores
library(dplyr)
# Datos: ítems de un cuestionario (ej: calidad de vida, síntomas)
# items <- data_cuestionario  # psych::bfi renombrado: item1…item25 (ver bloque de datos)
# 1. ADECUACIÓN DE LOS DATOS
# Test de esfericidad de Bartlett
cortest.bartlett(cor(items), n = nrow(items))
# Debe ser significativo (p < 0.05)
# KMO (Kaiser-Meyer-Olkin)
KMO(items)
# Overall MSA debe ser > 0.6 (idealmente > 0.8)
# Revisar MSA individual de cada ítem; eliminar ítems con MSA < 0.5

# ──────────────────────────────────────────────────────────────────────────────
# ■ Número de factores
# ──────────────────────────────────────────────────────────────────────────────

# 2. NÚMERO DE FACTORES
# Análisis paralelo de Horn (recomendado)
fa.parallel(items, fa = "fa", n.iter = 100,
 fm = "minres",
 main = "Análisis paralelo de Horn")
# Sugiere el número de factores donde el eigenvalue observado
# supera al de datos aleatorios
# MAP de Velicer + otros criterios simultáneamente
VSS(items, n = 8, rotate = "oblimin", fm = "minres")
# Reporta MAP, BIC, complejidad, etc.
# nFactors: múltiples criterios
library(nFactors)
ev <- eigen(cor(items))
# parallel() fue eliminada de nFactors; análisis paralelo ya realizado con fa.parallel()
nS <- nScree(x = ev$values)  # sin aparallel: usa criterios de scree observados
plotnScree(nS)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Extracción y rotación
# ──────────────────────────────────────────────────────────────────────────────

# 3. EXTRACCIÓN con K factores (ej: K = 3)
K <- 3
# Máxima verosimilitud + rotación oblicua (recomendado)
fa_result <- fa(
 items,
 nfactors = K,
 fm = "ml", # "ml", "minres", "pa" (principal axis)
 rotate = "oblimin", # rotación OBLICUA por defecto
 scores = "regression" # método de cálculo de scores
)
# Imprimir cargas (solo |loading| > 0.3)
print(fa_result$loadings, cutoff = 0.3, sort = TRUE)
# Comunalidades
fa_result$communality
# Unicidades
fa_result$uniquenesses
# Correlación entre factores (matriz Phi, si oblicua)
fa_result$Phi

# ──────────────────────────────────────────────────────────────────────────────
# ■ Bondad de ajuste (ML-EFA)
# ──────────────────────────────────────────────────────────────────────────────

# 4. ÍNDICES DE AJUSTE (solo con fm = "ml")
fa_result$STATISTIC # chi-cuadrado
fa_result$PVAL # p-valor del test
fa_result$RMSEA # < 0.06 buen ajuste
fa_result$TLI # > 0.95 buen ajuste
fa_result$BIC # para comparar modelos con distinto K
# Comparar modelos con distinto número de factores
fa_2 <- fa(items, nfactors = 2, fm = "ml", rotate = "oblimin")
fa_3 <- fa(items, nfactors = 3, fm = "ml", rotate = "oblimin")
fa_4 <- fa(items, nfactors = 4, fm = "ml", rotate = "oblimin")
# Comparar BIC
data.frame(
 K = c(2, 3, 4),
 BIC = c(fa_2$BIC, fa_3$BIC, fa_4$BIC),
 RMSEA = c(fa_2$RMSEA[1], fa_3$RMSEA[1], fa_4$RMSEA[1])
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización
# ──────────────────────────────────────────────────────────────────────────────

# Diagrama de la estructura factorial
fa.diagram(fa_result, cut = 0.3, simple = TRUE,
 main = "Estructura factorial")
# Heatmap de cargas
library(pheatmap)
loadings_mat <- unclass(fa_result$loadings)
pheatmap(loadings_mat,
 color = colorRampPalette(c(col_libro["navy"], "white", col_libro["ochre"]))(100),  # divergente: carga -/0/+
 display_numbers = TRUE, main = "Mapa de calor de las cargas factoriales (EFA)",
 cluster_rows = FALSE, cluster_cols = FALSE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Confirmación con CFA (lavaan)
# ──────────────────────────────────────────────────────────────────────────────

library(lavaan)
# Tras EFA exploratorio, confirmar la estructura con CFA
# (idealmente en una muestra INDEPENDIENTE)
modelo_cfa <- '
 factor1 =~ item1 + item2 + item5 + item6
 factor2 =~ item3 + item4 + item7 + item8
 factor3 =~ item9 + item10 + item11 + item12
'
fit_cfa <- cfa(modelo_cfa, data = data_validacion,
 estimator = "MLR") # robusto a no-normalidad
summary(fit_cfa, fit.measures = TRUE, standardized = TRUE)
# Reporta: CFI, TLI, RMSEA, SRMR
# Buen ajuste: CFI > 0.95, RMSEA < 0.06, SRMR < 0.08

# ──────────────────────────────────────────────────────────────────────────────
# ■ Fiabilidad de las escalas derivadas
# ──────────────────────────────────────────────────────────────────────────────

# Alpha de Cronbach y omega para cada factor
psych::alpha(items[, c("item1", "item2", "item5", "item6")])
# raw_alpha > 0.7 aceptable
# Omega (más robusto que alpha)
psych::omega(items[, c("item1", "item2", "item5", "item6")],
 nfactors = 1)
# omega_total y omega_hierarchical


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 27 · Análisis de conglomerados (cluster)
# ║  Método: cluster/factoextra/NbClust
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Preparación y tendencia al clustering
# ──────────────────────────────────────────────────────────────────────────────

library(cluster) # pam, silhouette, daisy
library(factoextra) # visualización, fviz_*
library(NbClust) # consenso de criterios
library(fpc) # clusterboot, validación
library(dplyr)
# Datos: variables clínicas de pacientes
clin_vars <- c("imc", "glucosa", "pcr", "pa_sistolica",
 "colesterol", "trigliceridos", "hba1c", "edad")
X <- data_pacientes[, clin_vars]
# IMPORTANTE: estandarizar (z-score)
X_std <- scale(X)
# 1. TENDENCIA AL CLUSTERING (estadístico de Hopkins)
factoextra::get_clust_tendency(X_std, n = 50, graph = FALSE)
# hopkins_stat cercano a 1 → estructura de clusters
# cercano a 0.5 → datos aleatorios (no clusterizar)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Número de clusters
# ──────────────────────────────────────────────────────────────────────────────

# 2. NÚMERO DE CLUSTERS (varios criterios)
# Método del codo (WSS)
factoextra::fviz_nbclust(X_std, kmeans, method = "wss") +
 labs(title = "Método del codo (within sum of squares)",
 x = "Número de clusters (k)",
 y = "Suma de cuadrados intra-cluster") +
 theme_libro()
# Coeficiente de silueta
factoextra::fviz_nbclust(X_std, kmeans, method = "silhouette") +
 labs(title = "Coeficiente de silueta promedio",
 x = "Número de clusters (k)",
 y = "Anchura media de silueta") +
 theme_libro()
# Gap statistic (Tibshirani)
set.seed(123)
gap_stat <- cluster::clusGap(X_std, FUN = kmeans, nstart = 25,
 K.max = 10, B = 100)
factoextra::fviz_gap_stat(gap_stat) +
 labs(title = "Estadístico gap (Tibshirani et al., 2001)") +
 theme_libro()
# NbClust: consenso de ~30 índices
nb <- NbClust(X_std, distance = "euclidean",
 min.nc = 2, max.nc = 10,
 method = "kmeans")
# Reporta el k propuesto por mayoría de índices

# ──────────────────────────────────────────────────────────────────────────────
# ■ k-means
# ──────────────────────────────────────────────────────────────────────────────

# 3. K-MEANS con k = 3
set.seed(123)
km <- kmeans(X_std, centers = 3,
 nstart = 25, # 25 inicios aleatorios (CRÍTICO)
 iter.max = 100,
 algorithm = "Hartigan-Wong")
# Asignación de clusters
data_pacientes$cluster_km <- factor(km$cluster)
# Tamaño de cada cluster
table(km$cluster)
# WSS por cluster
km$withinss
km$tot.withinss
# Visualización (proyección PCA)
factoextra::fviz_cluster(km, data = X_std,
 palette = unname(col_libro[c("navy", "rojo", "ochre")]),
 ellipse.type = "convex",
 main = "Conglomerados k-means en el espacio de componentes",
 ggtheme = theme_libro())

# ──────────────────────────────────────────────────────────────────────────────
# ■ k-medoids (PAM) — robusto a outliers
# ──────────────────────────────────────────────────────────────────────────────

# 4. PAM con distancia euclídea (o Gower para datos mixtos)
pam_result <- cluster::pam(X_std, k = 3, metric = "euclidean")
# Para datos MIXTOS (continuas + categóricas): distancia de Gower
gower_dist <- cluster::daisy(data_mixta, metric = "gower")
pam_gower <- cluster::pam(gower_dist, k = 3, diss = TRUE)
# Medoides (observaciones representativas)
pam_result$medoids
# Silueta
factoextra::fviz_silhouette(silhouette(pam_result)) +
 labs(title = "Diagrama de silueta por cluster (PAM)") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Clustering jerárquico
# ──────────────────────────────────────────────────────────────────────────────

# 5. CLUSTERING JERÁRQUICO (Ward)
dist_mat <- dist(X_std, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
# Dendrograma
dend_p <- fviz_dend(hc, k = 3,
 cex = 0.6,
 k_colors = unname(col_libro[c("navy", "rojo", "ochre")]),
 main = "Dendrograma jerárquico (enlace de Ward)",
 rect = TRUE)
dend_p$layers[[1]]$aes_params$linewidth <- 0.3
print(dend_p)
# Cortar el dendrograma en 3 clusters
clusters_hc <- cutree(hc, k = 3)
data_pacientes$cluster_hc <- factor(clusters_hc)
# Correlación cofenética (fidelidad del dendrograma)
coph <- cophenetic(hc)
cor(dist_mat, coph)
# > 0.75 indica buena representación

# ──────────────────────────────────────────────────────────────────────────────
# ■ Clustering basado en modelos (GMM)
# ──────────────────────────────────────────────────────────────────────────────

library(mclust)
# Gaussian Mixture Model: selecciona número de clusters Y forma por BIC
gmm <- Mclust(X_std)
summary(gmm)
# Reporta: modelo óptimo (forma de covarianza), número de clusters (G)
# BIC para distintos modelos y G
plot(gmm, what = "BIC", main = "Selección de modelo y número de clusters (BIC)")
# Asignación probabilística (soft clustering)
gmm$classification # cluster más probable
gmm$z # probabilidades de pertenencia

# ──────────────────────────────────────────────────────────────────────────────
# ■ Validación de estabilidad
# ──────────────────────────────────────────────────────────────────────────────

# 6. ESTABILIDAD por bootstrap (Hennig)
library(fpc)
set.seed(123)
boot <- clusterboot(X_std,
 B = 100,
 clustermethod = kmeansCBI,
 k = 3,
 seed = 123)
# Jaccard medio por cluster
boot$bootmean
# > 0.75 estable; 0.6-0.75 dudoso; < 0.6 inestable
# Validación EXTERNA: asociación con desenlace no usado
# (ej: mortalidad a 5 años por cluster)
library(survival)
fit <- coxph(Surv(tiempo, as.integer(evento == "si")) ~ cluster_km +
 edad + sexo, data = data_pacientes)
summary(fit)
# Si los clusters predicen el desenlace → relevancia clínica


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 28 · Análisis discriminante (LDA/QDA)
# ║  Método: MASS::lda()/qda()
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ LDA básico
# ──────────────────────────────────────────────────────────────────────────────

library(MASS) # lda, qda
library(klaR) # partimat, visualización
library(caret) # validación cruzada
library(biotools) # boxM (test de Box)
library(dplyr)
# Datos: clasificar pacientes en grupos diagnósticos
# Predictores continuos + variable grupo
predictors <- c("biomarcador1", "biomarcador2", "edad", "imc")
# IMPORTANTE: estandarizar no es estrictamente necesario para LDA
# (es invariante a escala), pero ayuda a interpretar las cargas
# 1. AJUSTAR LDA
lda_fit <- lda(grupo ~ biomarcador1 + biomarcador2 + edad + imc,
 data = data_pacientes,
 prior = c(1/3, 1/3, 1/3)) # priors (o proporcionales)
print(lda_fit)
# Reporta:
# - Probabilidades a priori
# - Medias por grupo
# - Coeficientes de las funciones discriminantes (LD1, LD2, ...)
# - Proporción de traza explicada por cada discriminante

# ──────────────────────────────────────────────────────────────────────────────
# ■ Verificación de supuestos
# ──────────────────────────────────────────────────────────────────────────────

# Normalidad multivariante por grupo
library(MVN)
for (g in unique(data_pacientes$grupo)) {
 cat("Grupo:", g, "\n")
 subset_g <- data_pacientes[data_pacientes$grupo == g, predictors]
# Sustituido por mvn() con subset= (ver código corregido en el script)
}
# Homogeneidad de covarianzas (test M de Box)
# RECORDAR: muy sensible a n; interpretar con cautela
boxM_result <- biotools::boxM(
 data_pacientes[, predictors],
 data_pacientes$grupo
)
print(boxM_result)
# Si se rechaza fuertemente Y hay n grande → considerar QDA

# ──────────────────────────────────────────────────────────────────────────────
# ■ QDA y comparación
# ──────────────────────────────────────────────────────────────────────────────

# 2. AJUSTAR QDA
qda_fit <- qda(grupo ~ biomarcador1 + biomarcador2 + edad + imc,
 data = data_pacientes)
# Comparar LDA vs QDA por validación cruzada
ctrl <- caret::trainControl(method = "cv", number = 10)
lda_cv <- caret::train(grupo ~ biomarcador1 + biomarcador2 + edad + imc,
 data = data_pacientes, method = "lda",
 trControl = ctrl)
qda_cv <- caret::train(grupo ~ biomarcador1 + biomarcador2 + edad + imc,
 data = data_pacientes, method = "qda",
 trControl = ctrl)
# Comparar accuracy de CV
resamps <- resamples(list(LDA = lda_cv, QDA = qda_cv))
summary(resamps)
# El modelo con mayor accuracy de CV (y menor gap train-CV) es preferible

# ──────────────────────────────────────────────────────────────────────────────
# ■ Predicción y matriz de confusión
# ──────────────────────────────────────────────────────────────────────────────

# 3. PREDICCIÓN
pred <- predict(lda_fit, newdata = data_test)
# pred$class : clase predicha
# pred$posterior: probabilidades posteriores por grupo
# pred$x : scores en las funciones discriminantes
# Matriz de confusión
conf <- caret::confusionMatrix(pred$class, data_test$grupo)
print(conf)
# Reporta: accuracy, kappa, sensibilidad/especificidad por clase
# Validación cruzada leave-one-out (incluida en lda con CV=TRUE)
lda_loo <- lda(grupo ~ biomarcador1 + biomarcador2 + edad + imc,
 data = data_pacientes, CV = TRUE)
table(Real = data_pacientes$grupo, Predicho = lda_loo$class)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización
# ──────────────────────────────────────────────────────────────────────────────

# Scores en las funciones discriminantes
lda_scores <- predict(lda_fit)$x
plot_data <- data.frame(
 LD1 = lda_scores[, 1],
 LD2 = lda_scores[, 2],
 grupo = data_pacientes$grupo
)
library(ggplot2)
ggplot(plot_data, aes(x = LD1, y = LD2, color = grupo)) +
 geom_point(size = 2, alpha = 0.6) +
 stat_ellipse(level = 0.95) +
 scale_color_manual(values = unname(col_libro[c("navy", "rojo", "ochre")])) +
 labs(title    = "Proyección sobre los ejes discriminantes de Fisher",
 subtitle = "Elipses de confianza al 95% por grupo",
 x        = "Primera función discriminante (LD1)",
 y        = "Segunda función discriminante (LD2)",
 color    = "Grupo") +
 theme_libro()
# Particiones de decisión (klaR)
klaR::partimat(grupo ~ biomarcador1 + biomarcador2,
 data = data_pacientes, method = "lda")

# ──────────────────────────────────────────────────────────────────────────────
# ■ RDA (regularizado, intermedio)
# ──────────────────────────────────────────────────────────────────────────────

library(klaR)
# RDA: compromiso regularizado entre LDA y QDA (Friedman 1989)
rda_fit <- klaR::rda(grupo ~.,
 data = data_pacientes[, c("grupo", predictors)])
# Estima gamma (hacia identidad) y lambda (LDA <-> QDA) por CV
print(rda_fit$regularization)
# lambda = 0 → QDA; lambda = 1 → LDA; intermedio = regularizado


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 29 · Análisis de correspondencias múltiples (MCA)
# ║  Método: FactoMineR::MCA()
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ MCA con FactoMineR
# ──────────────────────────────────────────────────────────────────────────────

library(FactoMineR) # MCA, HCPC
library(factoextra) # visualización
library(dplyr)
# Datos: variables categóricas (TODAS deben ser factores)
cat_vars <- c("tabaco", "diagnostico", "act_fisica",
 "grupo_edad", "nivel_educativo", "estado_civil")
# Asegurar que son factores
data_cat <- data_estudio %>%
 mutate(across(all_of(cat_vars), as.factor))
# 1. AJUSTAR MCA
mca <- MCA(
 data_cat[, cat_vars],
 ncp = 5, # dimensiones a retener
 graph = FALSE
)
# Resumen: eigenvalues e inercia
mca$eig
# columnas: eigenvalue, % de varianza, % acumulado
# RECORDAR: estos % crudos sobreestiman; usar correcciones

# ──────────────────────────────────────────────────────────────────────────────
# ■ Inercias corregidas (Benzécri / Greenacre)
# ──────────────────────────────────────────────────────────────────────────────

# Corrección de Benzécri y Greenacre
# (las inercias crudas de MCA son engañosas)
# Número de variables y eigenvalues
J <- length(cat_vars)
eig <- mca$eig[, 1]
# Corrección de Benzécri: solo eigenvalues > 1/J
benzecri <- ifelse(eig > 1/J,
 (J/(J-1) * (eig - 1/J))^2,
 0)
pct_benzecri <- benzecri / sum(benzecri) * 100
print(round(pct_benzecri[1:5], 1))
# Greenacre es más conservador (implementado en ca::mjca)
library(ca)
mjca_fit <- ca::mjca(data_cat[, cat_vars], lambda = "adjusted")
summary(mjca_fit)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización del mapa factorial
# ──────────────────────────────────────────────────────────────────────────────

# Mapa de categorías
factoextra::fviz_mca_var(
 mca,
 col.var = "contrib", # color por contribución (gradiente semántico)
 gradient.cols = unname(col_libro[c("gris", "ochre", "rojo")]),
 repel = TRUE,
 title = "Mapa factorial de categorías (MCA)"
) +
 theme_libro()
# Mapa de individuos coloreados por una variable
factoextra::fviz_mca_ind(
 mca,
 habillage = "diagnostico", # colorear por diagnóstico
 palette = unname(col_libro[c("verde", "rojo", "ochre")]),
 addEllipses = TRUE,
 repel = FALSE,
 title = "Mapa de individuos por diagnóstico (MCA)"
) +
 theme_libro()
# Biplot conjunto (individuos + categorías)
factoextra::fviz_mca_biplot(
 mca,
 repel = TRUE,
 title = "Biplot conjunto de individuos y categorías (MCA)",
 ggtheme = theme_libro()
)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Diagnósticos: contribuciones y cos²
# ──────────────────────────────────────────────────────────────────────────────

# Contribuciones de las categorías a cada dimensión
mca$var$contrib[, 1:2] # contribuciones a Dim1 y Dim2
# Visualizar las categorías que MÁS contribuyen
factoextra::fviz_contrib(mca, choice = "var", axes = 1, top = 15) +
 labs(title = "Contribución de categorías a la primera dimensión") +
 theme_libro()
# La línea roja marca la contribución promedio esperada
# Calidad de representación (cos²)
mca$var$cos2[, 1:2]
factoextra::fviz_cos2(mca, choice = "var", axes = 1:2, top = 15) +
 labs(title = "Calidad de representación (cos²) en el plano factorial") +
 theme_libro()
# Descripción automática de las dimensiones
dimdesc(mca, axes = 1:2)
# Identifica qué variables y categorías caracterizan cada dimensión

# ──────────────────────────────────────────────────────────────────────────────
# ■ Variables suplementarias
# ──────────────────────────────────────────────────────────────────────────────

# Proyectar variables suplementarias (no influyen en los ejes)
mca_sup <- MCA(
 data_cat,
 quali.sup = which(names(data_cat) == "mortalidad_5a"), # suplementaria
 quanti.sup = which(names(data_cat) == "edad_continua"),
 ncp = 5,
 graph = FALSE
)
# Permite relacionar los ejes con un desenlace sin sesgarlos

# ──────────────────────────────────────────────────────────────────────────────
# ■ MCA + clustering (HCPC)
# ──────────────────────────────────────────────────────────────────────────────

# Clustering jerárquico sobre las coordenadas del MCA
# (estrategia potente para identificar perfiles cualitativos)
hcpc <- HCPC(mca, nb.clust = -1, graph = FALSE)
# nb.clust = -1: número óptimo automático
# Perfiles de los clusters (categorías sobre/infra-representadas)
hcpc$desc.var$category
# v.test > 2 o < -2 indica categorías características del cluster
# Visualizar los clusters en el mapa factorial
factoextra::fviz_cluster(hcpc, geom = "point",
 palette = unname(col_libro[c("navy", "rojo", "verde", "ochre", "teal", "azul_med", "gris")])[seq_len(length(unique(hcpc$data.clust$clust)))],
 main = "Clasificación jerárquica sobre componentes principales (HCPC)",
 ggtheme = theme_libro())
# Individuos representativos de cada cluster (paragons)
hcpc$desc.ind$para

# ──────────────────────────────────────────────────────────────────────────────
# ■ FAMD para datos mixtos
# ──────────────────────────────────────────────────────────────────────────────

# Si hay variables continuas Y categóricas: FAMD
famd <- FactoMineR::FAMD(
 data_mixta, # mezcla de numéricas y factores
 ncp = 5,
 graph = FALSE
)
# Combina PCA (continuas) y MCA (categóricas) en un marco unificado
factoextra::fviz_famd_var(famd) +
 labs(title = "Contribución de variables al análisis factorial mixto (FAMD)") +
 theme_libro()


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 30 · Análisis de clases latentes (LCA)
# ║  Método: poLCA
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ LCA básico con poLCA
# ──────────────────────────────────────────────────────────────────────────────

library(poLCA)
library(dplyr)
library(ggplot2)
# Datos: indicadores categóricos (poLCA requiere factores con niveles 1,2,...)
# Ejemplo: comorbilidades binarias (1 = No, 2 = Sí)
indicadores <- c("hipertension", "diabetes", "obesidad",
 "epoc", "depresion", "artrosis", "dislipidemia")
# Convertir a formato poLCA (enteros desde 1)
data_lca <- data_cohorte %>%
 mutate(across(all_of(indicadores), ~ as.integer(as.factor(.))))
# Fórmula: indicadores a la izquierda, 1 a la derecha (sin covariables)
f <- cbind(hipertension, diabetes, obesidad, epoc,
 depresion, artrosis, dislipidemia) ~ 1
# Ajustar modelos de 1 a 6 clases.
# calc.se = FALSE: las SE no son necesarias para comparar K por BIC/AIC/G2,
# y evita la inestabilidad numérica de poLCA.se() ante soluciones de frontera
# (probabilidades de respuesta en 0/1 dentro de alguna clase), causa habitual
# de errores de dimensión en modelos con clases pequeñas o indicadores dispersos.
# tryCatch: un K inestable no debe interrumpir la comparación de los demás.
set.seed(123)
modelos <- list()
for (k in 1:6) {
 modelos[[k]] <- tryCatch(
 poLCA(f, data = data_lca, nclass = k,
 nrep = 20, # 20 inicios (evitar máximos locales)
 maxiter = 5000,
 calc.se = FALSE,
 verbose = FALSE),
 error = function(e) {
 message(sprintf("poLCA falló para K = %d: %s", k, conditionMessage(e)))
 NULL
 })
}
# Ks con ajuste válido (excluye los que hayan fallado, si los hubiera)
k_ok <- which(!sapply(modelos, is.null))
if (length(k_ok) < 6) {
 warning(sprintf("Modelos no ajustados: K = %s",
 paste(setdiff(1:6, k_ok), collapse = ", ")))
}

# ──────────────────────────────────────────────────────────────────────────────
# ■ Selección del número de clases
# ──────────────────────────────────────────────────────────────────────────────

# Comparar criterios de información (solo K con ajuste válido)
comparacion <- data.frame(
 K = k_ok,
 logLik = sapply(modelos[k_ok], function(m) m$llik),
 BIC = sapply(modelos[k_ok], function(m) m$bic),
 AIC = sapply(modelos[k_ok], function(m) m$aic),
 n_param = sapply(modelos[k_ok], function(m) m$npar),
 G2 = sapply(modelos[k_ok], function(m) m$Gsq) # bondad de ajuste
)
print(comparacion)
# Graficar BIC vs K (buscar el mínimo o el codo)
ggplot(comparacion, aes(x = K, y = BIC)) +
 geom_line(color = col_libro["navy"], linewidth = 1) +
 geom_point(color = col_libro["navy"], size = 3) +
 labs(title    = "Selección del número de clases latentes",
 subtitle = "BIC frente al número de clases (se prefiere el mínimo)",
 x        = "Número de clases (K)",
 y        = "BIC",
 caption  = "Criterio de información bayesiano; penaliza la complejidad del modelo") +
 theme_libro()
# Verificar replicación del máximo (los nrep inicios)
# poLCA reporta si la mejor solución se replicó
# Si no, aumentar nrep

# ──────────────────────────────────────────────────────────────────────────────
# ■ Examinar la solución elegida
# ──────────────────────────────────────────────────────────────────────────────

# Supongamos K = 3 según BIC + interpretabilidad.
# Reajuste único con calc.se = TRUE (por defecto): los modelos del bucle de
# selección se ajustaron con calc.se = FALSE, por lo que el modelo elegido
# se recalcula aquí para disponer de errores estándar en la interpretación.
mod3 <- poLCA(f, data = data_lca, nclass = 3,
 nrep = 20, maxiter = 5000, verbose = FALSE)
# Probabilidades condicionales (perfil de cada clase)
mod3$probs
# Lista por indicador: P(respuesta | clase)
# Prevalencia de cada clase
mod3$P
# proporción de la muestra en cada clase
# Entropía (calcular manualmente)
entropy <- function(p) sum(-p * log(p + 1e-12))
error_prior <- entropy(mod3$P)
error_post <- mean(apply(mod3$posterior, 1, entropy))
entropia_R2 <- 1 - error_post / error_prior
cat("Entropía:", round(entropia_R2, 3), "\n")
# > 0.80 indica buena separación de clases
# Asignación modal
clase_modal <- mod3$predclass
table(clase_modal)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Visualización de perfiles
# ──────────────────────────────────────────────────────────────────────────────

# Preparar datos para graficar las probabilidades condicionales
library(tidyr)
probs_df <- lapply(seq_along(mod3$probs), function(j) {
 p_si <- mod3$probs[[j]][, 2] # P(Sí | clase)
 data.frame(indicador = indicadores[j],
 clase = paste0("Clase ", 1:length(p_si)),
 prob = p_si)
}) %>% bind_rows()
# Gráfico de perfiles (líneas por clase)
ggplot(probs_df, aes(x = indicador, y = prob,
 color = clase, group = clase)) +
 geom_line(linewidth = 1) +
 geom_point(size = 2.5) +
 scale_color_manual(values = unname(col_libro[c("rojo", "navy", "verde")])) +
 labs(title    = "Perfiles de probabilidad de las clases latentes",
 subtitle = "Probabilidad de presencia de cada indicador por clase",
 x        = "Indicador",
 y        = "P(presencia | clase)",
 color    = "Clase latente") +
 theme_libro() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ──────────────────────────────────────────────────────────────────────────────
# ■ Covariables (predictores de clase)
# ──────────────────────────────────────────────────────────────────────────────

# LCA con covariables: modelo de regresión multinomial sobre la clase
# Fórmula con predictores a la derecha
f_cov <- cbind(hipertension, diabetes, obesidad, epoc,
 depresion, artrosis, dislipidemia) ~ edad + sexo + educacion
mod_cov <- poLCA(f_cov, data = data_lca, nclass = 3,
 nrep = 20, maxiter = 5000)
# Coeficientes: efecto de las covariables sobre la pertenencia a clase
mod_cov$coeff
# (referencia = clase 1; OR = exp(coeff))
# NOTA: para desenlaces distales, preferir el enfoque de 3 pasos
# con corrección del error de clasificación (paquete: BCH, o Mplus)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Enfoque de 3 pasos para desenlace distal
# ──────────────────────────────────────────────────────────────────────────────

# Método de 3 pasos manual (BCH) o vía MplusAutomation
# Paso 1: estimar LCA solo con indicadores (ya hecho: mod3)
# Paso 2: obtener pesos BCH a partir de las posteriores
# Paso 3: regresión del desenlace ponderada por BCH
# En R, el paquete 'BCH' o cálculo manual; alternativa robusta: Mplus
library(MplusAutomation)
# Generar input de Mplus con AUXILIARY = desenlace (DCAT/DCON) (3-step)
# Mplus implementa BCH y DU3STEP automáticamente
# Aproximación simple (solo si entropía MUY alta, > 0.9):
# usar clase modal como factor en el modelo del desenlace
data_cohorte$clase <- factor(mod3$predclass)
library(survival)
fit <- coxph(Surv(tiempo, as.integer(evento == "si")) ~ clase + edad + sexo,
 data = data_cohorte)
summary(fit)
# ADVERTENCIA: este atajo ignora el error de clasificación;
# válido solo con entropía muy alta. Preferir 3 pasos.


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 31 · Random Forest
# ║  Método: rpart/ranger/randomForest/randomForestSRC
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Árbol único (CART)
# ──────────────────────────────────────────────────────────────────────────────

library(rpart)
library(rpart.plot)
library(dplyr)
# Árbol de clasificación
arbol <- rpart(
 evento ~ edad + pa_sistolica + tabaco + colesterol + imc + glucosa,
 data = data_cohorte,
 method = "class",
 control = rpart.control(cp = 0.01, minsplit = 20)
)
# Visualizar
rpart.plot(arbol, type = 4, extra = 104,
 box.palette = "RdYlGn",
 main = "Árbol de clasificación (CART) — probabilidad de evento por hoja")
# Poda por complejidad de coste (cp óptimo por CV)
printcp(arbol)
cp_opt <- arbol$cptable[which.min(arbol$cptable[, "xerror"]), "CP"]
arbol_podado <- prune(arbol, cp = cp_opt)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Random Forest con ranger (rápido)
# ──────────────────────────────────────────────────────────────────────────────

library(ranger) # implementación rápida de RF
library(randomForest) # implementación clásica
# Asegurar que el desenlace es factor (clasificación)
data_cohorte$evento <- as.factor(data_cohorte$evento)
# Ajustar Random Forest
set.seed(123)
rf <- ranger(
 evento ~ edad + pa_sistolica + tabaco + colesterol +
 imc + glucosa + sexo + act_fisica,
 data = data_cohorte,
 num.trees = 1000,
 mtry = 3, # √p aproximado
 importance = "permutation", # MDA (más fiable que impurity)
 probability = TRUE, # predicciones probabilísticas
 respect.unordered.factors = "order",
 seed = 123
)
# Error OOB
rf$prediction.error
# AUC, accuracy, etc. a partir de las predicciones OOB

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ajuste de hiperparámetros (mtry)
# ──────────────────────────────────────────────────────────────────────────────

library(caret)
# Ajuste de mtry por validación cruzada
ctrl <- trainControl(method = "cv", number = 10,
 classProbs = TRUE,
 summaryFunction = twoClassSummary)
grid <- expand.grid(
 mtry = c(2, 3, 4, 5, 6),
 splitrule = "gini",
 min.node.size = c(1, 5, 10)
)
rf_tuned <- train(
 evento ~.,
 data = data_cohorte,
 method = "ranger",
 trControl = ctrl,
 tuneGrid = grid,
 metric = "ROC", # optimizar AUC, no accuracy
 num.trees = 1000
)
print(rf_tuned$bestTune)
# Alternativa rápida: tuneRF de randomForest (solo mtry, por OOB)
randomForest::tuneRF(x = data_cohorte[, predictores], y = data_cohorte$evento,
 ntreeTry = 500, stepFactor = 1.5)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Clases desbalanceadas
# ──────────────────────────────────────────────────────────────────────────────

# Desenlace raro: muestreo balanceado
n_eventos <- sum(data_cohorte$evento == "si")
rf_balanced <- ranger(
 evento ~.,
 data = data_cohorte,
 num.trees = 1000,
 mtry = 3,
 importance = "permutation",
 probability = TRUE,
 case.weights = ifelse(data_cohorte$evento == "si",
 nrow(data_cohorte) / n_eventos, 1), # ponderar
 seed = 123
)
# Evaluar con métricas apropiadas (NO solo accuracy)
library(pROC)
pred_oob <- rf_balanced$predictions[, "si"]
roc_obj <- roc(data_cohorte$evento, pred_oob)
auc(roc_obj)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Importancia de variables
# ──────────────────────────────────────────────────────────────────────────────

# Importancia por permutación
importancia <- ranger::importance(rf)
imp_df <- data.frame(
 variable = names(importancia),
 importancia = importancia
) %>% arrange(desc(importancia))
# Visualizar
library(ggplot2)
ggplot(imp_df, aes(x = reorder(variable, importancia), y = importancia)) +
 geom_col(fill = col_libro["navy"]) +
 coord_flip() +
 labs(title    = "Importancia de variables del Random Forest",
 subtitle = "Importancia por permutación (disminución media de la precisión)",
 x        = NULL,
 y        = "Importancia (permutación)",
 caption  = "Valores mayores indican variables más relevantes para la predicción") +
 theme_libro()
# Importancia condicional (corrige correlación entre predictores)
library(party)
# La fórmula con "." arrastraría comorbilidades (codificadas 1/2),
# 'educacion' y 'tiempo' (variable de seguimiento, no predictora);
# se usa el mismo conjunto de covariables que el resto de la ficha
cf <- cforest(evento ~ edad + pa_sistolica + tabaco + colesterol +
    imc + glucosa + sexo + act_fisica, data = data_cohorte,
 controls = cforest_unbiased(ntree = 500, mtry = 3))
varimp_cond <- varimp(cf, conditional = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Partial dependence y SHAP
# ──────────────────────────────────────────────────────────────────────────────

library(pdp)
library(iml)
# PDP univariante
pdp::partial(rf, pred.var = "edad", plot = TRUE,
 plot.engine = "ggplot2")
# PDP bivariante (interacción)
pd2 <- pdp::partial(rf, pred.var = c("edad", "pa_sistolica"))
plotPartial(pd2)
# Curvas ICE (más robustas que PDP)
pdp::partial(rf, pred.var = "edad", ice = TRUE, plot = TRUE,
 alpha = 0.1)
# SHAP values (interpretación local y global)
predict_ranger_prob <- function(model, newdata) predict(model, data = newdata)$predictions
predictor <- iml::Predictor$new(rf, data = data_cohorte[, predictores],
 y = data_cohorte$evento, predict.function = predict_ranger_prob)
shapley <- iml::Shapley$new(predictor,
 x.interest = data_cohorte[1, predictores])
plot(shapley) +
 labs(title = "Valores de Shapley — contribución de cada variable a la predicción individual") +
 theme_libro()

# ──────────────────────────────────────────────────────────────────────────────
# ■ Random Survival Forest
# ──────────────────────────────────────────────────────────────────────────────

library(randomForestSRC)
# Para datos de supervivencia (tiempo-a-evento con censura)
rsf <- rfsrc(
 Surv(tiempo, evento) ~ edad + pa_sistolica + tabaco + colesterol,
 data = data_cohorte,
 ntree = 1000,
 importance = TRUE
)
# Error OOB (índice C), importancia, predicción de supervivencia
print(rsf)
plot(rsf, main = "Random Survival Forest — error OOB y estimación del efecto")


# ══════════════════════════════════════════════════════════════════════════════
# ║  FICHA 32 · Regresión penalizada
# ║  Método: glmnet · Ridge/Lasso/Elastic Net, doble selección causal
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# ■ Preparación de datos
# ──────────────────────────────────────────────────────────────────────────────

library(glmnet) # regresión penalizada
library(caret) # ajuste de hiperparámetros
library(dplyr)
# glmnet requiere matriz de predictores (x) y vector respuesta (y)
# Variables categóricas → dummies con model.matrix
x <- model.matrix(~ . - 1, data = data_cohorte[, predictores])
y <- data_cohorte$evento
# Para desenlace binario: family = "binomial"
# Para continuo: family = "gaussian"
# Para supervivencia: family = "cox" con y = Surv(tiempo, evento)
# glmnet estandariza automáticamente (standardize = TRUE por defecto)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Lasso con validación cruzada
# ──────────────────────────────────────────────────────────────────────────────

# LASSO (alpha = 1)
set.seed(123)
cv_lasso <- cv.glmnet(
 x, y,
 family = "binomial",
 alpha = 1, # 1 = Lasso, 0 = Ridge
 nfolds = 10,
 type.measure = "deviance" # o "auc", "class", "mse"
)
# Gráfico del error de CV vs lambda
plot(cv_lasso, main = "Validación cruzada: error frente a lambda (Lasso)")
# muestra lambda.min y lambda.1se
# Lambdas óptimos
cv_lasso$lambda.min # mínimo error
cv_lasso$lambda.1se # más parsimonioso
# Coeficientes en lambda.1se (modelo parsimonioso)
coef_1se <- coef(cv_lasso, s = "lambda.1se")
# Las variables con coeficiente != 0 son las seleccionadas
selected <- rownames(coef_1se)[which(coef_1se != 0)]
print(selected)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Ridge
# ──────────────────────────────────────────────────────────────────────────────

# RIDGE (alpha = 0)
set.seed(123)
cv_ridge <- cv.glmnet(
 x, y,
 family = "binomial",
 alpha = 0, # Ridge
 nfolds = 10
)
# Ridge mantiene TODAS las variables (encoge, no anula)
coef_ridge <- coef(cv_ridge, s = "lambda.min")
# Todos los coeficientes != 0, pero encogidos
# Trayectorias de coeficientes
fit_ridge <- glmnet(x, y, family = "binomial", alpha = 0)
# Margen superior ampliado: plot.glmnet() dibuja un eje propio (nº de
# coeficientes activos) sobre el gráfico; sin este margen el título
# personalizado se solapa con esos números.
op <- par(mar = c(5, 4, 6, 2) + 0.1)
plot(fit_ridge, xvar = "lambda", label = TRUE,
 main = "Trayectorias de coeficientes — Ridge")
par(op)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Elastic Net (ajuste de α y λ)
# ──────────────────────────────────────────────────────────────────────────────

# ELASTIC NET: ajustar alpha y lambda conjuntamente con caret
ctrl <- trainControl(method = "cv", number = 10,
 classProbs = TRUE,
 summaryFunction = twoClassSummary)
grid <- expand.grid(
 alpha = seq(0, 1, by = 0.1), # de Ridge (0) a Lasso (1)
 lambda = 10^seq(-4, 0, length = 50)
)
set.seed(123)
enet <- train(
 x = x, y = y,
 method = "glmnet",
 trControl = ctrl,
 tuneGrid = grid,
 metric = "ROC" # optimizar AUC
)
# Mejores hiperparámetros
enet$bestTune # alpha y lambda óptimos
# Coeficientes del modelo final
coef(enet$finalModel, s = enet$bestTune$lambda)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Trayectorias y visualización
# ──────────────────────────────────────────────────────────────────────────────

# Trayectorias de coeficientes del Lasso
fit_lasso <- glmnet(x, y, family = "binomial", alpha = 1)
# vs lambda (escala log)
op <- par(mar = c(5, 4, 6, 2) + 0.1)
plot(fit_lasso, xvar = "lambda", label = TRUE,
 main = "Trayectorias de coeficientes — Lasso (vs. log lambda)")
par(op)
# vs fracción de desviación explicada
op <- par(mar = c(5, 4, 6, 2) + 0.1)
plot(fit_lasso, xvar = "dev", label = TRUE,
 main = "Trayectorias de coeficientes — Lasso (vs. devianza explicada)")
par(op)
# Número de variables activas por lambda
fit_lasso$df

# ──────────────────────────────────────────────────────────────────────────────
# ■ Lasso para Cox (supervivencia)
# ──────────────────────────────────────────────────────────────────────────────

library(survival)
# Penalización en modelos de supervivencia
y_surv <- Surv(data_cohorte$tiempo, as.integer(data_cohorte$evento == "si"))
cv_cox <- cv.glmnet(
 x, y_surv,
 family = "cox",
 alpha = 1,
 nfolds = 10
)
# Variables seleccionadas (predictores de supervivencia)
coef_cox <- coef(cv_cox, s = "lambda.1se")
selected_cox <- rownames(coef_cox)[which(coef_cox != 0)]
print(selected_cox)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Doble selección para inferencia causal
# ──────────────────────────────────────────────────────────────────────────────

library(hdm) # high-dimensional metrics
# Para estimar el efecto de una exposición ajustando por
# muchos confusores potenciales (Belloni et al. 2014)
# Doble selección: Lasso para Y~confusores Y para D~confusores
efecto <- rlassoEffect(
 x = confusores_matrix, # confusores potenciales
 y = desenlace,
 d = exposicion, # exposición de interés
 method = "double selection"
)
summary(efecto)
# Estima el efecto de la exposición con IC VÁLIDOS
# tras la selección de confusores (post-selección correcta)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Paquetes de CRAN (instalación conjunta)
# ──────────────────────────────────────────────────────────────────────────────

# 1) Paquetes de CRAN — copiar y pegar en la consola de R
pkgs_cran <- c(
  "AER", "biotools", "broom", "broom.mixed",
  "ca", "car", "caret", "cluster",
  "casebase", "cmprsk", "cobalt", "DHARMa",
  "dosresmeta", "dplyr", "DTRreg",
  "episensr", "epitools", "EValue", "factoextra",
  "FactoMineR", "flexsurv", "fpc", "gee",
  "geepack", "gfoRmula", "ggplot2", "ggsurvfit",
  "GJRM", "glmmTMB", "glmnet", "GPArotation",
  "hdm", "ieugwasr", "iml", "ipw",
  "ivDiag", "ivmodel", "ivtools", "klaR",
  "lavaan", "lme4", "lmerTest", "lmtest",
  "MatchIt", "mclust", "mecor", "mediation",
  "MendelianRandomization", "meta", "metafor", "mice",
  "miceadds", "modelsummary", "MplusAutomation", "mstate",
  "mvmeta", "MVN", "naniar",
  "NbClust", "nFactors", "parameters", "party", "patchwork",
  "pbkrtest", "pdp", "pec", "performance",
  "pheatmap", "poLCA", "pROC", "prodlim",
  "pscl", "psych", "randomForest", "randomForestSRC",
  "ranger", "regmedint", "remotes", "rms",
  "rpart", "rpart.plot", "rstpm2", "sandwich",
  "sensemakr", "simex", "survey", "survminer",
  "survRM2", "tableone", "tidycmprsk", "tidyr",
  "timeROC", "vcdExtra", "VIM", "WeightIt"
)
# Instala solo los que falten:
nuevos <- pkgs_cran[!(pkgs_cran %in% rownames(installed.packages()))]
if (length(nuevos)) install.packages(nuevos, dependencies = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# ■ Paquetes que no están en CRAN (desde GitHub)
# ──────────────────────────────────────────────────────────────────────────────

# 2) Paquetes desde GitHub (requieren el paquete remotes, ya incluido arriba)
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("rondolab/MR-PRESSO")
remotes::install_github("MathiasHarrer/dmetar")  # dmetar (Ficha 19)
# MVMR (Ficha 12): desde el R-universe de MRCIEU (binarios para R 4.4–4.6)
install.packages("MVMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
# topmodels (rootogram, Fichas 02 y 17): desde el R-universe de Zeileis
install.packages("topmodels", repos = c("https://zeileis.r-universe.dev", "https://cloud.r-project.org"))
# CMAverse (Ficha 18): análisis de mediación causal; no está en CRAN
remotes::install_github("BS1125/CMAverse")

# ──────────────────────────────────────────────────────────────────────────────
# ■ Comprobación rápida
# ──────────────────────────────────────────────────────────────────────────────

# 3) Verificar que todo carga sin errores
invisible(lapply(c(pkgs_cran, "TwoSampleMR", "MRPRESSO", "dmetar", "MVMR", "topmodels"),
  function(p) suppressMessages(require(p, character.only = TRUE))))
R.version.string  # confirma la versión de R instalada
