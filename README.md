# Manual Práctico de Bioestadística y Epidemiología — Código en R

Scripts en R que acompañan al *Manual Práctico de Bioestadística y Epidemiología*.
Cada script corresponde a un bloque temático del libro y contiene la sintaxis
anotada de las fichas de ese bloque, con datos de ejemplo (reales o simulados con
semilla fija) para que cada fragmento sea reproducible de forma independiente.

## Requisitos

- **R ≥ 4.4.0** (el material se desarrolló y validó con R 4.4.x).
- Los paquetes necesarios se instalan de una sola vez con el bloque de
  instalación del **Apéndice B** del libro. Consúltalo antes de ejecutar los
  scripts.
- Para fijar versiones exactas de forma persistente puede usarse
  [`renv`](https://rstudio.github.io/renv/) (`renv::init()` + `renv::snapshot()`).

## Estructura

El libro se organiza en 8 bloques temáticos (32 fichas):

| Bloque | Script | Contenido |
|:------:|--------|-----------|
| 1 | `Bloque 1 - Modelos de regresion clasicos.R` | Regresión logística, Poisson/binomial negativa, lineal |
| 2 | `Bloque 2 - Analisis de supervivencia.R` | Kaplan-Meier y log-rank, Cox, modelos paramétricos, riesgos competitivos |
| 3 | `Bloque 3 - Disenos epidemiologicos.R` | Casos-controles, caso-cohorte, cohorte |
| 4 | `Bloque 4 - Inferencia causal.R` | Variables instrumentales, aleatorización mendeliana, propensity score, G-methods |
| 5 | `Bloque 5 - Datos especiales.R` | Modelos mixtos, GEE, ZIP/ZINB, mediación |
| 6 | `Bloque 6 - Meta-analisis.R` | Efectos fijos y aleatorios, dosis-respuesta, meta-regresión |
| 7 | `Bloque 7 - Sesgos y analisis de sensibilidad.R` | E-value, errores de clasificación, imputación múltiple (MICE) |
| 8 | `Bloque 8 - Multivariantes y clasificacion.R` | PCA, EFA, conglomerados, LDA/QDA, MCA, LCA, random forest, regresión penalizada |

## Uso

Cada script está pensado como una hoja de consulta rápida: los fragmentos son
autocontenidos y pueden ejecutarse de forma aislada. Cada paquete se carga junto
al método que lo utiliza, de modo que al copiar un fragmento concreto se ve su
dependencia en ese mismo punto.

Antes de generar las gráficas, cargar una vez por sesión el módulo de identidad
visual (paleta `col_libro` y tema `theme_libro()`) que figura en el Apéndice B.

## Reproducibilidad

Los bloques que simulan datos fijan la semilla con `set.seed()` al inicio de cada
simulación. Los bloques que emplean conjuntos de datos reales (p. ej. `pbc`,
`lung`, `Insurance`, `Boston`) no requieren semilla. Guarda la salida de
`sessionInfo()` junto a tus resultados para documentar el entorno de ejecución.

## Guías de reporte

El material remite a las guías internacionales aplicables: **STROBE** (estudios
observacionales, Bloque 3), **PRISMA 2020** (revisiones sistemáticas y
meta-análisis, Bloque 6) y **STROBE-MR** (aleatorización mendeliana, Ficha 12).

## Cómo citar

[![DOI](https://zenodo.org/badge/1290744800.svg)](https://doi.org/10.5281/zenodo.21216092)

Cite este repositorio con el DOI de concepto anterior (`10.5281/zenodo.21216092`),
que resuelve siempre a la última versión publicada. Si se referencia el código
exacto empleado en un análisis, cite la versión específica correspondiente
(`v1.0.0`: `10.5281/zenodo.21216093`).

Los metadatos de citación (autores, ORCID) están en [`CITATION.cff`](CITATION.cff);
GitHub genera a partir de él la cita en formato APA o BibTeX (botón «Cite this
repository» en la barra lateral del repositorio).

## Autoría

- José María Huerta Castaño — ORCID [0000-0002-9637-3869](https://orcid.org/0000-0002-9637-3869)
- Sandra Milena Colorado Yohar — ORCID [0000-0002-6700-0780](https://orcid.org/0000-0002-6700-0780)

## Licencia

Ver el archivo [`LICENSE`](LICENSE).

---

> **Nota sobre los nombres de archivo.** Los scripts conservan nombres con
> espacios para coincidir con la nomenclatura del libro. Si prefieres rutas más
> robustas para línea de comandos y enlaces, considera renombrarlos sin espacios
> (p. ej. `bloque-1-modelos-regresion-clasicos.R`) y actualizar la tabla anterior.
