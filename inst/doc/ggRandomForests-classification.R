## ----setup, include=FALSE-----------------------------------------------------
# vignette dir under R CMD check, package root under an interactive knit
.fo <- c("_fig_optim.R", file.path("vignettes", "_fig_optim.R"))
if (any(file.exists(.fo))) source(.fo[file.exists(.fo)][1])
knitr::opts_chunk$set(
  fig.align = "center",
  fig.width = 7,
  fig.height = 5,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)
options(mc.cores = 1, rf.cores = 1)


## ----packages-----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(randomForestSRC)

if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}

theme_set(theme_bw())


## ----data-load----------------------------------------------------------------
data(iris)


## ----eda-petal, fig.height=4, fig.cap="Petal length vs. petal width, colored by species."----
ggplot(iris, aes(x = Petal.Length, y = Petal.Width, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set2")


## ----eda-sepal, fig.height=4, fig.cap="Sepal length vs. sepal width, colored by species."----
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set2")


## ----rfsrc-fit----------------------------------------------------------------
set.seed(42)
rfsrc_iris <- rfsrc(Species ~ ., data = iris,
                    ntree = 100, importance = TRUE, err.block = 5)
rfsrc_iris


## ----error-plot, fig.height=4, fig.cap="OOB misclassification rate vs. number of trees, overall and per class."----
plot(gg_error(rfsrc_iris))


## ----rfsrc-predicted, fig.height=5, fig.cap="OOB predicted class probabilities, faceted by the flower's true species."----
plot(gg_rfsrc(rfsrc_iris))


## ----vimp-plot, fig.cap="VIMP ranking, overall and per class."----------------
plot(gg_vimp(rfsrc_iris))


## ----varsel-------------------------------------------------------------------
md_iris <- max.subtree(rfsrc_iris)


## ----shap-fit-----------------------------------------------------------------
set.seed(43)
gg_shp <- gg_shap(rfsrc_iris, bg_n = 30, which.class = 3)


## ----shap-importance-plot, fig.cap="Mean absolute SHAP value per predictor, virginica probability."----
plot(gg_shp, type = "importance")


## ----shap-beeswarm-plot, fig.cap="SHAP beeswarm for virginica probability: every dot is one flower's contribution for one predictor."----
plot(gg_shp, type = "beeswarm")


## ----shap-dependence-plot, fig.cap="SHAP dependence for Petal.Width, virginica probability."----
plot(gg_shp, type = "dependence", xvar = "Petal.Width")


## ----vardep-panel, fig.cap="Variable dependence for Petal.Width, faceted by predicted-class probability."----
gg_v <- gg_variable(rfsrc_iris)

plot(gg_v, xvar = "Petal.Width")


## ----partial-dep, fig.cap="Partial dependence for petal measurements, setosa probability."----
pd <- gg_partial_rfsrc(rfsrc_iris, xvar.names = c("Petal.Length", "Petal.Width"))
plot(pd)


## ----roc-plot, fig.height=4, fig.cap="ROC curve for virginica vs. the rest."----
roc_virginica <- gg_roc(rfsrc_iris, which_outcome = 3)
plot(roc_virginica)


## ----auc-calc-----------------------------------------------------------------
calc_auc(roc_virginica)

