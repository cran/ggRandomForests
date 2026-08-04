## ----pkg-setup, include=FALSE-------------------------------------------------
# vignette dir under R CMD check, package root under an interactive knit
.fo <- c("_fig_optim.R", file.path("vignettes", "_fig_optim.R"))
if (any(file.exists(.fo))) source(.fo[file.exists(.fo)][1])
if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}


## ----error-demo---------------------------------------------------------------
library(randomForest)
set.seed(42)
rf_iris <- randomForest(Species ~ ., data = iris, ntree = 200, keep.forest = TRUE)
err_df <- ggRandomForests::gg_error(rf_iris, training = TRUE)
head(err_df)


## ----error-plot, fig.height=4-------------------------------------------------
plot(err_df)


## ----variable-demo------------------------------------------------------------
set.seed(99)
boston <- MASS::Boston
rf_boston <- randomForest(medv ~ ., data = boston, ntree = 150)
var_df <- ggRandomForests::gg_variable(rf_boston)
str(var_df[, c("lstat", "yhat")])


## ----variable-plot, fig.height=4----------------------------------------------
plot(var_df, xvar = "lstat")


## ----vimp-demo----------------------------------------------------------------
vimp_df <- ggRandomForests::gg_vimp(rf_boston)
head(vimp_df)
plot(vimp_df)


## ----quantile-demo------------------------------------------------------------
rm_breaks <- ggRandomForests::quantile_pts(boston$rm, groups = 6, intervals = TRUE)
rm_groups <- cut(boston$rm, breaks = rm_breaks)
table(rm_groups)

