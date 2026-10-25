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
# MASS is in Suggests; the Boston chunks below run only when it is installed
has_mass <- requireNamespace("MASS", quietly = TRUE)


## ----error-demo---------------------------------------------------------------
library(randomForest)
set.seed(42)
rf_iris <- randomForest(Species ~ ., data = iris, ntree = 100, keep.forest = TRUE)
err_df <- ggRandomForests::gg_error(rf_iris, training = TRUE)
head(err_df)


## ----error-plot, fig.height=4-------------------------------------------------
plot(err_df)

