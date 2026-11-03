## ----setup, include=FALSE-----------------------------------------------------
# vignette dir under R CMD check, package root under an interactive knit
.fo <- c("_fig_optim.R", file.path("vignettes", "_fig_optim.R"))
if (any(file.exists(.fo))) source(.fo[file.exists(.fo)][1])
knitr::opts_chunk$set(
  message    = FALSE,
  warning    = FALSE,
  fig.width  = 7,
  fig.height = 4.5,
  cache      = TRUE,
  cache.path = "explainability_cache/"
)
old_opts <- options(mc.cores = 1, rf.cores = 1)


## ----libs---------------------------------------------------------------------
library(ggplot2)

# Try the installed package first, fall back to pkgload::load_all() for the
# R CMD check vignette rebuild where the package isn't yet on .libPaths().
if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE,
                    attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}


## ----precomputed, include=FALSE, cache=FALSE----------------------------------
# gg_shap() on this forest costs about 320 seconds and gg_partial_varpro()
# about 21; everything else here runs in under three. Both are precomputed by
# precompute_explainability.R so the R CMD check vignette rebuild stays fast.
# If the file is absent, each chunk falls back to live computation.
.pc <- tryCatch({
  .cand <- c("explainability_precomputed.rds",
             file.path("vignettes", "explainability_precomputed.rds"))
  .hit <- .cand[file.exists(.cand)]
  if (length(.hit)) readRDS(.hit[[1]]) else list()
}, error = function(e) list())


## ----forest-------------------------------------------------------------------
data(Boston, package = "MASS")
set.seed(20260909L)
rfsrc_boston <- randomForestSRC::rfsrc(medv ~ ., data = Boston, ntree = 500)
rfsrc_boston


## ----correlations-------------------------------------------------------------
pred <- setdiff(names(Boston), "medv")
cm <- cor(Boston[, pred])
diag(cm) <- NA
head(sort(abs(cm["tax", ]), decreasing = TRUE), 3)


## ----gg-variable--------------------------------------------------------------
gg_dta <- gg_variable(rfsrc_boston)
plot(gg_dta, xvar = "lstat", alpha = 0.4) +
  labs(y = "Median value ($1000s)", x = "% lower status population")


## ----partial-dependence-------------------------------------------------------
pd <- gg_partial_rfsrc(rfsrc_boston,
                       xvar.names = c("crim", "lstat", "ptratio"),
                       n_eval = 25)
plot(pd)


## ----ale----------------------------------------------------------------------
ale <- gg_ale_rfsrc(rfsrc_boston,
                    xvar.names = c("crim", "lstat", "ptratio"),
                    n_eval = 25)
plot(ale)


## ----compare------------------------------------------------------------------
compare <- function(v) {
  a <- ale$continuous[ale$continuous$name == v, ]
  p <- pd$continuous[pd$continuous$name == v, ]
  lo <- max(min(a$x), min(p$x))
  hi <- min(max(a$x), max(p$x))
  grid <- seq(lo, hi, length.out = 100)
  ay <- approx(a$x, a$yhat, xout = grid)$y
  py <- approx(p$x, p$yhat, xout = grid)$y
  ay <- ay - mean(ay)
  py <- py - mean(py)
  data.frame(
    variable    = v,
    ale_range   = round(diff(range(ay)), 2),
    pd_range    = round(diff(range(py)), 2),
    pd_over_ale = round(diff(range(py)) / diff(range(ay)), 2)
  )
}
do.call(rbind, lapply(c("crim", "lstat", "ptratio"), compare))


## ----ale-interaction, fig.height=5--------------------------------------------
ale_int <- gg_ale_rfsrc(rfsrc_boston, xvar.names = "crim",
                        xvar2.name = "lstat", n_eval = 12)
plot(ale_int)


## ----shap---------------------------------------------------------------------
if (!is.null(.pc$shap_boston)) {
  shap_boston <- .pc$shap_boston
} else if (requireNamespace("kernelshap", quietly = TRUE)) {
  set.seed(20260909L)
  rows <- sort(sample(nrow(Boston), 60))
  shap_boston <- gg_shap(rfsrc_boston,
                         newdata = rfsrc_boston$xvar[rows, , drop = FALSE],
                         bg_n = 30)
} else {
  shap_boston <- NULL
}
if (!is.null(shap_boston)) plot(shap_boston)


## ----shap-note, echo=FALSE, results="asis"------------------------------------
# Without the precomputed bundle AND without kernelshap there is no figure
# above, and the prose below would describe a plot that is not on the page.
# Say so rather than leaving a silent gap.
if (is.null(shap_boston)) {
  cat("::: {.callout-note}\n",
      "The SHAP figure is not shown: neither the precomputed bundle ",
      "(`explainability_precomputed.rds`) nor the **kernelshap** package is ",
      "available in this build. Install kernelshap, or run ",
      "`vignettes/precompute_explainability.R`, to render it. The paragraphs ",
      "below describe what it shows.\n:::\n", sep = "")
}


## ----partial-varpro-----------------------------------------------------------
if (!is.null(.pc$pd_varpro)) {
  pd_varpro <- .pc$pd_varpro
} else {
  set.seed(20260909L)
  v_boston <- varPro::varpro(medv ~ ., data = Boston, split.weight = FALSE)
  pd_varpro <- gg_partial_varpro(object = v_boston)
}
plot(pd_varpro)


## -----------------------------------------------------------------------------
#| label: restore-options
#| include: false
#| cache: false
options(old_opts)

