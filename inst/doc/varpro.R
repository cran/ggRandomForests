## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  message    = FALSE,
  warning    = FALSE,
  fig.width  = 7,
  fig.height = 4.5,
  cache      = TRUE,
  cache.path = "varpro_cache/"
)
options(mc.cores = 1, rf.cores = 1)


## ----libs---------------------------------------------------------------------
library(ggplot2)

# Match the pattern used by the regression/survival vignettes: try the
# installed package first, fall back to pkgload::load_all() for the
# R CMD check vignette rebuild where the package isn't yet on
# .libPaths(). All varPro calls below are ::-qualified, so no
# library(varPro) is needed.
if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE,
                    attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}


## ----precomputed, include=FALSE, cache=FALSE----------------------------------
# The three gg_partial_varpro() calls and the Boston beta.varpro() fit are
# the only expensive parts of this vignette (~34s combined). They are
# precomputed by precompute_varpro.R and loaded here so the R CMD check
# vignette rebuild stays fast. If the file is absent (e.g. a fresh clone
# that hasn't run the script), each chunk falls back to live computation,
# so the vignette is reproducible either way.
.vp <- tryCatch({
  # R CMD check rebuilds the vignette in a copy whose working directory is the
  # vignette dir; an interactive knit from the package root sees vignettes/
  # instead. Try both locations, and fall back to live computation below if
  # neither is present or readRDS() fails.
  .cand <- c("varpro_precomputed.rds",
             file.path("vignettes", "varpro_precomputed.rds"))
  .hit <- .cand[file.exists(.cand)]
  if (length(.hit)) readRDS(.hit[[1]]) else list()
}, error = function(e) list())   # missing or unreadable -> live fallback below


## ----boston-fit---------------------------------------------------------------
data("Boston", package = "MASS")
set.seed(20260527L)
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
v_boston <- if (is.null(.vp$v_boston)) {
  varPro::varpro(medv ~ ., data = Boston, ntree = 50)
} else {
  .vp$v_boston
}
v_boston


## ----boston-gg-varpro---------------------------------------------------------
plot(gg_varpro(v_boston))


## ----boston-gg-partial-varpro-------------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
gg_pd <- if (is.null(.vp$pd_boston)) {
  gg_partial_varpro(object = v_boston)
} else {
  .vp$pd_boston
}
plot(gg_pd)


## ----boston-beta-fit, cache=TRUE----------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
b_boston <- if (is.null(.vp$b_boston)) {
  varPro::beta.varpro(v_boston)
} else {
  .vp$b_boston
}


## ----boston-gg-beta-varpro----------------------------------------------------
plot(gg_beta_varpro(v_boston, beta_fit = b_boston))


## ----boston-uvarpro, cache=TRUE-----------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
u_boston <- if (is.null(.vp$u_boston)) {
  varPro::uvarpro(Boston[, setdiff(names(Boston), "medv")], ntree = 50)
} else {
  .vp$u_boston
}


## ----boston-gg-udependent, eval = requireNamespace("ggraph", quietly = TRUE)----
plot(gg_udependent(u_boston))


## ----boston-isopro, cache=TRUE------------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
iso_boston <- if (is.null(.vp$iso_boston)) {
  varPro::isopro(data = Boston[, setdiff(names(Boston), "medv")],
                 method = "rnd", sampsize = 256, ntree = 50)
} else {
  .vp$iso_boston
}


## ----boston-gg-isopro---------------------------------------------------------
plot(gg_isopro(iso_boston))


## ----boston-ivarpro-fit, cache=TRUE-------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
iv_boston <- if (is.null(.vp$iv_boston)) {
  varPro::ivarpro(v_boston)
} else {
  .vp$iv_boston
}


## ----boston-gg-ivarpro-distribution-------------------------------------------
plot(gg_ivarpro(v_boston, ivarpro_fit = iv_boston))


## ----boston-gg-ivarpro-which-obs----------------------------------------------
plot(gg_ivarpro(v_boston, ivarpro_fit = iv_boston, which_obs = 1L))


## ----iris-fit-binary----------------------------------------------------------
iris_binary <- iris[iris$Species != "setosa", ]
iris_binary$Species <- droplevels(iris_binary$Species)
set.seed(20260527L)
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
v_iris_binary <- if (is.null(.vp$v_iris_binary)) {
  varPro::varpro(Species ~ ., data = iris_binary, ntree = 50)
} else {
  .vp$v_iris_binary
}


## ----iris-fit-multi-----------------------------------------------------------
set.seed(20260527L)
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
v_iris_multi <- if (is.null(.vp$v_iris_multi)) {
  varPro::varpro(Species ~ ., data = iris, ntree = 50)
} else {
  .vp$v_iris_multi
}


## ----iris-multi-gg-varpro-conditional-----------------------------------------
plot(gg_varpro(v_iris_multi, conditional = TRUE))


## ----iris-multi-gg-partial-varpro---------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
gg_pd_iris <- if (is.null(.vp$pd_iris_multi)) {
  gg_partial_varpro(object = v_iris_multi)
} else {
  .vp$pd_iris_multi
}
plot(gg_pd_iris)


## ----iris-binary-beta-fit, cache=TRUE-----------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
b_iris_binary <- if (is.null(.vp$b_iris_binary)) {
  varPro::beta.varpro(v_iris_binary)
} else {
  .vp$b_iris_binary
}


## ----iris-binary-gg-beta-varpro-----------------------------------------------
plot(gg_beta_varpro(v_iris_binary, beta_fit = b_iris_binary))


## ----iris-multi-beta-fit, cache=TRUE------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
b_iris_multi <- if (is.null(.vp$b_iris_multi)) {
  varPro::beta.varpro(v_iris_multi)
} else {
  .vp$b_iris_multi
}


## ----iris-multi-gg-beta-varpro------------------------------------------------
plot(gg_beta_varpro(v_iris_multi, beta_fit = b_iris_multi))


## ----pbc-data-----------------------------------------------------------------
library(survival)
data(pbc, package = "randomForestSRC")
pbc_small <- pbc[, c("days", "status", "age", "albumin", "bili",
                     "edema", "platelet")]
pbc_small <- na.omit(pbc_small)
set.seed(20260527L)
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
v_pbc <- if (is.null(.vp$v_pbc)) {
  varPro::varpro(Surv(days, status) ~ ., data = pbc_small, ntree = 50)
} else {
  .vp$v_pbc
}


## ----pbc-gg-varpro------------------------------------------------------------
plot(gg_varpro(v_pbc))


## ----pbc-gg-partial-varpro----------------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
gg_pd_pbc <- if (is.null(.vp$pd_pbc)) {
  gg_partial_varpro(object = v_pbc)
} else {
  .vp$pd_pbc
}
plot(gg_pd_pbc)


## ----pbc-isopro, cache=TRUE---------------------------------------------------
# Precomputed offline (see precompute_varpro.R); falls back to a live fit.
iso_pbc <- if (is.null(.vp$iso_pbc)) {
  varPro::isopro(data = pbc_small[, c("age", "albumin", "bili", "platelet")],
                 method = "rnd", sampsize = 256, ntree = 50)
} else {
  .vp$iso_pbc
}
plot(gg_isopro(iso_pbc))

