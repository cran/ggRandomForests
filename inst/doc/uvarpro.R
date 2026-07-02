## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  message    = FALSE,
  warning    = FALSE,
  fig.width  = 7,
  fig.height = 4.5,
  cache      = TRUE,
  cache.path = "uvarpro_cache/"
)
options(mc.cores = 1, rf.cores = 1)


## ----libs---------------------------------------------------------------------
library(ggplot2)

# Match the pattern used by the other vignettes: try the installed
# package first, fall back to pkgload::load_all() for the R CMD check
# vignette rebuild where the package isn't yet on .libPaths(). All varPro
# calls below are ::-qualified, so no library(varPro) is needed.
if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE,
                    attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}


## ----fit----------------------------------------------------------------------
set.seed(1)
u <- varPro::uvarpro(mtcars, ntree = 50)


## ----beta_fit-----------------------------------------------------------------
beta_fit <- varPro::get.beta.entropy(u)


## ----udep, eval = requireNamespace("ggraph", quietly = TRUE)------------------
plot(gg_udependent(u))


## ----beta---------------------------------------------------------------------
plot(gg_beta_uvarpro(u, beta_fit = beta_fit))


## ----sdep---------------------------------------------------------------------
plot(gg_sdependent(u, beta_fit = beta_fit))

