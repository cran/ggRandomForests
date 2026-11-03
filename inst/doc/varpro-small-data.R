## ----setup, include=FALSE-----------------------------------------------------
.fo <- c("_fig_optim.R", file.path("vignettes", "_fig_optim.R"))
if (any(file.exists(.fo))) source(.fo[file.exists(.fo)][1])
knitr::opts_chunk$set(message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 4.5)
old_opts <- options(mc.cores = 1, rf.cores = 1)


## ----libs---------------------------------------------------------------------
library(ggplot2)
if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE,
                    attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}


## ----precomputed, include=FALSE-----------------------------------------------
# Every fit in this vignette is precomputed by precompute_varpro_small.R and
# loaded here, so the R CMD check rebuild grows no forest. If the file is
# missing, the script is sourced (without saving) and the fits run live.
.cand <- c("varpro_small_precomputed.rds",
           file.path("vignettes", "varpro_small_precomputed.rds"))
.hit <- .cand[file.exists(.cand)]
if (length(.hit)) {
  vs <- readRDS(.hit[[1]])
} else {
  .scr <- c("precompute_varpro_small.R", file.path("vignettes", "precompute_varpro_small.R"))
  .varpro_small_no_save <- TRUE
  .env <- new.env()
  sys.source(.scr[file.exists(.scr)][[1]], envir = .env)
  vs <- .env$varpro_small_objects
}
cg <- vs$cov_guided_v
cg_real <- cg[!cg$noise, ]
co <- vs$cov_open_v
co_real <- co[!co$noise, ]
pct <- function(x) sprintf("%.1f%%", 100 * x)
z2 <- function(x) sprintf("%.2f", x)


## ----data-sketch, eval=FALSE--------------------------------------------------
# library(survival)
# vp_guided <- varPro::varpro(Surv(days, status) ~ ., data = d, ntree = 1000)
# vp_open   <- varPro::varpro(Surv(days, status) ~ ., data = d, ntree = 1000,
#                             split.weight = FALSE, mtry = 4)


## ----coverage-guided----------------------------------------------------------
cg_real[, c("variable", "eligible", "n_trees", "coverage", "z", "tested")]


## ----importance-plot, fig.height = 5------------------------------------------
plot(vs$gg_guided_all)
# Live: plot(gg_varpro(vp_guided, max.tree = vp_guided$rf$ntree))


## ----trim-table---------------------------------------------------------------
trim_tab <- data.frame(
  trim = names(vs$trim_z),
  terms_above_0 = vapply(vs$trim_z, function(t) sum(t$z > 0), integer(1)),
  terms_above_0.79 = vapply(vs$trim_z, function(t) sum(t$z > 0.79), integer(1)),
  rank_agreement = vapply(vs$trim_z, function(t) {
    b <- vs$trim_z[[1]]
    round(stats::cor(t$z, b$z[match(t$term, b$term)], method = "spearman"), 2)
  }, numeric(1)))
trim_tab


## ----cutoff-table-------------------------------------------------------------
vz <- cg_real[!is.na(cg_real$z), ]
cuts <- c(0, 0.5, 0.7, 0.79, 0.9, 1, 1.5)
data.frame(cutoff = cuts,
           variables_kept = vapply(cuts, function(k) sum(vz$z > k), integer(1)))
band <- vz[vz$z > 0.5 & vz$z <= 1, ]
band[order(-band$z), c("variable", "z", "coverage")]


## ----noise-open---------------------------------------------------------------
co_noise <- co[co$noise, ]
data.frame(kind = c("continuous", "categorical"),
           noise_tested = c(sum(co_noise$tested & co_noise$kind == "continuous"),
                            sum(co_noise$tested & co_noise$kind == "categorical")),
           largest_noise_z = c(max(co_noise$z[co_noise$kind == "continuous"], na.rm = TRUE),
                               max(co_noise$z[co_noise$kind == "categorical"], na.rm = TRUE)))


## ----rare-precomputed, include=FALSE------------------------------------------
# Summaries from precompute_varpro_rare.R. If the file is missing, the script
# runs one small replicate at one setting, so the check stays fast.
.cand <- c("varpro_rare_precomputed.rds",
           file.path("vignettes", "varpro_rare_precomputed.rds"))
.hit <- .cand[file.exists(.cand)]
if (length(.hit)) {
  vr <- readRDS(.hit[[1]])
} else {
  .scr <- c("precompute_varpro_rare.R",
            file.path("vignettes", "precompute_varpro_rare.R"))
  .varpro_rare_fallback <- TRUE
  .env <- new.env()
  sys.source(.scr[file.exists(.scr)][[1]], envir = .env)
  vr <- .env$varpro_rare_objects
}
rare_table <- function(x) {
  if (!nrow(x)) return(x)
  num <- c("events", "seconds", "tested", "truth_top6", "cov_b1", "cov_b2",
           "z_b1", "z_x4", "z_n1")
  out <- stats::aggregate(x[num], x["setting"], mean)
  out <- out[match(unique(x$setting), out$setting), ]
  out$events <- round(out$events)
  out$seconds <- round(out$seconds)
  out[c("cov_b1", "cov_b2")] <- lapply(out[c("cov_b1", "cov_b2")],
                                      function(v) sprintf("%.0f%%", 100 * v))
  out[c("tested", "truth_top6")] <- lapply(out[c("tested", "truth_top6")],
                                          round, 1)
  out[c("z_b1", "z_x4", "z_n1")] <- lapply(out[c("z_b1", "z_x4", "z_n1")],
                                          round, 2)
  rownames(out) <- NULL
  out
}
rs <- rare_table(vr$small)
rl <- rare_table(vr$large)
# One value from a summary table, or a dash when the fallback did not run it
rv <- function(tab, setting, col) {
  v <- tab[[col]][tab$setting == setting]
  if (length(v)) as.character(v) else "(not computed)"
}


## ----rare-settings, eval=FALSE------------------------------------------------
# p <- 30
# defaults   <- varpro(Surv(time, status) ~ ., d)
# open_all   <- varpro(Surv(time, status) ~ ., d, split.weight = FALSE,
#                      sparse = FALSE)
# open_sqrtp <- varpro(Surv(time, status) ~ ., d, split.weight = FALSE,
#                      sparse = FALSE, mtry = ceiling(sqrt(p)))
# open_more  <- varpro(Surv(time, status) ~ ., d, split.weight = FALSE,
#                      sparse = FALSE, mtry = ceiling(sqrt(p)),
#                      ntree = 1000, max.tree = 1000, max.rules.tree = 300,
#                      nodesize.external = round(2 * nrow(d) / sum(d$status)))


## ----rare-small---------------------------------------------------------------
rs


## ----rare-large---------------------------------------------------------------
rl


## ----rare-recipe, eval=FALSE--------------------------------------------------
# p  <- ncol(d) - 2
# vp <- varpro(Surv(time, status) ~ ., d, split.weight = FALSE, sparse = FALSE,
#              mtry = ceiling(sqrt(p)))
# gg <- gg_varpro(vp, max.tree = vp$rf$ntree)
# varpro_coverage(vp, by = "variable")


## ----partial-local-marginal, fig.height = 3.8---------------------------------
one_var <- function(p, v) {
  p$continuous <- p$continuous[p$continuous$name == v, , drop = FALSE]
  p$categorical <- p$categorical[0, , drop = FALSE]
  p
}
patchwork::wrap_plots(
  plot(one_var(vs$pd_local, "bili")) + labs(title = "varPro (local)"),
  plot(one_var(vs$pd_marginal, "bili")) + labs(title = "randomForestSRC (marginal)"),
  nrow = 1)


## ----partial-refits, fig.height = 3.6-----------------------------------------
refits <- do.call(rbind, lapply(vs$pd_refits, function(p) {
  cbind(p$continuous, refit = factor(p$refit))
}))
full <- vs$pd_marginal$continuous
full <- full[full$name == "bili", ]
ggplot(refits, aes(x, yhat, group = refit)) +
  geom_line(colour = "grey60") +
  geom_line(data = full, aes(x, yhat), inherit.aes = FALSE, linewidth = 1) +
  labs(x = "bili", y = sprintf("Survival at %d days", vs$tau),
       title = "Marginal partial: full fit (black) and five subsample refits (grey)")


## ----partial-factor-----------------------------------------------------------
ct <- vs$pd_local$categorical
lvl <- c("stage3", "stage4", "edema0", "edema1")
level_median <- function(l, at) {
  median(ct$parametric[ct$name == l & ct$variable == at], na.rm = TRUE)
}
data.frame(level = lvl,
           partial_at_0 = vapply(lvl, level_median, numeric(1), at = 0),
           partial_at_1 = vapply(lvl, level_median, numeric(1), at = 1),
           row.names = NULL)


## -----------------------------------------------------------------------------
#| label: restore-options
#| include: false
#| cache: false
options(old_opts)

