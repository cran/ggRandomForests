## ----setup, include=FALSE-----------------------------------------------------
# vignette dir under R CMD check, package root under an interactive knit
.fo <- c("_fig_optim.R", file.path("vignettes", "_fig_optim.R"))
if (any(file.exists(.fo))) source(.fo[file.exists(.fo)][1])
knitr::opts_chunk$set(
  fig.align = "center",
  fig.width = 7,
  fig.height = 4.5,
  message = FALSE,
  warning = FALSE,
  comment = "#>"
)
old_opts <- options(mc.cores = 1, rf.cores = 1)


## ----packages, include=FALSE--------------------------------------------------
.has_rhf_extractor <- requireNamespace(
  "ggRandomForests",
  quietly = TRUE
) && "gg_rhf" %in% getNamespaceExports("ggRandomForests")
if (.has_rhf_extractor) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(
    export_all = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE
  )
} else {
  stop(
    "Install ggRandomForests (or pkgload for dev builds) to render this vignette."
  )
}


## ----load-bundle, include=FALSE, cache=FALSE----------------------------------
.required_bundle_names <- c(
  "data", "fit", "auct_cumulative", "auct_incident", "importance",
  "tune_risk", "tune_iauc", "seed", "settings", "versions"
)
.required_setting_names <- c(
  "formula", "fit", "auct_cumulative", "auct_incident",
  "importance_cache", "importance_time_index", "tune_risk", "tune_iauc"
)
.valid_rhf_bundle <- function(x) {
  if (!is.list(x)) return(FALSE)
  if (!identical(names(x), .required_bundle_names)) return(FALSE)
  checks <- c(
    is.data.frame(x$data),
    inherits(x$fit, "rhf"),
    inherits(x$auct_cumulative, "auct.rhf"),
    inherits(x$auct_incident, "auct.rhf"),
    inherits(x$importance, "importance.rhf"),
    inherits(x$tune_risk, "tune.treesize.rhf"),
    inherits(x$tune_iauc, "tune.treesize.rhf"),
    is.integer(x$seed),
    length(x$seed) == 1L,
    is.list(x$settings),
    identical(names(x$settings), .required_setting_names),
    is.character(x$versions),
    length(x$versions) > 0L
  )
  all(checks)
}

.bundle_paths <- c(
  "rhf_precomputed.rds",
  file.path("vignettes", "rhf_precomputed.rds")
)
bundle <- NULL
for (.bundle_path in .bundle_paths[file.exists(.bundle_paths)]) {
  .candidate_bundle <- tryCatch(
    readRDS(.bundle_path),
    error = function(e) NULL
  )
  if (.valid_rhf_bundle(.candidate_bundle)) {
    bundle <- .candidate_bundle
    break
  }
}
if (is.null(bundle)) {
  stop(
    paste(
      "The RHF vignette needs a valid precomputed bundle.",
      "Run vignettes/precompute_rhf.R from the package root."
    )
  )
}


## ----subject-intervals--------------------------------------------------------
rhf_data <- bundle$data
one_subject <- rhf_data[rhf_data$id == rhf_data$id[1L],
                        c("id", "start", "stop", "event", "xtd")]
knitr::kable(one_subject, digits = 3)


## ----data-summary-------------------------------------------------------------
records_per_subject <- as.integer(table(rhf_data$id))
events_by_subject <- tapply(rhf_data$event, rhf_data$id, max)
data_summary <- data.frame(
  measure = c(
    "Subjects",
    "Counting-process records",
    "Records per subject",
    "Subjects with an event",
    "xtd range"
  ),
  value = c(
    length(unique(rhf_data$id)),
    nrow(rhf_data),
    sprintf(
      "%d to %d (median %g)",
      min(records_per_subject),
      max(records_per_subject),
      stats::median(records_per_subject)
    ),
    sum(events_by_subject),
    sprintf("%.3f to %.3f", min(rhf_data$xtd), max(rhf_data$xtd))
  )
)
knitr::kable(data_summary, col.names = c("Measure", "Value"))


## ----preparation, eval=FALSE--------------------------------------------------
# set.seed(20260825L)
# sim <- randomForestRHF::hazard.simulation(1)
# data <- sim$dta
# formula <- stats::as.formula("Surv(id, start, stop, event) ~ .")
# with(data, stopifnot(identical(xtd, (x.4 + x.5) * stop)))
# 
# fit <- randomForestRHF::rhf(
#   formula, data, ntree = 50L, seed = -1L
# )
# auct_cumulative <- randomForestRHF::auct.rhf(
#   fit, marker = "cumhaz", method = "cumulative", verbose = FALSE
# )
# auct_incident <- randomForestRHF::auct.rhf(
#   fit, marker = "hazard", method = "incident",
#   riskset = "subject", verbose = FALSE
# )
# cache <- randomForestRHF::varpro.cache.rhf(
#   fit, max.rules.tree = 30L, max.tree = 20L, verbose = FALSE
# )
# time_index <- unique(as.integer(round(seq.int(
#   1L, cache$K, length.out = 5L
# ))))
# importance <- randomForestRHF::importance.rhf(
#   fit, cache = cache, time.index = time_index, verbose = FALSE
# )
# tune_risk <- randomForestRHF::tune.treesize.rhf(
#   formula, data, ntree = 20L, perf = "risk", lower = 2L,
#   upper = 6L, max.evals = 5L, seed = 20260825L, verbose = FALSE,
#   forest = FALSE
# )
# tune_iauc <- randomForestRHF::tune.iAUC.rhf(
#   formula, data, ntree = 20L, lower = 2L, upper = 6L,
#   max.evals = 5L, seed = 20260825L, verbose = FALSE,
#   forest = FALSE
# )


## ----fit-summary--------------------------------------------------------------
fit_summary <- data.frame(
  measure = c(
    "Object class",
    "Trees",
    "Subjects",
    "Working time points",
    "Saved ensemble source"
  ),
  value = c(
    paste(class(bundle$fit), collapse = ", "),
    bundle$fit$ntree,
    length(bundle$fit$ensemble.id),
    length(bundle$fit$time.interest),
    "out-of-bag hazard and cumulative hazard"
  )
)
knitr::kable(fit_summary, col.names = c("Fit component", "Value"))


## ----extract-rhf--------------------------------------------------------------
rhf_curves <- ggRandomForests::gg_rhf(bundle$fit)
head(as.data.frame(rhf_curves))


## ----hazard-plot--------------------------------------------------------------
case_ids <- c(1, 2, 3, 11, 14, 22)
plot(rhf_curves, idx = case_ids)


## ----cumulative-hazard-plot---------------------------------------------------
plot(rhf_curves, idx = case_ids, hazard.only = FALSE)


## ----auct-cumulative----------------------------------------------------------
auct_cumulative <- gg_auct(
  bundle$fit,
  marker = "chf",
  method = "cumulative",
  auct_fit = bundle$auct_cumulative
)
plot(auct_cumulative)


## ----auct-incident------------------------------------------------------------
auct_incident <- gg_auct(
  bundle$fit,
  marker = "haz",
  method = "incident",
  riskset = "subject",
  auct_fit = bundle$auct_incident
)
plot(auct_incident)


## ----rhf-importance, fig.width=8, fig.height=5.25-----------------------------
rhf_priority <- gg_rhf_importance(
  bundle$fit,
  importance_fit = bundle$importance
)
plot(rhf_priority) +
  ggplot2::theme(plot.margin = ggplot2::margin(l = 80))


## ----tune-risk----------------------------------------------------------------
risk_tuning <- gg_tune_rhf(bundle$tune_risk)
risk_tuning
plot(risk_tuning)


## ----tune-iauc----------------------------------------------------------------
iauc_tuning <- gg_tune_rhf(bundle$tune_iauc)
iauc_tuning
plot(iauc_tuning)


## ----extractor-support--------------------------------------------------------
support <- data.frame(
  extractor = c(
    "gg_rhf()",
    "gg_auct()",
    "gg_rhf_importance()",
    "gg_tune_rhf()"
  ),
  upstream_input = c(
    "rhf fit",
    "rhf fit plus an optional auct.rhf object",
    "rhf fit plus an optional importance.rhf object",
    "tune.treesize.rhf object"
  ),
  returned_class = c(
    "gg_rhf, data.frame",
    "gg_auct, data.frame",
    "gg_rhf_importance, data.frame",
    "gg_tune_rhf, data.frame"
  ),
  scale = c(
    paste(
      "Hazard is a local event rate; cumulative hazard is accumulated",
      "event pressure. Neither is a probability."
    ),
    paste(
      "AUC(t) and iAUC are ranking measures from 0 to 1; interpretation",
      "depends on the case definition and marker."
    ),
    paste(
      "Time-local rule-release priority; larger values mean a larger",
      "working-response contrast, not statistical significance."
    ),
    paste(
      "OOB risk is minimized or OOB iAUC is maximized; the selected",
      "tree size is marked."
    )
  ),
  check.names = FALSE
)
knitr::kable(
  support,
  col.names = c(
    "Extractor", "Upstream input class or object",
    "Returned class", "Output scale and interpretation"
  )
)


## ----bundle-versions----------------------------------------------------------
version_table <- data.frame(
  package = names(bundle$versions),
  version = unname(bundle$versions),
  row.names = NULL
)
knitr::kable(
  version_table,
  col.names = c("Software", "Saved version")
)


## -----------------------------------------------------------------------------
#| label: restore-options
#| include: false
#| cache: false
options(old_opts)

