## ----setup, include=FALSE-----------------------------------------------------
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
library(tidyr)
library(randomForestSRC)
library(survival)

if (requireNamespace("ggRandomForests", quietly = TRUE)) {
  library(ggRandomForests)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
} else {
  stop("Install ggRandomForests (or pkgload for dev builds) to render this vignette.")
}

theme_set(theme_bw())

# Plotting constants
event_marks  <- c(1, 4)
event_labels <- c("Censored", "Death")
event_colors <- c("steelblue", "firebrick")


## ----data-load----------------------------------------------------------------
data("pbc", package = "randomForestSRC")


## ----data-clean---------------------------------------------------------------
pbc <- pbc |>
  mutate(
    years     = days / 365.25,
    age       = age / 365.25,
    treatment = factor(
      ifelse(treatment == 1, "DPCA",
             ifelse(treatment == 2, "Placebo", NA)),
      levels = c("DPCA", "Placebo")
    )
  ) |>
  select(-days)

# Low-cardinality numerics (including binary 0/1) to factor.
# NOTE: do NOT convert to logical — partial.rfsrc() fails with logical
# predictors in survival forests (randomForestSRC <= 3.5.1).
# Exclude the response columns (status, years) from conversion.
resp_cols <- c("status", "years")
for (nm in setdiff(names(pbc), resp_cols)) {
  v <- pbc[[nm]]
  if (is.numeric(v) && !is.factor(v) && length(unique(v[!is.na(v)])) <= 5) {
    pbc[[nm]] <- factor(v)
  }
}


## ----variable-labels----------------------------------------------------------
# Human-readable labels for plot axes
st_labs <- c(
  status       = "Death Event",
  treatment    = "Treatment",
  age          = "Age (years)",
  sex          = "Female",
  ascites      = "Ascites",
  hepato       = "Hepatomegaly",
  spiders      = "Spiders",
  edema        = "Edema (0, 0.5, 1)",
  bili         = "Serum Bilirubin (mg/dl)",
  chol         = "Serum Cholesterol (mg/dl)",
  albumin      = "Albumin (gm/dl)",
  copper       = "Urine Copper (ug/day)",
  alk.phos     = "Alkaline Phosphatase (U/liter)",
  ast          = "SGOT (U/ml)",
  trig         = "Triglycerides (mg/dl)",
  platelet     = "Platelets (per cubic ml/1000)",
  prothrombin  = "Prothrombin Time (sec)",
  stage        = "Histologic Stage",
  years        = "Follow-up Time (years)"
)


## ----eda-categorical, fig.height=4, fig.cap="EDA for categorical variables. Bar color indicates factor level; white = missing."----
cls <- sapply(pbc, class)
cnt_idx <- which(cls %in% c("numeric", "integer"))
fct_idx <- setdiff(seq_along(pbc), cnt_idx)
fct_idx <- union(fct_idx, which(names(pbc) == "years"))

dta_cat <- suppressWarnings(
  pbc[, fct_idx] |>
    pivot_longer(-years, names_to = "variable", values_to = "value",
                 values_transform = list(value = as.character))
)

ggplot(dta_cat, aes(x = years, fill = value)) +
  geom_histogram(color = "black", binwidth = 1) +
  labs(y = "", x = st_labs["years"]) +
  scale_fill_brewer(palette = "RdBu", na.value = "white") +
  facet_wrap(~variable, scales = "free_y", nrow = 2) +
  theme(legend.position = "none")


## ----eda-continuous, fig.height=4, fig.cap="EDA for continuous variables. Points colored by death event; rug marks indicate missing values."----
cnt_idx <- union(cnt_idx, which(names(pbc) == "status"))
dta_cont <- pbc[, cnt_idx] |>
  pivot_longer(c(-years, -status),
               names_to = "variable", values_to = "value")

ggplot(dta_cont |> filter(!is.na(value)),
       aes(x = years, y = value, color = factor(status), shape = factor(status))) +
  geom_point(alpha = 0.4) +
  geom_rug(data = dta_cont |> filter(is.na(value)), color = "grey50") +
  labs(y = "", x = st_labs["years"], color = "Death", shape = "Death") +
  scale_color_manual(values = event_colors) +
  scale_shape_manual(values = event_marks) +
  facet_wrap(~variable, scales = "free_y", ncol = 4) +
  theme(legend.position = c(0.8, 0.2))


## ----km-survival, fig.cap="KM survival estimates by treatment group with 95% confidence bands."----
pbc_trial <- pbc |> filter(!is.na(treatment))
pbc_test  <- pbc |> filter(is.na(treatment))

gg_km <- gg_survival(interval = "years", censor = "status",
                     by = "treatment", data = pbc_trial,
                     conf.int = 0.95)

plot(gg_km) +
  labs(y = "Survival Probability", x = "Time (years)",
       color = "Treatment", fill = "Treatment") +
  theme(legend.position = c(0.2, 0.2)) +
  coord_cartesian(ylim = c(0, 1.01))


## ----km-cumhaz, fig.cap="KM cumulative hazard estimates by treatment group."----
plot(gg_km, type = "cum_haz") +
  labs(y = "Cumulative Hazard", x = "Time (years)",
       color = "Treatment", fill = "Treatment") +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(ylim = c(-0.02, 1.22))


## ----km-bili, fig.width=6, fig.cap="KM survival stratified by bilirubin groups."----
pbc_bili <- pbc_trial |>
  mutate(bili_grp = cut(bili, breaks = c(0, 0.8, 1.3, 3.4, 29)))

plot(gg_survival(interval = "years", censor = "status",
                 by = "bili_grp", data = pbc_bili),
     error = "none") +
  labs(y = "Survival Probability", x = "Time (years)", color = "Bilirubin")


## ----rfsrc-fit----------------------------------------------------------------
# Step 1: impute missing values via random forest proximity
pbc_imputed <- impute.rfsrc(Surv(years, status) ~ .,
                             data    = pbc_trial,
                             ntree   = 500,
                             nimpute = 2)

# Step 2: grow the survival forest on the complete imputed data
rfsrc_pbc <- rfsrc(Surv(years, status) ~ .,
                   data      = pbc_imputed,
                   nsplit    = 10,
                   tree.err  = TRUE,
                   importance = TRUE)


## ----error-plot, fig.height=4, fig.cap="OOB prediction error vs. number of trees."----
plot(gg_error(rfsrc_pbc))


## ----rfsrc-predicted, fig.cap="OOB predicted survival curves. Blue = censored, red = death events."----
gg_rsf <- plot(gg_rfsrc(rfsrc_pbc), alpha = 0.2) +
  scale_color_manual(values = event_colors) +
  theme(legend.position = "none") +
  labs(y = "Survival Probability", x = "Time (years)") +
  coord_cartesian(ylim = c(-0.01, 1.01))
gg_rsf


## ----rfsrc-by-treatment, fig.cap="Median predicted survival by treatment group with 95% confidence bands."----
plot(gg_rfsrc(rfsrc_pbc, by = "treatment")) +
  theme(legend.position = c(0.2, 0.2)) +
  labs(y = "Survival Probability", x = "Time (years)") +
  coord_cartesian(ylim = c(-0.01, 1.01))


## ----predict-test, fig.cap="Predicted survival for non-trial patients (test set)."----
rfsrc_pbc_test <- predict(rfsrc_pbc, newdata = pbc_test,
                          na.action = "na.impute",
                          importance = TRUE)

plot(gg_rfsrc(rfsrc_pbc_test), alpha = 0.2) +
  scale_color_manual(values = event_colors) +
  theme(legend.position = "none") +
  labs(y = "Survival Probability", x = "Time (years)") +
  coord_cartesian(ylim = c(-0.01, 1.01))


## ----vimp-plot, fig.width=6, fig.cap="Variable importance ranking. Blue = positive VIMP, red = negative."----
plot(gg_vimp(rfsrc_pbc), lbls = st_labs) +
  theme(legend.position = c(0.8, 0.2)) +
  labs(fill = "VIMP > 0")


## ----varsel-------------------------------------------------------------------
md_pbc <- max.subtree(rfsrc_pbc)


## ----select-vars--------------------------------------------------------------
xvar     <- c("bili", "albumin", "copper", "prothrombin", "age")
xvar_cat <- "edema"
xvar_all <- c(xvar, xvar_cat)


## ----vardep-bili, fig.height=4, fig.cap="Variable dependence of survival at 1 and 3 years on bilirubin."----
gg_v <- gg_variable(rfsrc_pbc, time = c(1, 3),
                    time.labels = c("1 Year", "3 Years"))

plot(gg_v, xvar = "bili", alpha = 0.4) +
  labs(y = "Survival", x = st_labs["bili"]) +
  theme(legend.position = "none") +
  scale_color_manual(values = event_colors, labels = event_labels) +
  scale_shape_manual(values = event_marks, labels = event_labels) +
  coord_cartesian(ylim = c(-0.01, 1.01))


## ----vardep-panel, fig.height=5, fig.cap="Variable dependence at 1 and 3 years for continuous predictors."----
plot(gg_v, xvar = xvar[-1], panel = TRUE, alpha = 0.4) +
  labs(y = "Survival") +
  theme(legend.position = "none") +
  scale_color_manual(values = event_colors, labels = event_labels) +
  scale_shape_manual(values = event_marks, labels = event_labels) +
  coord_cartesian(ylim = c(-0.05, 1.05))


## ----vardep-categorical, fig.height=4, fig.cap="Variable dependence on edema (categorical)."----
plot(gg_v, xvar = xvar_cat, alpha = 0.4) +
  labs(y = "Survival") +
  theme(legend.position = "none") +
  scale_color_manual(values = event_colors, labels = event_labels) +
  scale_shape_manual(values = event_marks, labels = event_labels) +
  coord_cartesian(ylim = c(-0.01, 1.02))


## ----partial-dep, error=TRUE, fig.height=5, fig.cap="Partial dependence of survival at approximately 1 and 3 years on continuous predictors."----
try({
# partial.rfsrc() requires times that match the model's time.interest grid;
# gg_partial_rfsrc() snaps the requested values to the nearest observed times.
ti   <- rfsrc_pbc$time.interest
t1yr <- ti[which.min(abs(ti - 1))]
t3yr <- ti[which.min(abs(ti - 3))]

pd <- gg_partial_rfsrc(rfsrc_pbc, xvar.names = xvar,
                       partial.time = c(t1yr, t3yr))

# Quick S3 plot — survival forests produce time-series curves per predictor value
plot(pd)
})


## ----partial-dep-custom, error=TRUE, fig.height=5, fig.cap="Partial dependence (custom styling)."----
try({
ggplot(pd$continuous, aes(x = x, y = yhat,
                          color = factor(round(time, 2)),
                          group = factor(time))) +
  geom_line(linewidth = 1) +
  facet_wrap(~name, scales = "free_x") +
  labs(y = "Predicted Survival", x = "", color = "Year") +
  scale_color_manual(values = setNames(c("steelblue", "firebrick"),
                                       as.character(round(c(t1yr, t3yr), 2)))) +
  theme_bw()
})


## ----coplot-edema, fig.width=7, fig.height=4, fig.cap="Variable dependence on bilirubin, conditional on edema group."----
gg_v1 <- gg_variable(rfsrc_pbc, time = 1)
gg_v1$edema <- paste("edema =", gg_v1$edema)

plot(gg_v1, xvar = "bili", alpha = 0.5) +
  labs(y = "Survival at 1 Year", x = st_labs["bili"]) +
  theme(legend.position = "none") +
  scale_color_manual(values = event_colors, labels = event_labels) +
  scale_shape_manual(values = event_marks, labels = event_labels) +
  facet_grid(~edema) +
  coord_cartesian(ylim = c(-0.01, 1.01))


## ----coplot-albumin, fig.width=7, fig.height=5, fig.cap="Variable dependence on bilirubin, conditional on albumin groups."----
albumin_cts <- quantile_pts(gg_v1$albumin, groups = 6, intervals = TRUE)
gg_v1$albumin_grp <- cut(gg_v1$albumin, breaks = albumin_cts)
levels(gg_v1$albumin_grp) <- paste("albumin =",
                                    levels(gg_v1$albumin_grp))

plot(gg_v1, xvar = "bili", alpha = 0.5) +
  labs(y = "Survival at 1 Year", x = st_labs["bili"]) +
  theme(legend.position = "none") +
  scale_color_manual(values = event_colors, labels = event_labels) +
  scale_shape_manual(values = event_marks, labels = event_labels) +
  facet_wrap(~albumin_grp) +
  coord_cartesian(ylim = c(-0.01, 1.01))


## ----surface-data, error=TRUE, cache=FALSE------------------------------------
try({
# Create grid of albumin values
alb_grid <- quantile_pts(pbc_trial$albumin, groups = 25)

# For each albumin value, compute partial dependence on bili at ~1 year
surface_list <- lapply(alb_grid, function(alb_val) {
  newx <- pbc_trial[, rfsrc_pbc$xvar.names, drop = FALSE]
  newx$albumin <- alb_val
  pd_alb <- gg_partial_rfsrc(rfsrc_pbc, xvar.names = "bili",
                              newx = newx, partial.time = t1yr)
  df <- pd_alb$continuous
  df$albumin <- alb_val
  df
})

surface_df <- bind_rows(surface_list)
})


## ----plotly-surface, error=TRUE, fig.cap="Interactive partial dependence surface: survival as a function of bilirubin and albumin."----
try({
if (!exists("surface_df")) {
  message("surface_df not available — skipping plotly surface (see surface-data chunk error above).")
} else if (requireNamespace("plotly", quietly = TRUE)) {
  # Reshape for surface
  library(plotly)

  surface_wide <- surface_df |>
    select(bili = x, albumin, survival = yhat) |>
    arrange(albumin, bili)

  # Create matrix form
  bili_vals <- sort(unique(surface_wide$bili))
  alb_vals  <- sort(unique(surface_wide$albumin))
  z_matrix  <- matrix(surface_wide$survival,
                       nrow = length(alb_vals),
                       ncol = length(bili_vals),
                       byrow = TRUE)

  plot_ly(x = bili_vals, y = alb_vals, z = z_matrix) |>
    add_surface(colorscale = "Viridis", showscale = TRUE) |>
    layout(
      scene = list(
        xaxis = list(title = "Bilirubin"),
        yaxis = list(title = "Albumin"),
        zaxis = list(title = "Survival")
      )
    )
} else {
  message("Install the plotly package for interactive 3D surfaces.")
  # Fallback: contour plot with ggplot2
  ggplot(surface_df, aes(x = x, y = albumin, fill = yhat)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Survival") +
    labs(x = "Bilirubin", y = "Albumin") +
    theme_bw()
}
})


## ----brier-overall, fig.width=7, fig.height=3.5, fig.cap="Time-resolved Brier score for the PBC survival forest."----
gg_bs <- gg_brier(rfsrc_pbc)
plot(gg_bs)


## ----brier-envelope, fig.width=7, fig.height=3.5, fig.cap="Brier score with 15--85% per-subject envelope."----
plot(gg_bs, envelope = TRUE)


## ----brier-crps, fig.width=7, fig.height=3.5, fig.cap="Running CRPS for the PBC survival forest."----
plot(gg_bs, type = "crps")


## ----brier-crps-scalar--------------------------------------------------------
attr(gg_bs, "crps_integrated")

