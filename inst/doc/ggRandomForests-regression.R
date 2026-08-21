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
library(tidyr)
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
data(Boston, package = "MASS")
Boston$chas <- as.logical(Boston$chas) # nolint: object_name_linter


## ----variable-labels----------------------------------------------------------
st_labs <- c(
  crim    = "Crime rate by town",
  zn      = "Residential land zoned > 25k sq ft (%)",
  indus   = "Non-retail business acres (%)",
  chas    = "Borders Charles River",
  nox     = "Nitrogen oxides (10 ppm)",
  rm      = "Rooms per dwelling",
  age     = "Units built before 1940 (%)",
  dis     = "Distance to employment centers",
  rad     = "Highway accessibility index",
  tax     = "Property tax rate per $10,000",
  ptratio = "Pupil-teacher ratio",
  black   = "Proportion of Black residents",
  lstat   = "Lower status population (%)",
  medv    = "Median home value ($1000s)"
)


## ----eda, fig.cap="EDA: each predictor vs. median home value, colored by Charles River indicator."----
dta <- Boston |>
  pivot_longer(c(-medv, -chas), names_to = "variable", values_to = "value")

ggplot(dta, aes(x = medv, y = value, color = chas)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = medv, y = value), color = "gray30",
              inherit.aes = FALSE, se = FALSE) +
  labs(y = "", x = st_labs["medv"]) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~variable, scales = "free_y", ncol = 3)


## ----rfsrc-fit----------------------------------------------------------------
rfsrc_Boston <- rfsrc(medv ~ ., data = Boston, # nolint: object_name_linter
                      ntree = 100, importance = TRUE, err.block = 5)
rfsrc_Boston


## ----error-plot, fig.height=4, fig.cap="OOB mean squared error vs. number of trees."----
gg_e <- gg_error(rfsrc_Boston)
gg_e <- gg_e |> filter(!is.na(error))
class(gg_e) <- c("gg_error", class(gg_e))
plot(gg_e)


## ----rfsrc-predicted, fig.height=4, fig.cap="OOB predicted median home values. Points are jittered; boxplot shows the distribution."----
plot(gg_rfsrc(rfsrc_Boston), alpha = 0.5) +
  coord_cartesian(ylim = c(5, 49))


## ----vimp-plot, fig.cap="VIMP ranking. Longer blue bars indicate more important variables."----
plot(gg_vimp(rfsrc_Boston), lbls = st_labs)


## ----varsel-------------------------------------------------------------------
md_Boston <- max.subtree(rfsrc_Boston) # nolint: object_name_linter


## ----select-vars--------------------------------------------------------------
xvar <- md_Boston$topvars


## ----shap-fit-----------------------------------------------------------------
set.seed(42)
shap_sample <- Boston[sample(nrow(Boston), 25), setdiff(colnames(Boston), "medv")]
gg_shp <- gg_shap(rfsrc_Boston, newdata = shap_sample, bg_n = 30)


## ----shap-importance-plot, fig.cap="Mean absolute SHAP value per predictor."----
plot(gg_shp, type = "importance")


## ----shap-beeswarm-plot, fig.cap="SHAP beeswarm: every dot is one tract's contribution for one predictor."----
plot(gg_shp, type = "beeswarm")


## ----shap-dependence-plot, fig.cap="SHAP dependence for lstat."---------------
plot(gg_shp, type = "dependence", xvar = "lstat")


## ----vardep-panel, fig.cap="Variable dependence for top predictors (minimal depth rank order)."----
gg_v <- gg_variable(rfsrc_Boston)

plot(gg_v, xvar = xvar, panel = TRUE, alpha = 0.5) +
  labs(y = st_labs["medv"], x = "")


## ----vardep-chas, fig.height=4, fig.cap="Variable dependence for Charles River (categorical)."----
plot(gg_v, xvar = "chas", alpha = 0.4) +
  labs(y = st_labs["medv"])


## ----partial-dep, fig.cap="Partial dependence for top predictors."------------
pd <- gg_partial_rfsrc(rfsrc_Boston, xvar.names = xvar)

# Quick S3 plot — works out of the box for the standard regression case
plot(pd)


## ----partial-dep-custom, fig.cap="Partial dependence (custom styling)."-------
ggplot(pd$continuous, aes(x = x, y = yhat)) +
  geom_line(color = "steelblue", linewidth = 1) +
  facet_wrap(~name, scales = "free_x") +
  labs(y = st_labs["medv"], x = "") +
  theme_bw()


## ----coplot-chas, fig.height=4, fig.cap="Variable dependence on lstat, conditional on Charles River."----
gg_v$chas_label <- ifelse(gg_v$chas, "Borders Charles River",
                          "Does not border")

plot(gg_v, xvar = "lstat", alpha = 0.5) +
  labs(y = st_labs["medv"], x = st_labs["lstat"]) +
  theme(legend.position = "none") +
  facet_wrap(~chas_label)


## ----coplot-rm, fig.cap="Predicted medv vs. lstat, conditional on rooms-per-dwelling groups."----
rm_pts <- quantile_pts(rfsrc_Boston$xvar$rm, groups = 6, intervals = TRUE)
gg_v$rm_grp <- cut(rfsrc_Boston$xvar$rm, breaks = rm_pts)
levels(gg_v$rm_grp) <- paste("rm in", levels(gg_v$rm_grp))

plot(gg_v, xvar = "lstat", alpha = 0.5) +
  labs(y = st_labs["medv"], x = st_labs["lstat"]) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set3") +
  facet_wrap(~rm_grp)


## ----coplot-lstat, fig.cap="Predicted medv vs. rooms, conditional on lower-status groups."----
lstat_pts <- quantile_pts(rfsrc_Boston$xvar$lstat, groups = 6,
                          intervals = TRUE)
gg_v$lstat_grp <- cut(rfsrc_Boston$xvar$lstat, breaks = lstat_pts)
levels(gg_v$lstat_grp) <- paste("lstat in", levels(gg_v$lstat_grp))

plot(gg_v, xvar = "rm", alpha = 0.5) +
  labs(y = st_labs["medv"], x = st_labs["rm"]) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set3") +
  facet_wrap(~lstat_grp)


## ----surface-data, cache=TRUE-------------------------------------------------
rm_grid <- quantile_pts(rfsrc_Boston$xvar$rm, groups = 6)

surface_list <- lapply(rm_grid, function(rm_val) {
  newx <- rfsrc_Boston$xvar
  newx$rm <- rm_val
  pd_rm <- gg_partial_rfsrc(rfsrc_Boston, xvar.names = "lstat", newx = newx)
  df <- pd_rm$continuous
  df$rm <- rm_val
  df
})

surface_df <- bind_rows(surface_list)


## ----pd-surface, fig.cap="Partial dependence surface: median home value as a function of lstat and rm. Fill color is the predicted median value."----
ggplot(surface_df, aes(x = x, y = rm, fill = yhat)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Median Value\n($1000s)") +
  labs(x = "Lower Status (%)", y = "Rooms per Dwelling") +
  theme_bw()

