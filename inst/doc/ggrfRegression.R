## ----setup, include = FALSE, cache = FALSE, echo=FALSE----------------------------------
## Not displayed ##
library("knitr")

# set global chunk options for knitr. These can be changed in the header for 
# each individual R code chunk
opts_chunk$set(#fig.path = "regression/",
               fig.align = "center",
               fig.pos = "!htb",
               fig.show = "hold",
               # fig.height = 3,
               # fig.width = 4,
               size = "footnotesize",
               prompt = TRUE,
               highlight = FALSE,
               comment = NA,
               
               # Change echo to TRUE if you want to see all the code examples
               echo = FALSE, 
               
               results = FALSE,
               message = FALSE,
               warning = FALSE,
               error = FALSE)

# Setup the R environment
options(object.size = Inf, expressions = 100000, memory = Inf, 
        replace.assign = TRUE, width = 90)
options(mc.cores = 1, rf.cores = -1)

## ----libraries, echo=TRUE---------------------------------------------------------------

################## Load packages ##################
library("ggplot2")         # Graphics engine
library("RColorBrewer")    # Nice color palettes
library("plot3D")          # for 3d surfaces. 
library("dplyr")           # Better data manipulations
library("tidyr")           # gather variables into long format
library("parallel")        # mclapply for multicore processing

# Analysis packages.
library("randomForestSRC") # random forests for survival, regression and 
# classification
library("ggRandomForests") # ggplot2 random forest figures (This!)

################ Default Settings ##################
theme_set(theme_bw())     # A ggplot2 theme with white background

## Set open circle for censored, and x for events 
event_marks <- c(1, 4)
event_labels <- c(FALSE, TRUE)

## We want red for death events, so reorder this set.
str_col <- brewer.pal(3, "Set1")[c(2, 1, 3)]


## ----vignette, eval=FALSE, echo=TRUE----------------------------------------------------
#  vignette("randomForestSRC-Regression", package = "ggRFVignette")

## ----datastep, echo=TRUE----------------------------------------------------------------
# Load the Boston Housing data
data(Boston, package = "MASS")

# Set modes correctly. For binary variables: transform to logical
Boston$chas <- as.logical(Boston$chas)

## ----cleanup, results="asis", echo=FALSE------------------------------------------------
cls <- sapply(Boston, class) 
# 
lbls <- 
  #crim
  c("Crime rate by town.",
    # zn
    "Proportion of residential land zoned for lots over 25,000 sq.ft.",
    # indus
    "Proportion of non-retail business acres per town.",
    # chas
    "Charles River (tract bounds river).",
    # nox
    "Nitrogen oxides concentration (10 ppm).",
    # rm
    "Number of rooms per dwelling.",
    # age
    "Proportion of units built prior to 1940.",
    # dis
    "Distances to Boston employment center.",
    # rad
    "Accessibility to highways.",
    # tax
    "Property tax rate per $10,000.",
    # ptratio
    "Pupil teacher ratio by town.",
    # black
    "Proportion of blacks by town.",
    # lstat
    "Lower status of the population (percent).",
    # medv
    "Median value of homes ($1000s).")

# Build a table for data description
dta_labs <- data.frame(cbind(Variable = names(cls), 
                            Description = lbls, 
                            type = cls))

# Build a named vector for labeling figures later/
st_labs <- as.character(dta_labs$Description)
names(st_labs) <- names(cls)

# Print the descriptive table.
kable(dta_labs, 
      row.names = FALSE, 
      caption  = "`Boston` housing data dictionary.",
      booktabs = FALSE)

## ----eda, fig.cap="EDA variable plots. Points indicate variable value against the median home value variable. Points are colored according to the chas variable.", fig.width=7, fig.height=5----
# Use tidyr::gather to transform the data into long format.
dta <- gather(Boston, variable, value, -medv, -chas)

# plot panels for each covariate colored by the logical chas variable.
ggplot(dta) +
  geom_point(alpha = 0.4, aes(x = medv, y = value, color = chas)) +
  geom_smooth(aes(x = medv, y = value), se = FALSE) + 
  labs(y = "", x = st_labs["medv"]) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~ variable, scales = "free_y", ncol = 3)

## ----randomforest, echo=TRUE------------------------------------------------------------
# Load the data, from the call:
rfsrc_boston <- rfsrc(medv ~ ., data = Boston, 
                      importance = TRUE, err.block = 5)

# print the forest summary
rfsrc_boston

## ----error, echo=TRUE, fig.cap="Random forest generalization error. OOB error convergence along the number of trees in the forest."----
# Plot the OOB errors against the growth of the forest.
gg_e <- gg_error(rfsrc_boston)
gg_e <- gg_e %>% filter(!is.na(error))
class(gg_e) <- c("gg_error", class(gg_e))
plot(gg_e)

## ----rfsrc, echo=TRUE, fig.cap="OOB predicted median home values. Points are jittered to help visualize predictions for each observation. Boxplot indicates the distribution of the predicted values."----
# Plot predicted median home values.
plot(gg_rfsrc(rfsrc_boston), alpha = .5) +
  coord_cartesian(ylim = c(5, 49))

## ----vimp, echo=TRUE, fig.cap="Random forest VIMP plot. Bars are colored by sign of VIMP, longer blue bars indicate more important variables.", fig.width=7, fig.height=5----
# Plot the VIMP rankings of independent variables.
plot(gg_vimp(rfsrc_boston), lbls = st_labs)

## ----minimaldepth, echo=TRUE, fig.cap="Minimal Depth variables in rank order, most important at the top. Vertical dashed line indicates the maximal minimal depth for important variables.", fig.width=7, fig.height=5----
# Load the data, from the call:
varsel_boston <- var.select(rfsrc_boston)

# Save the gg_minimal_depth object for later use.
gg_md <- gg_minimal_depth(varsel_boston)

# plot the object
plot(gg_md, lbls = st_labs)

## ----minimalvimp, echo=TRUE, fig.cap="Comparing Minimal Depth and Vimp rankings. Points on the red dashed line are ranked equivalently, points below have higher VIMP, those above have higher minimal depth ranking. Variables are colored by the sign of the VIMP measure."----
# gg_minimal_depth objects contain information about
# both minimal depth and VIMP.
plot(gg_minimal_vimp(gg_md))

## ----variable, echo=FALSE, fig.cap="Variable dependence plot. Individual case predictions are marked with points. Loess smooth curve indicates the trend as the variables increase with shaded 95\\% confidence band.", fig.width=7, fig.height=5----
# Create the variable dependence object from the random forest
gg_v <- gg_variable(rfsrc_boston)

# We want the top ranked minimal depth variables only,
# plotted in minimal depth rank order. 
xvar <- gg_md$topvars

# plot the variable list in a single panel plot
plot(gg_v, xvar = xvar, panel = TRUE, alpha = .5) +
  labs(y = st_labs["medv"], x = "")

## ----chas, echo=TRUE, fig.cap="Variable dependence for Charles River logical variable."----
plot(gg_v, xvar = "chas", alpha = 0.4) +
  labs(y = st_labs["medv"])

# , points=FALSE, se=FALSE, notch=TRUE

## ----partial, echo=FALSE, fig.cap="Partial dependence panels. Risk adjusted variable dependence for variables in minimal depth rank order.", fig.width=7, fig.height=5, cache=TRUE----

# Load the data, from the call:
partial_boston <- plot.variable(rfsrc_boston,
                                xvar = gg_md$topvars,
                                partial = TRUE, sorted = FALSE,
                                show.plots = FALSE)

# generate a list of gg_partial objects, one per xvar.
gg_p <- gg_partial(partial_boston)

# plot the variable list in a single panel plot
plot(gg_p, panel = TRUE) +
  labs(y = st_labs["medv"], x = "") +
  geom_smooth(method = "loess", formula = "y ~ x", se = FALSE)

## ----interactions, echo=FALSE, fig.cap="Minimal depth variable interactions. Reference variables are marked with red cross in each panel. Higher values indicate lower interactivity with reference variable.", fig.width=7, fig.height=5, cache=TRUE----
# Load the data, from the call:
interaction_boston <- find.interaction(rfsrc_boston)

# Plot the results in a single panel.
plot(gg_interaction(interaction_boston), 
     xvar = gg_md$topvars, panel = TRUE)

## ----coplots, echo=FALSE, fig.cap="Variable Coplots. Predicted median home values as a function of percentage of lower status population, stratified by average number of rooms groups.", fig.width=7, fig.height=5----
# Find the rm variable points to create 6 intervals of roughly 
# equal size population
rm_pts <- quantile_pts(rfsrc_boston$xvar$rm, groups = 6, 
                       intervals = TRUE)

# Pass these variable points to create the 6 (factor) intervals
rm_grp <- cut(rfsrc_boston$xvar$rm, breaks = rm_pts)

# Append the group factor to the gg_variable object
gg_v$rm_grp <- rm_grp

# Modify the labels for descriptive panel titles 
levels(gg_v$rm_grp) <- paste("rm in ", 
                             levels(gg_v$rm_grp), sep = "")

# Create a variable dependence (co)plot, 
# faceted on group membership.
plot(gg_v, xvar = "lstat") +
  #   method = "loess", span=1.5, se = FALSE) + 
  labs(y = st_labs["medv"], x = st_labs["lstat"]) + 
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Set3") + 
  # geom_smooth(method = "loess", formula = "y ~ x", se=FALSE) +
  facet_wrap(~rm_grp)

## ----coplots2, echo=FALSE, fig.cap="Variable Coplots. Predicted median home value as a function of average number of rooms, stratified by percentage of lower status groups.", fig.width=7, fig.height=5----
# Find the lstat variable points to create 6 intervals of roughly 
# equal size population
lstat_pts <- quantile_pts(rfsrc_boston$xvar$lstat, 
                          groups = 6, intervals = TRUE)

# Pass these variable points to create the 6 (factor) intervals
lstat_grp <- cut(rfsrc_boston$xvar$lstat, breaks = lstat_pts)

# Append the group factor to the gg_variable object
gg_v$lstat_grp <- lstat_grp

# Modify the labels for descriptive panel titles 
levels(gg_v$lstat_grp) <- paste("lstat in ", 
                                levels(gg_v$lstat_grp), 
                                " (%)", sep = "")

# Create a variable dependence (co)plot, faceted on group membership.
plot(gg_v, xvar = "rm", alpha = .5) +
  #method = "loess", span = 1.5, , se = FALSE) + 
  labs(y = st_labs["medv"], x = st_labs["rm"]) + 
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Set3") + 
  # geom_smooth() +
  # scale_shape_manual(values = event_marks, labels = event_labels)+ 
  facet_wrap(~lstat_grp)

## ----prtl-copl, echo=TRUE, cache=TRUE---------------------------------------------------
partial_coplot_boston <- gg_partial_coplot(rfsrc_boston, 
                                           xvar = "lstat", 
                                           groups = rm_grp)

## ----prtl-coplots, echo=TRUE, fig.cap="Partial Coplots. Risk adjusted predicted median value as a function of Lower Status, conditional on groups of average number of rooms.", fig.width=7, fig.height=5----
# Load the stored partial coplot data.

# # Partial coplot
ggplot(partial_coplot_boston, aes(x = lstat, y = yhat, col = group)) +
  geom_point() +
  geom_smooth(method = "loess", formula = "y ~ x",
              se = FALSE, alpha = 0.25) +
  labs(x = st_labs["lstat"], y = st_labs["medv"],
       color = "Room", shape = "Room") +
  scale_color_brewer(palette = "Set1")

## ----prtl-copl2, echo=TRUE, cache=TRUE--------------------------------------------------
partial_coplot_boston2 <- gg_partial_coplot(rfsrc_boston, 
                                            xvar = "rm", 
                                            groups = lstat_grp)

## ----prtl-coplots2, echo=TRUE, fig.cap="Partial Coplots. Risk adjusted predicted median value as a function of average number of rooms, conditional on groups of percentage of lower status population.", fig.width=7, fig.height=5----
# Partial coplot
#plot(partial_coplot_boston2)+ ## again plot.gg_partial_coplot
ggplot(partial_coplot_boston2, aes(x = rm, y = yhat, col = group)) +
  geom_point() +
  geom_smooth(method = "loess", formula = "y ~ x", se = FALSE) +
  labs(x = st_labs["rm"], y = st_labs["medv"], 
       color = "Lower Status", shape = "Lower Status") +
  scale_color_brewer(palette = "Set1")

## ----def-pts, echo=TRUE-----------------------------------------------------------------
# Find the quantile points to create 50 cut points
rm_pts <- quantile_pts(rfsrc_boston$xvar$rm, 
                       groups = 49, intervals = TRUE)

## ----prtl-surface, echo=TRUE, cache=TRUE------------------------------------------------
# Generate the gg_partial_coplot data object
system.time(partial_boston_surf <- lapply(rm_pts, function(ct) {
  rfsrc_boston$xvar$rm <- ct
  plot.variable(rfsrc_boston, xvar = "lstat", time = 1,
                npts = 50, show.plots = TRUE, 
                partial = TRUE)
}))
#     user   system  elapsed 
#    49.702   0.776  51.006 

## ----contour3d, echo=FALSE, fig.cap="Partial coplot contour plot. Contours of median home value along the lstat/rm plane.", fig.width=7, fig.height=5----

# Instead of groups, we want the raw rm point values,
# To make the dimensions match, we need to repeat the values
# for each of the 50 points in the lstat direction
rm_tmp <- do.call(c, lapply(rm_pts, 
                           function(grp) rep(grp, 50)))

# Convert the list of plot.variable output to 
partial_surf <- do.call(rbind, lapply(partial_boston_surf, gg_partial))

# attach the data to the gg_partial_coplot
partial_surf$rm <- rm_tmp

# ggplot2 contour plot of x, y and z data.
ggplot(partial_surf, 
       aes(x = lstat, y = rm, z = yhat, colour = ..level..)) +
  stat_contour() +
  labs(x = st_labs["lstat"], y = st_labs["rm"], 
       color = "Median Home Values") +
  scale_colour_gradientn(colours = topo.colors(10))

## ----surface3d, echo=FALSE, fig.cap="Partial plot surface.", fig.width=7, fig.height=5----
# Modify the figure margins to make the figure larger
par(mai = c(0, 0, 0, 0))

# Transform the gg_partial_coplot object into a list of three named matrices
# for surface plotting with plot3D::surf3D
suppressWarnings(
  srf <- surface_matrix(partial_surf, c("lstat", "rm", "yhat"))
)

# Generate the figure.
surf3D(x = srf$x, y = srf$y, z = srf$z, col = topo.colors(10),
       colkey = FALSE, border = "black", bty = "b2", 
       shade = 0.5, expand = 0.5, 
       lighting = TRUE, lphi = -50,
       xlab = "Lower Status", 
       ylab = "Average Rooms", 
       zlab = "Median Value"
)

