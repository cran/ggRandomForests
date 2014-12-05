## ----setup, include = FALSE, cache = FALSE, echo = FALSE--------------------------------
library(knitr)
knitr::render_sweave() 
# set global chunk options for knitr. These can be changed in the header for each individual R code chunk
opts_chunk$set(fig.path = 'figure/rfs-', 
        fig.align = 'center', 
        fig.pos = "!htpb", 
        fig.show = 'hold', 
        fig.height = 3, 
        fig.width = 4, 
        size = 'footnotesize', 
        prompt = TRUE, 
        highlight = FALSE, 
        comment = NA, 
        echo = FALSE, results = FALSE, message = FALSE, warning = FALSE, 
        error = FALSE, dev = 'pdf', prompt = TRUE)

# Setup the R environment
options(object.size = Inf, expressions = 100000, memory = Inf, 
    replace.assign = TRUE, width = 90, prompt = "R> ")

#################
# Load_packages #
#################
library(ggplot2) # Graphics engine for generating all types of plots

library(dplyr) # Better data manipulations
library(tidyr)

library(ggRandomForests)

# Analysis packages.
library(randomForestSRC) 
library(RColorBrewer)

library(xtable)

options(mc.cores = 1, rf.cores = 1)

#########################################################################
# Default computation settings
#########################################################################
theme_set(theme_bw())
event.marks <- c(1, 4)
event.labels <- c(FALSE, TRUE)
strCol <- brewer.pal(3, "Set1")
strCol <- strCol[c(2, 1, 3)]
alpha <- .3
## Set the event/censor marks. Want open circle for censored, and 
## x for events 
event.marks = c(1, 4)
event.labels = c(FALSE, TRUE)
strCol <- brewer.pal(3, "Set1")
strCol <- c(strCol[2], strCol[1])


## ----datastep---------------------------------------------------------------------------
data(pbc, package = "randomForestSRC")

## Set modes correctly. For binary variables: transform to logical
## Check for range of 0, 1
## There is probably a better way to do this.
for(ind in 1:dim(pbc)[2]){
 if(!is.factor(pbc[, ind])){
  if(length(unique(pbc[which(!is.na(pbc[, ind])), ind]))<= 2) {
   if(sum(range(pbc[, ind], na.rm = TRUE) ==  c(0, 1)) ==  2){
    pbc[, ind] <- as.logical(pbc[, ind])
   }
  }
 }else{
  if(length(unique(pbc[which(!is.na(pbc[, ind])), ind]))<= 2) {
   if(sum(sort(unique(pbc[, ind])) ==  c(0, 1)) ==  2){
    pbc[, ind] <- as.logical(pbc[, ind])
   }
   if(sum(sort(unique(pbc[, ind])) ==  c(FALSE, TRUE)) ==  2){
    pbc[, ind] <- as.logical(pbc[, ind])
   }
  }
 }
 if(!is.logical(pbc[, ind]) & 
    length(unique(pbc[which(!is.na(pbc[, ind])), ind]))<= 5) {
  pbc[, ind] <- factor(pbc[, ind])
 }
}
# Convert age to years
pbc$age <- pbc$age/364.24
pbc$years <- pbc$days/364.24
pbc <- pbc %>% select(-days)

cls <- sapply(pbc, class) 

labels <- c("censoring indicator", "1 = D-penicillamine, 2 = placebo", 
"age in years", "0 = female, 1 = male", "presence of asictes", "presence of hepatomegaly", 
"presence of spiders", "presence of edema", "serum bilirubin in mg/dl", 
"serum cholesterol in mg/dl", "albumin in gm/dl", "urine copper in ug/day", 
"alkaline phosphatase in U/liter", "SGOT in U/ml", "triglicerides in mg/dl", 
"platelets per cubic ml/1000", "prothrombin time in seconds", "histologic stage of disease", 
"survival time in years")

dta.labs <- data.frame(cbind(names = colnames(pbc), label = labels, type = cls))

## ----dta-table, results = "asis"--------------------------------------------------------
rws <- seq(1, (nrow(dta.labs)), by = 2)
col <- rep("\\rowcolor[gray]{0.95}", length(rws))

print(xtable(dta.labs %>% select(-names), 
       caption = "PBC Data field descriptions", 
       label = "T:dataLabs", 
       digits = 3), 
   size = 'footnotesize', # fontsize
   booktabs = TRUE, 
   add.to.row = list(pos = as.list(rws), command = col), 
   command =  c('\\toprule' , 
          '\\midrule' , 
          '\\bottomrule ')
)

## ----rfsrc, echo = TRUE, eval = FALSE---------------------------------------------------
#  pbc_rf <- rfsrc(Surv(years, status) ~ ., data = pbc,
#           ntree = 2000,
#           na.action = "na.impute",
#           fast.restore = TRUE)

## ----read-forest, echo = FALSE, results = FALSE-----------------------------------------
data(pbc_rf, package="ggRandomForests")
pbc_rf

## ----errorPlot, fig.cap = "RSF prediction error estimates"------------------------------
ggerr <- gg_error(pbc_rf)
plot(ggerr)+
  coord_cartesian(ylim=c(.09,.31))

## ----rfsrc-plot, fig.cap = "PBC Survival"-----------------------------------------------
ggRFsrc <- plot.gg_rfsrc(pbc_rf, alpha = .2) + 
 scale_color_manual(values = strCol) + 
 theme(legend.position = "none") + 
 labs(y = "Survival Probability", x = "time (years)")
ggRFsrc

## ----rf-vimp, echo = TRUE, fig.cap = "Variable Importance"------------------------------
plot.gg_vimp(pbc_rf) + 
 theme(legend.position = "none")

## ----mindepth-view, eval = FALSE, echo = TRUE-------------------------------------------
#  pbc_vs <- var.select(pbc_rf)
#  ggMindepth <- gg_minimal_depth(pbc_vs)
#  print(ggMindepth)

## ----mindepth-load----------------------------------------------------------------------
data(pbc_vs, package="ggRandomForests")
ggMindepth <- gg_minimal_depth(pbc_vs)
ggMindepth

## ----mindepth-plot, echo=TRUE, fig.cap = "Minimal Depth Plot"---------------------------
plot(ggMindepth)

## ----rfsrc-plot3Mnth, echo = TRUE, fig.cap = "PBC Survival"-----------------------------
ggRFsrc + 
 geom_vline(aes(xintercept = c(1, 3)), linetype = "dashed") + 
 coord_cartesian(x = c(0, 4))

## ----variable-plotbili, echo = TRUE, fig.cap = "Variable dependence Survival vs. Bilirubin", fig.height = 4----
xvar <- pbc_vs$topvars[1:6]
ind = 1
ggrf <- gg_variable(pbc_rf, time = c(1, 3), time.labels = c("1 Year", "3 Years"))

plot(ggrf, x_var = xvar[ind], se=FALSE, alpha=.3) + 
 labs(y = "Survival") + 
 theme(legend.position = "none") + 
 scale_color_manual(values = strCol, labels = event.labels) + 
 scale_shape_manual(values = event.marks, labels = event.labels)

## ----variable-plotCombines, echo = TRUE, fig.cap = "Variable dependence Survival panel", fig.width = 7, fig.height = 4----
plot(ggrf, x_var = xvar[c(2,3,5,6)], panel = TRUE, 
     se=FALSE, alpha=.3, 
     method="glm", formula=y~poly(x,2)) + 
 labs(y = "Survival") + 
 theme(legend.position = "none") + 
 scale_color_manual(values = strCol, labels = event.labels) + 
 scale_shape_manual(values = event.marks, labels = event.labels)+
  coord_cartesian(y=c(1,102))

## ----pbc-partial, echo = TRUE, eval = FALSE---------------------------------------------
#  # Calculate the 1 year partial dependence
#  pbc_prtl <- plot.variable(pbc_rf, surv.type = "surv",
#               time = 364.25,
#             xvar.names = xvar, partial = TRUE,
#             show.plots = FALSE)
#  
#  # Calculate the 3 year partial dependence
#  pbc_prtl.3 <- plot.variable(pbc_rf, surv.type = "surv",
#                time = 3*364.25,
#             xvar.names = xvar, partial = TRUE,
#             show.plots = FALSE)
#  
#  # Create gg_partial objects
#  ggPrtl <- gg_partial(pbc_prtl)
#  ggPrtl.3 <- gg_partial(pbc_prtl.3)
#  
#  # Combine the objects to get multiple time curves
#  # along variables on a single figure.
#  pbc_ggpart <- combine(ggPrtl, ggPrtl.3,
#             labels = c("1 Year", "3 Years"))
#  

## ----pbc-partial-load-------------------------------------------------------------------
data("pbc_prtl", package="ggRandomForests")
data("pbc_ggpart", package="ggRandomForests")

## ----pbc-partial-bili, echo = TRUE, fig.cap = "Risk adjusted Survival"------------------
plot(pbc_ggpart[["bili"]], se = FALSE) + 
 theme(legend.position = c(.8, .5)) + 
 labs(y = "Survival", 
    x = dta.labs[which(rownames(dta.labs) ==  "bili"), "label"])

## ----pbc-partial-panel, echo = TRUE, fig.cap = "Risk adjusted Survival - panel plot", fig.width = 7, fig.height = 4----
pbc_ggpart$bili <- pbc_ggpart$edema <- NULL
plot(pbc_ggpart, se = FALSE, panel = TRUE) + 
  theme(legend.position=c(.6,.6))+
 labs(x = "", y = "Survival")

## ----interaction-show, echo = TRUE, eval = FALSE----------------------------------------
#  pbc_interaction <- find.interaction(pbc_rf)
#  
#  ggint <- gg_interaction(pbc_interaction)
#  plot(ggint, x_var = "bili") +
#   labs(y = "Interactive Minimal Depth")

## ----interaction, fig.cap = "Minimal Depth interaction for Surgical Date"---------------
 #pDat.int <- find.interaction(pbc_rf)
data("pbc_interaction", package="ggRandomForests")

plot(gg_interaction(pbc_interaction), x_var = "bili") + 
 labs(y = "Interactive Minimal Depth")

## ----interactionPanel, echo = TRUE, fig.cap = "Risk adjusted Survival - panel plot", fig.width = 7, fig.height = 4----
plot(gg_interaction(pbc_interaction), x_var = xvar[2:5]) + 
 labs(y = "Interactive Minimal Depth") + 
 theme(legend.position = "none")

## ----var_dep, echo = TRUE, fig.cap = "Bilirubin Variable Dependence at 1 year."---------
ggrf <- gg_variable(pbc_rf, time = 1)

ggvar <- ggrf
ggvar$treatment <- as.numeric(ggvar$treatment)
ggvar$treatment[which(ggvar$treatment==1)] <- "D-pen" 
ggvar$treatment[which(ggvar$treatment==2)] <- "placebo" 
ggvar$treatment <- factor(ggvar$treatment)

ggvar$stage <- paste("stage=", ggvar$stage, sep="")

var_dep <- plot(ggvar, x_var = "bili", smooth = TRUE, 
        method = "glm", alpha = .5, se = FALSE) + 
 labs(y = "Survival", x = dta.labs[which(rownames(dta.labs) ==  xvar[ind]), "label"]) + 
 theme(legend.position = "none") + 
 scale_color_manual(values = strCol, labels = event.labels) + 
 scale_shape_manual(values = event.marks, labels = event.labels)

show(var_dep)

## ----coplot_bilirubin, echo = TRUE, fig.cap = "Conditional Variable Dependence. Interactions between bilrubin with treatment and stage variables.", fig.width = 7, fig.height = 4----
var_dep + 
 facet_grid(treatment~stage)

