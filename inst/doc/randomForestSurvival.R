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
library(parallel)

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
pbc$treatment <- as.numeric(pbc$treatment)
pbc$treatment[which(pbc$treatment==1)] <- "DPCA"
pbc$treatment[which(pbc$treatment==2)] <- "placebo"
pbc$treatment <- factor(pbc$treatment)

cls <- sapply(pbc, class) 

labels <- c("censoring indicator", 
            "Treament", 
            "age in years", 
            "Female", 
            "Asictes", 
            "Hepatomegaly", 
            "Spiders", 
            "Edema", 
            "serum bilirubin (mg/dl)", 
            "serum cholesterol (mg/dl)", 
            "albumin (gm/dl)", 
            "urine copper (ug/day)", 
            "alkaline phosphatase (U/liter)", 
            "SGOT (U/ml)", 
            "triglicerides (mg/dl)", 
            "platelets per cubic ml/1000", 
            "prothrombin time (sec)", 
            "histologic stage", 
            "survival time (years)")

dta.labs <- data.frame(cbind(names = colnames(pbc), label = labels, type = cls))

st.labs <- as.character(dta.labs$label)
names(st.labs) <- rownames(dta.labs)

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


## ----gg_survival, fig.cap="Kaplan-Meier pbc data survival estimates comparing the treatment with placebo. Mean survival with shaded 95\\% condfidence band.", echo=TRUE----
gg_dta <- gg_survival(interval="years",censor="status", 
                      strat="treatment", data=pbc )
plot(gg_dta, se=.95) +
  labs(y="Survival Probability", x="Observation Time (years)", 
       color="Treatment", fill="Treatment")+
  theme(legend.position=c(.2,.2))

## ----gg_survival-bili, fig.cap="Kaplan-Meier pbc data survival estimates comparing Bilirubin measures.", echo=TRUE, fig.width=5----
pbc.alt <- pbc
pbc.alt$bili_grp <- cut(pbc.alt$bili, 
                        breaks=c(0,.8,1.3,3.4,max(pbc.alt$bili)))

gg_dta <- gg_survival(interval="years",censor="status", 
                      strat="bili_grp", data=pbc.alt )

plot(gg_dta, error="none") +
  labs(y="Survival Probability", 
       x="Observation Time (years)", 
       color="Bilirubin")

## ----xtab, results="asis"---------------------------------------------------------------
fleming.table <- data.frame(matrix(ncol=3, nrow=5))
rownames(fleming.table) <- c("Age", "log(Albumin)", "log(Bilirubin", "Edema", "log(Prothrombin Time)")
colnames(fleming.table) <- c("Coef.", "Std. Err.", "Z stat.")
fleming.table[,1] <- c(0.0333, -3.0553,0.8792, 0.7847, 3.0157) 
fleming.table[,2] <- c(0.00866, 0.72408,0.09873,0.29913,1.02380) 
fleming.table[,3] <- c(3.84,-4.22,8.9,2.62,2.95) 

rws <- seq(1, (nrow(fleming.table)), by = 2)
col <- rep("\\rowcolor[gray]{0.95}", length(rws))

print(xtable(fleming.table, 
             caption = "Regression model with log transformations of continuous variables, 312 randomized cases with PBC.", 
             label = "T:FHmodel", 
             digits = 4), 
      size = 'footnotesize', # fontsize
      booktabs = TRUE, 
      add.to.row = list(pos = as.list(rws), command = col), 
      command =  c('\\toprule' , 
                   '\\midrule' , 
                   '\\bottomrule ')
      )

## ----rfsrc, echo = TRUE, eval = FALSE---------------------------------------------------
#  pbc_rf <- rfsrc(Surv(years, status) ~ ., data = pbc,
#                  nsplit = 10,
#                  na.action = "na.impute")
#  pbc_rf

## ----read-forest, echo = FALSE, results = FALSE-----------------------------------------
data(pbc_rf, package="ggRandomForests")
pbc_rf

## ----errorPlot, fig.cap = "Random forest prediction error estimates as a function of the number of trees in the forest."----
ggerr <- gg_error(pbc_rf)
plot(ggerr)+
  coord_cartesian(ylim=c(.09,.31))

## ----rfsrc-plot, fig.cap = "Random forest predicted survival. Blue lines correspond to censored observations, red lines correspond to patients who experienced the event (death)."----
ggRFsrc <- plot.gg_rfsrc(pbc_rf, alpha=.2) + 
  scale_color_manual(values = strCol) + 
  theme(legend.position = "none") + 
  labs(y = "Survival Probability", x = "time (years)")+
  coord_cartesian(ylim=c(-.01,1.01))

ggRFsrc

## ----rfsrc-mean, fig.cap = "Mean value random forest predicted survival with shaded 95\\% confidence band."----
gg_src <- plot.gg_rfsrc(pbc_rf,  se=.95) + 
  scale_color_manual(values = strCol) + 
  theme(legend.position = "none") + 
  labs(y = "Survival Probability", x = "time (years)")+
  coord_cartesian(ylim=c(-.01,1.01))
gg_src

## ----rf-vimp, echo = TRUE, fig.cap = "Random forest variable Importance (VIMP). Blue bars indicate important variables (positive VIMP), red indicates noise variables (negative VIMP).", fig.width=5----
plot.gg_vimp(pbc_rf, lbls = st.labs) + 
  theme(legend.position = c(.8,.2))+
  labs(fill="VIMP > 0")+
  scale_fill_brewer(palette="Set1")

## ----mindepth-view, eval = FALSE, echo = TRUE-------------------------------------------
#  pbc_vs <- var.select(pbc_rf)
#  ggMindepth <- gg_minimal_depth(pbc_vs, lbls = st.labs)
#  print(ggMindepth)

## ----mindepth-load----------------------------------------------------------------------
data(pbc_vs, package="ggRandomForests")
ggMindepth <- gg_minimal_depth(pbc_vs)
ggMindepth

## ----mindepth-plot, echo=TRUE, fig.cap = "Minimal Depth variable selection. Low minimal depth indicates important variables. The dashed line is the threshold of maximum value for variable selection.", fig.width=5----
plot(ggMindepth, lbls = st.labs)

## ----depthVimp, fig.cap="Comparing Minimal Depth and Vimp rankings. Points on the red dashed line are ranked equivalently, points below have higher VIMP, those above have higher minimal depth ranking. Variables are colored by the sign of the VIMP measure.", fig.width=6----
gg.mv <- gg_minimal_vimp(pbc_vs)
plot(gg.mv, lbls = st.labs)+
  scale_y_continuous(breaks=seq(0,20,2))

## ----rfsrc-plot3Mnth, echo = TRUE, fig.cap = "Random forest OOB predicted patient survival. Red curves correspond to patients which have died, blue corresponds to alive (or censored) cases. Vertical dashed lines indicate the 1 and 3 year survival estimates."----
ggRFsrc + 
  geom_vline(aes(xintercept = c(1, 3)), linetype = "dashed") + 
  coord_cartesian(x = c(0, 4))

## ----variable-plotbili, echo = TRUE, fig.cap = "Bilirubin variable dependence at 1 and 3 years. Individual cases are marked with blue circles (alive or censored) and red xs (dead). Loess smooth curve with shaded 95\\% confidence band indicates the survival trend with increasing bilirubin.", fig.height = 4----
xvar <- pbc_vs$topvars[1:6]
ind = 1
ggrf <- gg_variable(pbc_rf, time = c(1, 3), 
                    time.labels = c("1 Year", "3 Years"))

plot(ggrf, x_var = xvar[ind], se=.95, alpha=.3) + 
  labs(y = "Survival") + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)

## ----variable-plotCombines, echo = TRUE, fig.cap = "Variable dependence plots at 1 and 3 years for continuous variables age, albumin, copper and prothrombin. Individual cases are marked with blue circles (alive or censored) and red xs (dead). Loess smooth curve indicates the survival trend with increasing variable value.", fig.width = 7, fig.height = 4----
plot(ggrf, x_var = xvar[c(2,3,5,6)], panel = TRUE, 
     se=FALSE, alpha=.3, 
     method="glm", formula=y~poly(x,2)) + 
  labs(y = "Survival") + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)+
  coord_cartesian(y=c(1,102))

## ----variable-plotEdema, echo = TRUE, fig.cap = "Variable dependence plots at 1 and 3 years for categorical edema variable. Individual cases are marked with blue circles (alive or censored) and red xs (dead). Boxes indicate distributional properties of observations in each group."----
plot(ggrf, x_var = "edema", notch=TRUE) + 
  labs(y = "Survival") + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)+
  coord_cartesian(y=c(1,102))

## ----pbc-partial, echo = TRUE, eval = FALSE---------------------------------------------
#  # Calculate the 1, 3 and 5 year partial dependence
#  pbc_prtl_time <- lapply(c(1,3,5), function(tm){
#    plot.variable(pbc_rf, surv.type = "surv",
#                            time = tm,
#                            xvar.names = xvar, partial = TRUE,
#                            show.plots = FALSE)
#  })
#  

## ----pbc-partial-load-------------------------------------------------------------------
data("pbc_prtl_time", package="ggRandomForests")

## ----pbc-partial-bili, echo = TRUE, fig.cap = "Partial dependence plot of (risk adjusted) predicted survival probability as a function of serum bilirubin at 1 year (red circle) and 3 years (blue triangle). Loess smooth curves indicates the trend."----
# Convert all partial plots to gg_partial objects
gg_dta <- lapply(pbc_prtl_time, gg_partial)

# Combine the objects to get multiple time curves 
# along variables on a single figure.
pbc_ggpart <- combine(gg_dta[[1]],gg_dta[[2]], 
                      lbls = c("1 Year", "3 Years"))

plot(pbc_ggpart[["bili"]], se = FALSE) + 
  theme(legend.position = c(.8, .5)) + 
  labs(y = "Survival", 
       x = dta.labs["bili", "label"],
       color="Time", shape="Time")+
  scale_color_brewer(palette="Set1")

## ----pbc-partial-panel, echo = TRUE, fig.cap = "Partial dependence plot of (risk adjusted) predicted survival probability as a function continuous variables prothrombin time, albumin, age and urin copper at 1 year (red circle) and 3 years (blue triangle).", fig.width = 7, fig.height = 4----
ggpart <- pbc_ggpart
ggpart$bili <- ggpart$edema <- NULL
plot(ggpart, se = FALSE, panel = TRUE) + 
  labs(x = "", y = "Survival") +
  scale_color_brewer(palette="Set1")

## ----pbc-partial-edema, echo = TRUE, fig.cap = "Partial dependence plot of (risk adjusted) predicted survival probability as a function of edema (categorical variable) at 1 year (red circle) and 3 years (blue triangle). Points indicate risk adjusted prediction for all patients within each edema group. Box plots indicate distributional properties within each group."----
plot(pbc_ggpart$edema, notch=TRUE, alpha=.3, outlier.shape = NA) + 
  labs(x = "Edema", y = "Survival (%)")+
  scale_color_brewer(palette="Set1")+
  facet_grid(~group)+
  theme(legend.position="none")

## ----interaction-show, echo = TRUE, eval = FALSE----------------------------------------
#  pbc_interaction <- find.interaction(pbc_rf)
#  
#  ggint <- gg_interaction(pbc_interaction)
#  plot(ggint, x_var = "bili") +
#    labs(y = "Interactive Minimal Depth")

## ----interaction, fig.cap = "Minimal depth variable interaction with bilirubin (marked with red cross). Higher values indicate lower interactivity with target variable."----
#pDat.int <- find.interaction(pbc_rf)
data("pbc_interaction", package="ggRandomForests")

plot(gg_interaction(pbc_interaction), x_var = "bili") + 
  labs(y = "Interactive Minimal Depth")

## ----interactionPanel, echo = TRUE, fig.cap = "Minimal depth variable interaction panel with prothrombin time, albumin, urine copper and edema. Higher values indicate lower interactivity with target variable.", fig.width = 7, fig.height = 4----
plot(gg_interaction(pbc_interaction), x_var = xvar[2:5]) + 
  labs(y = "Interactive Minimal Depth") + 
  theme(legend.position = "none")

## ----var_dep, echo = TRUE, fig.cap = "Variable dependence plot. Survival at 1 year against bilirubin. Individual cases are marked with blue circles (alive or censored) and red x (dead). Loess smooth curve indicates the trend as bilirubin  increases."----
ggrf <- gg_variable(pbc_rf, time = 1)

ggvar <- ggrf
ggvar$treatment <- as.numeric(ggvar$treatment)
ggvar$treatment[which(ggvar$treatment==1)] <- "DPCA" 
ggvar$treatment[which(ggvar$treatment==2)] <- "placebo" 
ggvar$treatment <- factor(ggvar$treatment)

ggvar$stage <- paste("stage=", ggvar$stage, sep="")

var_dep <- plot(ggvar, x_var = "bili", smooth = TRUE, 
                method = "loess", span=1.5,alpha = .5, se = FALSE) + 
  labs(y = "Survival", 
       x = dta.labs["bili", "label"]) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)

show(var_dep)

## ----coplot_bilirubin, echo = TRUE, fig.cap = "Variable dependence coplot. Survival at 1 year against bilirubin, stratified by treatment and histological stage.", fig.width = 7, fig.height = 4----
var_dep + 
  facet_grid(treatment~stage)

## ----age-coplot, fig.cap="Variable dependence coplot. Survival at 1 year against bilirubin, stratified by continous variable age.", fig.width = 7, fig.height = 4, echo=TRUE----
age_grp <- cut(pbc_rf$xvar$age, breaks=seq(0,100,10))
ggvar$age_grp <- paste("age=",age_grp, sep="")

var_dep <- plot(ggvar, x_var = "bili", smooth = TRUE, 
                method = "loess", span=1.5,alpha = .5, se = FALSE) + 
  labs(y = "Survival", x = dta.labs["bili", "label"]) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)+ 
  facet_wrap(~age_grp)

var_dep

## ----bilirubin-age-partial, echo=TRUE, eval=FALSE---------------------------------------
#  
#  # Get the training data to work with...
#  dta.train  <- pbc_rf$xvar
#  dta.train$age_grp <- age_grp
#  
#  # Create a series of coplot subsets....
#  lng <- length(levels(age_grp))
#  sbst <- mclapply(1:lng, function(ind){
#    st <- which(dta.train$age_grp==levels(age_grp)[ind])
#    if(length(st) == 0) NULL
#    else st
#    })
#  
#  lvl <- levels(age_grp)
#  # Collapse the subset list to interesting items
#  # (those with observations)
#  # If you work backwards, you do extra tests, but it
#  # cuts the correct items. Cute.
#  for(ind in lng:1){
#    if(is.null(sbst[[ind]])){
#      sbst[[ind]] <- NULL
#  
#      # reset the levels, so we can label things later
#      lvl <- lvl[-ind]
#      }
#  }
#  
#  pDat.partlist <- lapply(1:length(sbst), function(ind){
#    plot.variable(pbc_rf, surv.type="surv", time=1,
#                             subset = sbst[[ind]],
#                          xvar.names="bili", partial=TRUE,
#                          show.plots = FALSE)
#    })
#  
#  gg_part <- mclapply(pDat.partlist, gg_partial)
#  
#  # Flip y-axis
#  cls <- class(gg_part)
#  class(gg_part) <-  c("gg_partial_list", cls)
#  
#  for(ind in 1:length(gg_part)){
#  gg_part[[ind]]$age <- lvl[ind]
#  }
#  pbc_prtl_bili_age <- do.call(rbind, gg_part)
#  pbc_prtl_bili_age$age <- paste("Age=", gg_merge$age)
#  pbc_prtl_bili_age$age <- factor(gg_merge$age)
#  
#  ggpl <- ggplot(pbc_prtl_bili_age,
#                 aes(x=bili, y=yhat, shape=age, color=age))+
#    geom_point()+geom_smooth(se=FALSE)+
#    labs(x="Surgical Date", y="Survival 1 year")+
#    scale_color_brewer(palette="Set1")
#  ggpl

## ----bili-age, fig.cap="Partial (risk adjusted) variable dependence coplot. Survival at 1 year against bilirubin, stratified by age groups. Points mark risk adjusted estimates, loess smooth indicates predicted trend within each age group as a function of bilirubin.", fig.width = 7, fig.height = 4, echo=FALSE----
data(pbc_prtl_bili_age, package="ggRandomForests")
ggpl <- ggplot(pbc_prtl_bili_age, aes(x=bili, y=yhat, shape=age, color=age))+
  geom_point()+geom_smooth(se=FALSE)+
  labs(x="Surgical Date", y="Survival 1 year")+
  scale_color_brewer(palette="Set1")
ggpl

## ----albumin-coplot, fig.cap="Variable dependence coplot. Survival at 1 year against bilirubin, stratified by continous variable albumin.", fig.width = 7, fig.height = 4, echo=TRUE----
albumin_grp <- cut(pbc_rf$xvar$albumin, breaks=c(0,seq(3,3.5,.5),5))
ggvar$albumin_grp <- paste("albumin=",albumin_grp, sep="")

var_dep <- plot(ggvar, x_var = "bili", smooth = TRUE, 
                method = "loess", span=1.5,alpha = .5, se = FALSE) + 
  labs(y = "Survival", x = dta.labs["bili", "label"]) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)+ 
  facet_wrap(~albumin_grp)

var_dep

## ----bilirubin-albumin-partial, echo=TRUE, eval=FALSE-----------------------------------
#  
#  # Get the training data to work with...
#  dta.train  <- pbc_rf$xvar
#  dta.train$albumin_grp <- ggvar$albumin_grp
#  
#  # Create a series of coplot subsets....
#  lng <- length(levels(albumin_grp))
#  sbst <- mclapply(1:lng, function(ind){
#    st <- which(dta.train$albumin_grp==levels(albumin_grp)[ind])
#    if(length(st) == 0) NULL
#    else st
#    })
#  
#  lvl <- levels(albumin_grp)
#  # Collapse the subset list to interesting items
#  # (those with observations)
#  # If you work backwards, you do extra tests, but it
#  # cuts the correct items. Cute.
#  for(ind in lng:1){
#    if(is.null(sbst[[ind]])){
#      sbst[[ind]] <- NULL
#  
#      # reset the levels, so we can label things later
#      lvl <- lvl[-ind]
#      }
#  }
#  
#  pDat.partlist <- mclapply(1:length(sbst), function(ind){
#    plot.variable(pbc_rf, surv.type="surv", time=1,
#                             subset = sbst[[ind]],
#                          xvar.names="bili", partial=TRUE,
#                          show.plots = FALSE)
#    })
#  
#  gg_part <- mclapply(pDat.partlist, gg_partial)
#  
#  # Flip y-axis
#  cls <- class(gg_part)
#  class(gg_part) <-  c("gg_partial_list", cls)
#  
#  for(ind in 1:length(gg_part)){
#  gg_part[[ind]]$albumin <- lvl[ind]
#  }
#  gg_merge <- do.call(rbind, gg_part)
#  gg_merge$albumin <- paste("albumin=", gg_merge$albumin)
#  gg_merge$albumin <- factor(gg_merge$albumin)
#  
#  ggpl <- ggplot(gg_merge, aes(x=bili, y=yhat, shape=albumin, color=albumin))+
#    geom_point()+geom_smooth(se=FALSE)+
#    labs(x="Surgical Date", y="Survival 1 year")+
#    scale_color_brewer(palette="Set1")
#  ggpl

## ----bili-albumin, fig.cap="Partial (risk adjusted) variable dependence coplot. Survival at 1 year against bilirubin, stratified by albumin groups. Points mark risk adjusted estimates, loess smooth indicates predicted trend within each age group as a function of bilirubin.", fig.width = 7, fig.height = 4, echo=FALSE----
data(pbc_prtl_bili_albumin, package="ggRandomForests")
ggpl <- ggplot(pbc_prtl_bili_albumin, aes(x=bili, y=yhat, 
                                          shape=albumin, color=albumin))+
  geom_point()+geom_smooth(se=FALSE)+
  labs(x="Surgical Date", y="Survival 1 year")+
  scale_color_brewer(palette="Set1")
ggpl

## ----prothrombin-coplot, fig.cap="Variable dependence coplot. Survival at 1 year against bilirubin, stratified by continous variable prothrombin time.", fig.width = 7, fig.height = 4, echo=TRUE----
ggvar$prothrombin_grp <- cut(pbc_rf$xvar$prothrombin, breaks=c(8.9,10,11,12,18))

var_dep <- plot(ggvar[which(!is.na(ggvar$prothrombin_grp)),],
                x_var = "bili", smooth = TRUE, 
                method = "loess", span=1.5,alpha = .5, se = FALSE) + 
  labs(y = "Survival", x = dta.labs["bili", "label"]) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = strCol, labels = event.labels) + 
  scale_shape_manual(values = event.marks, labels = event.labels)+ 
  facet_wrap(~prothrombin_grp)

var_dep

## ----bilirubin-prothrombin-partial, echo=TRUE, eval=FALSE-------------------------------
#  
#  # Get the training data to work with...
#  dta.train  <- pbc_rf$xvar
#  dta.train$prothrombin_grp <- ggvar$prothrombin_grp
#  
#  # Create a series of coplot subsets....
#  lng <- length(levels(prothrombin_grp))
#  sbst <- mclapply(1:lng, function(ind){
#    st <- which(dta.train$prothrombin_grp==levels(prothrombin_grp)[ind])
#    if(length(st) == 0) NULL
#    else st
#    })
#  
#  lvl <- levels(prothrombin_grp)
#  # Collapse the subset list to interesting items
#  # (those with observations)
#  # If you work backwards, you do extra tests, but it
#  # cuts the correct items. Cute.
#  for(ind in lng:1){
#    if(is.null(sbst[[ind]])){
#      sbst[[ind]] <- NULL
#  
#      # reset the levels, so we can label things later
#      lvl <- lvl[-ind]
#      }
#  }
#  
#  pDat.partlist <- mclapply(1:length(sbst), function(ind){
#    plot.variable(pbc_rf, surv.type="surv", time=1,
#                             subset = sbst[[ind]],
#                          xvar.names="bili", partial=TRUE,
#                          show.plots = FALSE)
#    })
#  
#  gg_part <- mclapply(pDat.partlist, gg_partial)
#  
#  # Flip y-axis
#  cls <- class(gg_part)
#  class(gg_part) <-  c("gg_partial_list", cls)
#  
#  for(ind in 1:length(gg_part)){
#  gg_part[[ind]]$prothrombin <- lvl[ind]
#  }
#  gg_merge <- do.call(rbind, gg_part)
#  gg_merge$prothrombin <- paste("prothrombin=", gg_merge$prothrombin)
#  gg_merge$prothrombin <- factor(gg_merge$prothrombin)
#  
#  ggpl <- ggplot(gg_merge, aes(x=bili, y=yhat, shape=prothrombin, color=prothrombin))+
#    geom_point()+geom_smooth(se=FALSE)+
#    labs(x="Surgical Date", y="Survival 1 year")+
#    scale_color_brewer(palette="Set1")
#  ggpl

## ----bili-prothrombin, fig.cap="Partial (risk adjusted) variable dependence coplot. Survival at 1 year against bilirubin, stratified by prothrombin groups. Points mark risk adjusted estimates, loess smooth indicates predicted trend within each age group as a function of bilirubin.", fig.width = 7, fig.height = 4, echo=FALSE----
data(pbc_prtl_bili_prothrombin, package="ggRandomForests")
ggpl <- ggplot(pbc_prtl_bili_prothrombin, aes(x=bili, y=yhat, 
                                          shape=prothrombin, color=prothrombin))+
  geom_point()+geom_smooth(se=FALSE)+
  labs(x="Surgical Date", y="Survival 1 year")+
  scale_color_brewer(palette="Set1")
ggpl

