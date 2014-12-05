####**********************************************************************
####**********************************************************************
####
####  ----------------------------------------------------------------
####  Written by:
####  ----------------------------------------------------------------
####    John Ehrlinger, Ph.D.
####    Assistant Staff
####    Dept of Quantitative Health Sciences
####    Learner Research Institute
####    Cleveland Clinic Foundation
####
####    email:  john.ehrlinger@gmail.com
####    URL:    https://github.com/ehrlinger/ggRandomForests
####  ----------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************
#'
#' Plot a \code{\link{gg_minimal_depth}} object for random forest variable ranking.
#' 
#' @param x \code{\link{gg_minimal_depth}} object created from a \code{randomForestSRC::rfsrc} object
#' @param selection should we restrict the plot to only include variables selected by the 
#' minimal depth criteria (boolean).
#' @param list_vars add text list of ranked variables. Only used if selection = TRUE (boolean)
#' @param type select type of y axis labels c("named","rank")
#' @param ... optional arguments passed to \code{\link{gg_minimal_depth}}
#' 
#' @return \code{ggplot} object
#' 
#' @seealso \code{randomForestSRC::var.select} \code{\link{gg_minimal_depth}}
#' @export plot.gg_minimal_depth
#' 
#' @references
#' Breiman L. (2001). Random forests, Machine Learning, 45:5-32.
#' 
#' Ishwaran H. and Kogalur U.B. (2007). Random survival forests for R, 
#' Rnews, 7(2):25-31.
#' 
#' Ishwaran H. and Kogalur U.B. (2014). Random Forests for Survival, 
#' Regression and Classification (RF-SRC), R package version 1.5.
#' 
#' @examples
#' \dontrun{
#' #' ## Examples from RFSRC package... 
#' ## ------------------------------------------------------------
#' ## classification example
#' ## ------------------------------------------------------------
#' ## You can build a randomForest
#' # iris_rf <- rfsrc(Species ~ ., data = iris)
#' # iris_vs <- var.select(iris_rf)
#' # ... or load a cached randomForestSRC object
#' data(iris_vs, package="ggRandomForests")
#' 
#' # Get a data.frame containing minimaldepth measures
#' ggrf.obj<- gg_minimal_depth(iris_vs)
#' 
#' # Plot the gg_mkinimal_depth object
#' plot(ggrf.obj)
#' 
#' ## ------------------------------------------------------------
#' ## Regression example
#' ## ------------------------------------------------------------
#' # airq_rf <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")
#' # airq_vs <- var.select(airq_rf)
#' # ... or load a cached randomForestSRC object
#' data(airq_vs, package="ggRandomForests")
#' 
#' # Get a data.frame containing error rates
#' ggrf.obj<- gg_minimal_depth(airq_vs)
#' 
#' # Plot the gg_error object
#' plot(ggrf.obj)
#' 
#' ## ------------------------------------------------------------
#' ## Survival example
#' ## ------------------------------------------------------------
#' ## veteran data
#' ## randomized trial of two treatment regimens for lung cancer
#' # data(veteran, package = "randomForestSRC")
#' # veteran_rf <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100)
#' # veteran_vs <- var.select(veteran_rf)
#' # Load a cached randomForestSRC object
#' data(veteran_vs, package="ggRandomForests")
#' 
#' ggrf.obj <- gg_minimal_depth(veteran_vs)
#' plot(ggrf.obj)
#' 
#' }
#' @importFrom ggplot2 ggplot geom_line theme aes_string labs coord_cartesian geom_text annotate geom_hline coord_flip geom_vline
### error rate plot
plot.gg_minimal_depth <- function(x, selection=FALSE, 
                                  list_vars = TRUE, 
                                  type=c("named","rank"),
                                  ...){
  object <- x
  if(!inherits(x, "gg_minimal_depth")){
    object <- gg_minimal_depth(x, ...)
  }
  type=match.arg(type)
  
  xl <-c(0,ceiling(max(object$varselect$depth))+1)
  sel.th = object$md.obj$threshold
  
  if(selection){
    modelSize<-object$modelsize
    
    # Labels for the top md vars.
    md.labs <- object$topvars
    
    ## Number the variables
    for(ind in 1:length(md.labs)){
      md.labs[ind] <- paste(ind, md.labs[ind], sep=". ")
    }
    vSel <- object$varselect[1:modelSize,]
    vSel$rank <- 1:nrow(vSel)
    
    ## Reorder the minimal depth to place most "important" at top of figure
    vSel$names <- factor(vSel$names, 
                         levels=rev(levels(vSel$names )))
    gDta <- ggplot(vSel)
    gDta <- switch(type,
                   rank = gDta +
                     geom_point(aes_string(y="rank", x="depth", label="rank"))+
                     coord_cartesian(xlim=xl) + 
                     geom_text(aes_string(y="rank", x="depth"-.7, label="rank"), 
                               size=3, hjust=0),
                   named = gDta +
                     geom_point(aes_string(y="depth", x="names"))+
                     coord_cartesian(ylim=xl)
    )
    
    if(list_vars){
      
      # We will need to modify this as we get more examples.
      x.pt <- ceiling(object$md.obj$threshold)
      
      # Populate the text box with min.depth ordered labels
      yPlace= object$modelsize -1
      
      for(ind in 1:length(md.labs)){
        gDta <- gDta +
          annotate(x=x.pt,y=yPlace,geom="text",label=md.labs[ind], hjust=0,size=3)
        yPlace <- yPlace - .95
      }
    }
    
  }else{ 
    vSel <- object$varselect
    vSel$rank <- 1:dim(vSel)[1]
    vSel$names <- factor(vSel$names, 
                         levels=rev(levels(vSel$names )))
    gDta <- ggplot(vSel)
    gDta <- switch(type,
                   rank = gDta +
                     geom_point(aes_string(y="rank", x="depth"))+
                     coord_cartesian(xlim=xl),
                   named  =gDta +
                     geom_point(aes_string(y="depth", x="names"))+
                     coord_cartesian(ylim=xl)
    )
    
  }
  
  if(type=="named"){
    gDta <- gDta+
      geom_hline(yintercept=sel.th, lty=2)+
      labs(y="Minimal Depth of a Variable", x="")+
      coord_flip() 
  }else{
    gDta <- gDta+
      labs(y="Rank", x="Minimal Depth of a Variable")+
      geom_vline(xintercept=sel.th, lty=2)
  }
  return(gDta)
}