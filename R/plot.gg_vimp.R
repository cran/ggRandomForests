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
#' Plot a \code{\link{gg_vimp}} object, extracted variable importance of a 
#' \code{randomForestSRC::rfsrc} object
#' 
#' @param x \code{\link{gg_vimp}} object created from a \code{randomForestSRC::rfsrc} object
#' @param n_var restrict the plot to only nvar variable importance measures
#' @param barcolor Vector of booleans to color vimp bars. By default the "positive" column, TRUE-> "blue". 
#' @param label A vector of alternative variable names.
#' @param relative should we plot vimp or relative vimp. Defaults to vimp.
#' @param ... optional arguments passed to gg_vimp if necessary
#' 
#' @return \code{ggplot} object
#' 
#' @export plot.gg_vimp
#' 
#' @seealso \code{\link{gg_vimp}}
#' 
#' @references
#' Breiman L. (2001). Random forests, Machine Learning, 45:5-32.
#' 
#' Ishwaran H. and Kogalur U.B. (2007). Random survival forests for 
#' R, Rnews, 7(2):25-31.
#' 
#' Ishwaran H. and Kogalur U.B. (2013). Random Forests for Survival, 
#' Regression and Classification (RF-SRC), R package version 1.4.
#' 
#' @examples
#' \dontrun{
#' #' ## ------------------------------------------------------------
#' ## classification example
#' ## ------------------------------------------------------------
#' # iris_rf <- rfsrc(Species ~ ., data = iris)
#' data(iris_rf, package="ggRandomForests")
#' ggrf <- gg_vimp(iris_rf)
#' plot(ggrf)
#'  
#' ## ------------------------------------------------------------
#' ## regression example
#' ## ------------------------------------------------------------
#' 
#' # airq.obj <- rfsrc(Ozone ~ ., airquality)
#' data(airq_rf, package="ggRandomForests")
#' ggrf <- gg_vimp(airq_rf)
#' plot(ggrf)
#' 
#' ## ------------------------------------------------------------
#' ## survival example
#' ## ------------------------------------------------------------
#' data(veteran_rf, package="ggRandomForests")
#' ggrf <- gg_vimp(veteran_rf)
#' plot(ggrf)
#'}
#'
#' @importFrom ggplot2 ggplot geom_bar aes_string labs coord_flip facet_grid
### error rate plot
plot.gg_vimp<- function(x, n_var, barcolor, label, relative, ...){
  object  <- x
  if(!inherits(object, "gg_vimp")) object<- gg_vimp(object, ...)
  if(missing(n_var)) n_var <- dim(object)[1]
  if(n_var > dim(object)[1]) n_var <- dim(object)[1]
  
  if(!missing(barcolor)){
    # We have an alternative coloring
    if(length(barcolor)==n_var){
      object$positive[1:n_var]  <- barcolor
    }  
  }
  if(!missing(label)){
    # We have an alternative coloring
    if(length(label)==n_var){
      object$vars <- as.character(object$vars)
      object$vars[1:n_var]  <- label
      object$vars  <- factor(object$vars, levels=object$vars[order(object$vimp)])
    }  
  }
  
  vimp.plt<-ggplot(object[1:n_var,])
  
  if(missing(relative) | is.null(object$rel_vimp)){
    if(length(unique(object$positive))>1){
      vimp.plt<-vimp.plt+
        geom_bar(aes_string(y="vimp", x="vars", fill="positive"), 
                 stat="identity", width=.5, color="black")
    }else{
      vimp.plt<-vimp.plt+
        geom_bar(aes_string(y="vimp", x="vars"), 
                 stat="identity", width=.5, color="black")
    }
    vimp.plt<-vimp.plt+labs(x="", y="Variable Importance")
      
  }else{
    if(length(unique(object$positive))>1){
      vimp.plt<-vimp.plt+
        geom_bar(aes_string(y="rel_vimp", x="vars", fill="positive"), 
                 stat="identity", width=.5, color="black")
    }else{
      vimp.plt<-vimp.plt+
        geom_bar(aes_string(y="rel_vimp", x="vars"), 
                 stat="identity", width=.5, color="black")
    }   
    vimp.plt<-vimp.plt+ 
      labs(x="", y="Relative Variable Importance") 
  }
  
  if(is.null(object$set))
    vimp.plt<-vimp.plt+ 
    coord_flip()
  else  
    vimp.plt<-vimp.plt+ 
    coord_flip()+facet_grid(~set)
  
  return(vimp.plt)
}