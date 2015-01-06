####**********************************************************************
####**********************************************************************
####  ----------------------------------------------------------------
####  Written by:
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
#' Minimal Depth Variable Interaction data object (\code{randomForestSRC::find.interaction}). 
#' 
#' Converts the matrix returned from
#' \code{randomForestSRC::find.interaction} to a \code{data.frame} and add attributes 
#' for S3 identification. 
#' If passed  a \code{randomForestSRC::rfsrc} object, \code{gg_interaction} first runs 
#' the \code{randomForestSRC::find.interaction} 
#' function with all optional arguments.
#'
#' @param object a \code{randomForestSRC::rfsrc} object or the output from the
#' \code{randomForestSRC::find.interaction} function call.
#' @param ... optional extra arguments passed to \code{randomForestSRC::find.interaction}.
#' 
#' @return \code{gg_interaction} object
#' 
#' @seealso \code{randomForestSRC::rfsrc} 
#' \code{randomForestSRC::find.interaction} 
#' \code{randomForestSRC::max.subtree} 
#' \code{randomForestSRC::var.select} 
#' \code{randomForestSRC::vimp}
#' \code{\link{plot.gg_interaction}} 
#' 
#' @export gg_interaction gg_interaction.rfsrc 
#' @aliases gg_interaction
#' 
#' @importFrom randomForestSRC find.interaction
#' 
#' @references
#' Ishwaran H. (2007). Variable importance in binary regression trees and 
#' forests, Electronic J. Statist., 1:519-537.
#' 
#' Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and Lauer M.S. (2010).
#' High-dimensional variable selection for survival data. J. Amer. Statist. 
#' Assoc., 105:205-217.
#' 
#' Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2011). Random survival 
#' forests for high-dimensional data. Statist. Anal. Data Mining, 4:115-132.
#' 
#' @examples
#' ## Examples from randomForestSRC package... 
#' ## ------------------------------------------------------------
#' ## find interactions, classification setting
#' ## ------------------------------------------------------------
#' ## -------- iris data
#' ## iris.obj <- rfsrc(Species ~., data = iris)
#' ## TODO: VIMP interactions not handled yet....
#' ## find.interaction(iris.obj, method = "vimp", nrep = 3)
#' ## interaction_iris <- find.interaction(iris.obj)
#' data(interaction_iris, package="ggRandomForests")
#' gg_dta <- gg_interaction(interaction_iris)
#' 
#' plot(gg_dta, xvar="Petal.Width")
#' plot(gg_dta, panel=TRUE)
#' 
#' ## ------------------------------------------------------------
#' ## find interactions, regression setting
#' ## ------------------------------------------------------------
#' ## -------- air quality data
#' ## airq.obj <- rfsrc(Ozone ~ ., data = airquality)
#' ##
#' ## TODO: VIMP interactions not handled yet....
#' ## find.interaction(airq.obj, method = "vimp", nrep = 3)
#' ## interaction_airq <- find.interaction(airq.obj)
#' data(interaction_airq, package="ggRandomForests")
#' gg_dta <- gg_interaction(interaction_airq)
#' 
#' \dontrun{
#' plot(gg_dta, xvar="Temp")
#' plot(gg_dta, xvar="Solar.R")
#' }
#' plot(gg_dta, panel=TRUE)
#' 
#' ## -------- Boston data
#' data(interaction_Boston, package="ggRandomForests")
#' gg_dta <- gg_interaction(interaction_Boston)
#' 
#' plot(gg_dta, panel=TRUE)
#' 
#' \dontrun{
#' ## -------- mtcars data
#' data(interaction_mtcars, package="ggRandomForests")
#' gg_dta <- gg_interaction(interaction_mtcars)
#' 
#' plot(gg_dta, panel=TRUE)
#' }
#' 
#' ## ------------------------------------------------------------
#' ## find interactions, survival setting
#' ## ------------------------------------------------------------
#' ## -------- pbc data
#' ## data(pbc, package = "randomForestSRC") 
#' ## pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, nsplit = 10)
#' ## interaction_pbc <- find.interaction(pbc.obj, nvar = 8)
#' data(interaction_pbc, package="ggRandomForests")
#' gg_dta <- gg_interaction(interaction_pbc)
#' 
#' plot(gg_dta, xvar="bili")
#' plot(gg_dta, panel=TRUE)
#' 
#' ## -------- veteran data
#' data(interaction_veteran, package="ggRandomForests")
#' gg_dta <- gg_interaction(interaction_veteran)
#' 
#' plot(gg_dta, panel=TRUE)
#' 
#' 

# gg_interaction <- function (object, ...) {
#   UseMethod("gg_interaction", object)
# }

gg_interaction.rfsrc <- function(object, ...){
  if (inherits(object, "rfsrc")) {
    # If we called this with a rfsrc object, we need to run find.interaction.
    warning("Forest object means we assume max.subtree method for finding interactions.\nThis may take some time.")
    
    object <- find.interaction(object,...)
  }
  if(!inherits(object, "matrix")){
    stop("gg_interaction expects a rfsrc or find.interaction object.")
  }
  
  # make the matrix a data.frame
  gg_dta <- data.frame(object)
  
  # Check to make sure it's the right type of structre...
  if(nrow(gg_dta)!= ncol(gg_dta)){
    stop("gg_interaction expects a find.interaction object.")
  }
  class(gg_dta) <- c("gg_interaction",  class(gg_dta))
  
  invisible(gg_dta)
}
gg_interaction <- gg_interaction.rfsrc