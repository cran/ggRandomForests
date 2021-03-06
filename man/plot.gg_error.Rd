% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.gg_error.R
\name{plot.gg_error}
\alias{plot.gg_error}
\title{Plot a \code{\link{gg_error}} object}
\usage{
\method{plot}{gg_error}(x, ...)
}
\arguments{
\item{x}{gg_error object created from a \code{\link[randomForestSRC]{rfsrc}} object}

\item{...}{extra arguments passed to \code{ggplot} functions}
}
\value{
\code{ggplot} object
}
\description{
A plot of the cumulative OOB error rates of the random forest as a 
function of number of trees.
}
\details{
The gg_error plot is used to track the convergence of the 
randomForest. This figure is a reproduction of the error plot
from the \code{\link[randomForestSRC]{plot.rfsrc}} function.
}
\examples{
\dontrun{
## Examples from RFSRC package... 
## ------------------------------------------------------------
## classification example
## ------------------------------------------------------------
## ------------- iris data
## You can build a randomForest
# rfsrc_iris <- rfsrc(Species ~ ., data = iris)
# ... or load a cached randomForestSRC object
data(rfsrc_iris, package="ggRandomForests")

# Get a data.frame containing error rates
gg_dta<- gg_error(rfsrc_iris)

# Plot the gg_error object
plot(gg_dta)

## ------------------------------------------------------------
## Regression example
## ------------------------------------------------------------
## ------------- airq data
# rfsrc_airq <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")
# ... or load a cached randomForestSRC object
data(rfsrc_airq, package="ggRandomForests")

# Get a data.frame containing error rates
gg_dta<- gg_error(rfsrc_airq)

# Plot the gg_error object
plot(gg_dta)

## ------------- Boston data
data(rfsrc_Boston, package="ggRandomForests")

# Get a data.frame containing error rates
gg_dta<- gg_error(rfsrc_Boston)

# Plot the gg_error object
plot(gg_dta)

## ------------- mtcars data
data(rfsrc_mtcars, package="ggRandomForests")

# Get a data.frame containing error rates
gg_dta<- gg_error(rfsrc_mtcars)

# Plot the gg_error object
plot(gg_dta)

## ------------------------------------------------------------
## Survival example
## ------------------------------------------------------------
## ------------- veteran data
## randomized trial of two treatment regimens for lung cancer
# data(veteran, package = "randomForestSRC")
# rfsrc_veteran <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100)

# Load a cached randomForestSRC object
data(rfsrc_veteran, package="ggRandomForests")

gg_dta <- gg_error(rfsrc_veteran)
plot(gg_dta)

## ------------- pbc data
# Load a cached randomForestSRC object
data(rfsrc_pbc, package="ggRandomForests")

gg_dta <- gg_error(rfsrc_pbc)
plot(gg_dta)
}
}
\references{
Breiman L. (2001). Random forests, Machine Learning, 45:5-32.

Ishwaran H. and Kogalur U.B. (2007). Random survival forests for R, Rnews, 
7(2):25-31.

Ishwaran H. and Kogalur U.B. (2013). Random Forests for Survival, Regression 
and Classification (RF-SRC), R package version 1.4.
}
\seealso{
\code{\link{gg_error}} \code{\link[randomForestSRC]{rfsrc}}
 \code{\link[randomForestSRC]{plot.rfsrc}}
}

