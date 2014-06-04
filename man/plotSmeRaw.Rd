% File sme/man/plotSmeRaw.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later

\name{plotSmeRaw}
\title{Visualisation of raw data for a smoothing-splines mixed-effects model}
\usage{
plotSmeRaw(x,xlab="Time",ylab="Y",mainTitle="",showModelFits=TRUE,showRawLines=FALSE,...)
}
\alias{plotSmeRaw}
\arguments{
  \item{x}{an object inheriting from class \code{sme}, representing a fitted smoothing-splines
  mixed-effects model}
  \item{xlab}{the label for the x-axis, defaulting to \code{"Time"}}
  \item{ylab}{the label for the y-axis, defaulting to \code{"Y"}}
  \item{mainTitle}{the main title for the plot, defaulting to none}
  \item{showModelFits}{if \code{TRUE}, the fitted individual curve for each subject will be overlaid
  on   their observations in their panel}
  \item{showRawLines}{if \code{TRUE}, the raw observations will be connected by lines}
  \item{\dots}{additional arguments passed through to \code{xyplot}}
}
\description{
  A visualisation of the raw data used to fit the smoothing-splines mixed-effects model is obtained.
  This is presented as a trellis plot where each panel shows the observations for one of the
  subjects as a scatter plot (with optional lines connecting the observations if \code{showRawLines}
  is \code{TRUE}). If \code{showModelFits} is \code{TRUE}, then these observations are overlaid
  with the fitted individual curve for that subject.
}
\value{
  A trellis plot as described above.
}
\references{
  Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
}
\author{Maurice Berk \email{maurice@mauriceberk.com}}
\seealso{\code{\link{plot.sme}}, \code{\link{plotSmeModel}}, \code{\link{plotSmeDiagnostic}}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  plotSmeRaw(fit,showModelFits=FALSE)
  plotSmeRaw(fit,showModelFits=FALSE,showRawLines=TRUE)
  plotSmeRaw(fit)
}