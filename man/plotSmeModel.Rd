% File sme/man/plotSmeModel.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later

\name{plotSmeModel}
\title{Visualisation of a fitted smoothing-splines mixed-effects model}
\usage{
plotSmeModel(x,xlab="Time",ylab="Y",showIndividuals=T,showConfidenceBands=F,
col.meanCurve="red",...)
}
\alias{plotSmeModel}
\arguments{
  \item{x}{an object inheriting from class \code{sme}, representing a fitted smoothing-splines
  mixed-effects model}
  \item{xlab}{the label for the x-axis, defaulting to \code{"Time"}}
  \item{ylab}{the label for the y-axis, defaulting to \code{"Y"}}
  \item{showIndividuals}{if \code{TRUE} then each of the fitted subject specific curves will be
  displayed in the plot as dashed lines, otherwise only the mean curve will be shown}
  \item{showConfidenceBands}{if \code{TRUE}, then the 95\% confidence band for the fitted mean curve
  will  be shown}
  \item{col.meanCurve}{the color to be used for the fitted mean curve and,
  if \code{showConfidenceBands} is \code{TRUE}, the confidence band as well}
  \item{\dots}{additional arguments passed through to \code{plot} when producing the scatter plot
  of observations}
}
\description{
  The observations used to fit the given smoothing-splines mixed-effects model are shown as a
  scatter plot and the fitted mean curve is overlaid. Optionally, the fitted individual curves are
  also shown. Also, the 95\% confidence band for the fitted mean curve can optionally be displayed.
}
\value{
  A plot visualising the fitted smoothing-splines mixed-effects model as described above.
}
\references{
  Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
}
\author{Maurice Berk \email{maurice@mauriceberk.com}}
\seealso{\code{\link{plot.sme}}, \code{\link{plotSmeRaw}}, \code{\link{plotSmeDiagnostic}}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  plotSmeModel(fit,showIndividuals=FALSE)
  plotSmeModel(fit,showConfidenceBands=TRUE)
}