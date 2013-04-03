% File sme/man/plot.sme.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later

\name{plot.sme}
\title{Visualisation of smoothing-splines mixed-effects models}
\usage{
\method{plot}{sme}(x,type="model",...)
}
\alias{plot.sme}
\arguments{
  \item{x}{an object inheriting from class \code{sme}, representing a fitted smoothing-splines
  mixed-effects model}
  \item{type}{the type of plot desired, currently either \code{"model"} for visualising the model
  fit, \code{"raw"} for visualising the raw data or \code{"diagnostic"} for visualising the model
  residuals}
  \item{...}{further arguments to \code{plotSmeModel}, \code{plotSmeRaw} or
  \code{plotSmeDiagnostic}}
}
\description{
  A visualisation of the smoothing-splines mixed-effects model fit is obtained. For
  \code{type="model"}, a scatter plot of the observations overlaid with the fitted mean and
  (optionally) individual curves. Confidence limits for the mean curve can optionally also be shown.
  See \code{plot.sme.model} for further details. For \code{type="raw"}, a trellis plot of the raw
  data for each subject, with the fitted individual curves optionally overlaid. See
  \code{plot.sme.raw} for further details. For \code{type="diagnostic"}, a diagnostic plot of the
  model residuals. See \code{plot.sme.diagnostic} for further details.
}
\value{
  Depending on \code{type}, either a plot of the fitted model, a plot of the raw data or a
  diagnostic plot of the model residuals
}
\references{
  Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
}
\author{Maurice Berk \email{maurice@mauriceberk.com}}
\seealso{\code{\link{plotSmeModel}}, \code{\link{plotSmeRaw}},
\code{\link{plotSmeDiagnostic}}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  plot(fit,type="model")
  plot(fit,type="raw")
  plot(fit,type="diagnostic")
}
