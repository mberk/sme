% File sme/man/plotSmeDiagnostic.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later

\name{plotSmeDiagnostic}
\title{Diagnostic plot of the residuals for a fitted smoothing-splines mixed-effects model}
\usage{
plotSmeDiagnostic(x)
}
\alias{plotSmeDiagnostic}
\arguments{
  \item{x}{an object inheriting from class \code{sme}, representing a fitted smoothing-splines
  mixed-effects model}
}
\description{
  This type of plot is heavily inspired by Wu and Zhang (2006). A plot containing four panels is
  produced: (1) standardised residuals against fitted values, (2) standardised residuals against
  time, (3) standardised residuals against response and (4) a normal Q-Q plot of the standardised
  residuals. These plots can be used to assess how well the fitted model has accounted for the
  temporal correlation and between-subject variance in the data.
}
\value{
  A four panel plot as described above.
}
\references{
  Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
  
  Wu, H. and Zhang, J.-T. (2006). \emph{Nonparametric Regression Methods for Longitudinal Data Analysis}. Wiley
}
\author{Maurice Berk \email{maurice@mauriceberk.com}}
\seealso{\code{\link{plot.sme}}, \code{\link{plotSmeRaw}}, \code{\link{plotSmeModel}}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  plotSmeDiagnostic(fit)
}