% File sme/man/sme.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{sme}
\title{Smoothing-splines mixed-effects models}
\usage{
sme(object,tme,ind,verbose=F,lambda.mu=NULL,lambda.v=NULL,maxIter=500,
knots=NULL,zeroIntercept=F,deltaEM=1e-3,deltaNM=1e-3,criteria="AICc",
initial.lambda.mu=10000,initial.lambda.v=10000,normalizeTime=FALSE,\dots)
}
\alias{sme}
\alias{sme.default}
\arguments{
  \item{object}{either a vector of observations, a \code{data.frame} object or a \code{list} of
  vectors of observations. The method functions \code{sme.data.frame} and \code{sme.list} are
  documented separately}
  \item{tme}{either a vector of time points corresponding to the observations given in \code{object}
  or a list of vectors of time points in the case of \code{sme.list}. Ignored in the case of
  \code{sme.data.frame}}
  \item{ind}{a factor (or a vector that can be coerced to a factor) of subject identifiers
  corresponding to the observations given in \code{object} or a list of factors of subject
  identifiers in the case of \code{sme.list}. Ignored   in the case of \code{sme.data.frame}}
  \item{verbose}{if \code{TRUE}, debug information will be output while fitting the model}
  \item{lambda.mu}{smoothing parameter used for the fixed-effect function. If \code{NULL}, the
  optimal values for this and \code{lambda.v} will be found according to \code{criteria} using
  Nelder-Mead search}
  \item{lambda.v}{smoothing parameter used for the random-effects functions. If \code{NULL}, the
  optimal values for this and \code{lambda.mu} will be found according to \code{criteria} using
  Nelder-Mead search}
  \item{maxIter}{maximum number of iterations to be performed for the EM algorithm}
  \item{knots}{location of spline knots. If \code{NULL}, an incidence matrix representation will be
  used. See `Details'}
  \item{zeroIntercept}{experimental feature. If \code{TRUE}, the fitted values of the fixed- and
  random-effects functions at the intercept will be zero}
  \item{deltaEM}{convergence tolerance for the EM algorithm}
  \item{deltaNM}{(relative) convergence tolerance for the Nelder-Mead optimisation}
  \item{criteria}{one of \code{"AICc"}, \code{"AIC"}, \code{"BICN"} or \code{"BICn"} indicating
  which criteria to use to score a particular combination of \code{lambda.mu} and \code{lambda.v} in
  the Nelder-Mead search}
  \item{initial.lambda.mu}{value to initialise the smoothing parameter for the fixed-effects to in
  the Nelder-Mead search. See details below}
  \item{initial.lambda.v}{value to initialise the smoothing parameter for the random-effects to in
  the Nelder-Mead search. See details below}
  \item{normalizeTime}{should time be normalized to lie in $[0,1]$? See details below}
  \item{\dots}{additional arguments to \code{sme.data.frame} or \code{sme.list}}
}
\description{
  This generic function fits a smoothing-splines mixed-effects model
}
\value{
  An object of class \code{sme} representing the smoothing-splines mixed-effects model fit. See
  \code{smeObject} for the components of the fit and \code{plot.sme} for visualisation options
}
\details{
Prior to package version 0.9, starting values for the smoothing parameters in the Nelder-Mead search
were fixed to $10000$ for both \code{lambda.mu} and \code{lambda.v}. As it turns out, the
appropriate scale for the smoothing parameters depends on the scale for \code{tme} and so \code{tme}
will now automatically be rescaled to lie in $[0,1]$ and much smaller initial values for the
smoothing parameters will be used, although these can now optionally changed to achieve best
results. To reproduce results obtained using previous versions of the package, set
\code{initial.lambda.mu=10000}, \code{initial.lambda.v=10000} and \code{normalizeTime=FALSE}.

The default behaviour is to use an incidence matrix representation for the smoothing-splines. This
works well in most situations but may incur a high computational cost when the number of distinct
time points is large, as may be the case for irregularly sampled data. Alternatively, a basis
projection can be used by giving a vector of \code{knots} of length (much) less than the number of
distinct time points.
}
\references{
  Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
}
\author{Maurice Berk \email{maurice.berk01@imperial.ac.uk}}
\seealso{\code{\link{smeObject}}, \code{\link{sme.data.frame}}, \code{\link{sme.list}}, \code{\link{plot.sme}}}
\examples{
  data(MTB)
  fit.AIC <- sme(MTB[MTB$variable==6031,c("y","tme","ind")],criteria="AIC")
  fit.BICN <- sme(MTB[MTB$variable==6031,c("y","tme","ind")],criteria="BICN")
  fit.BICn <- sme(MTB[MTB$variable==6031,c("y","tme","ind")],criteria="BICn")
  fit.AICc <- sme(MTB[MTB$variable==6031,c("y","tme","ind")],criteria="AICc")

  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")],lambda.mu=1e5,lambda.v=1e5)

  data(inflammatory)
  system.time(fit <- sme(inflammatory,knots=c(29.5,57,84.5),deltaEM=0.1,deltaNM=0.1))
}
