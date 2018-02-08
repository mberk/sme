% File sme/man/sme.data.frame.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later
\name{sme.list}
\title{Carry out mulitple independent smoothing-splines mixed-effects model fits simultaneously}
\usage{
\method{sme}{list}(object,tme,ind,verbose=F,lambda.mu=NULL,lambda.v=NULL,maxIter=500,
knots=NULL,zeroIntercept=F,deltaEM=1e-3,deltaNM=1e-3,criteria="AICc",
initial.lambda.mu=10000,initial.lambda.v=10000,normalizeTime=FALSE,numberOfThreads=-1,
\dots)
}
\alias{sme.list}
\arguments{
  \item{object}{a list of vectors of observations}
  \item{tme}{a list of vectors of time points corresponding to the observations in \code{object}}
  \item{ind}{a list of factors (or vectors that can be coerced to factors) of subject identifiers
  corresponding to the observations in \code{object}}
  \item{verbose}{if \code{TRUE}, debug information will be output while fitting the model(s)}
  \item{lambda.mu}{either a single smoothing parameter to be used for the fixed-effect function for
  all fits, or a vector of smoothing parameters, one for each fit, or \code{NULL} if Nelder-Mead
  search should be used to find the optimal values for this and \code{lambda.v} for all fits}
  \item{lambda.v}{either a single smoothing parameter to be used for the random-effects functions
  for all fits, or a vector of smoothing parameters, one for each fit, or \code{NULL} if Nelder-Mead
  search should be used to find the optimal values for this and \code{lambda.v} for all fits}
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
  \item{numberOfThreads}{The number of threads to use to fit the multiple smoothing-splines
  mixed-effects models simultaneously. When \code{numberOfThreads=-1}, as is the default, the
  OpenMP system will handle thread creation dynamically}
  \item{\dots}{additional arguments, currently not used}
}
\description{
  Carry out multiple independent smoothing-splines mixed-effects model fits simultaneously
}
\value{
  A list of objects of class \code{sme}. See \code{smeObject} for the components of the fit and
  \code{plot.sme} for visualisation options
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
\author{Maurice Berk \email{maurice@mauriceberk.com}}
\seealso{\code{\link{smeObject}}, \code{\link{sme}}, \code{\link{sme.data.frame}},
\code{\link{plot.sme}}}
