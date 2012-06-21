% File sme/man/smeObject.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{smeObject}
\title{Fitted sme Object}
\alias{smeObject}
\description{
  An object returned by the \code{sme} function, inheriting from class
  \code{sme}, representing a fitted smoothing-splines mixed-effects model}
\value{
  An \code{sme} object must contain the following components.
  \item{call}{a list containing an image of the \code{sme} call that produced the object.}
  \item{coefficients}{an \emph{n+1} by \emph{p} matrix where \emph{n} is the number of subjects and
  \emph{p} is the number of knots used in the spline basis. The first row of the matrix corresponds
  to the fitted values of the mean curve at the knots, and the remaining \emph{n} rows correspond
  to the fitted values of the individual deviations from the mean curve at the knots.}
  \item{fitted}{a vector of model fitted values corresponding to the original observations.}
  \item{logLik}{the log-likelihood of the fitted model.}
  \item{residuals}{a vector of model residuals.} 
  \item{nobs}{the number of observations.}
  \item{data}{a data frame of the original data used to fit the model with variables \code{y},
  \code{tme} and \code{ind} corresponding to observations, time points and subject identifiers
  respectively.}
  \item{df}{a vector consisting of named components \code{mu} and \code{v} corresponding to the
  degrees of freedom of the mean and subject curves respectively.}
  \item{smoothingParameters}{a vector of named components \code{mu} and \code{v} corresponding to
  the smoothing parameters for the mean and subject curves respectively.}
  \item{parameters}{a list with named components \code{sigmaSquared} and \code{D} corresponding to
  the error variance and (unregularized) variance of the random-effects respectively.}
  \item{iterations}{the number of iterations of the EM algorithm that ran before convergence (or the
  limit was reached)}
  \item{info}{a numeric code indicating diagnostic information from the EM algorithm with zero
  indicating a successful run.}
  In some instances, an \code{sme} object may also contain the following components.
  \item{knots}{a vector of the knots used in the spline basis when they do not coincide with the
  unique design time points.}
}

\author{Maurice Berk \email{maurice.berk01@imperial.ac.uk}}

\seealso{\code{\link{sme}}}
