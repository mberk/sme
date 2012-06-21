% File sme/man/BICn.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{BICn}
\alias{BICn}
\title{An alternative BIC for longitudinal models}
\description{Calculates the BIC as \emph{-2 * log-likelihood + log(n) * npar} for a longitudinal
model where \emph{npar} is the number of parameters in the fitted-model and \emph{n} is the
\emph{number of subjects}}

\usage{BICn(object,\dots)}

\arguments{
  \item{object}{a fitted longitudinal model object}
  \item{\dots}{some methods for this generic function may require additional arguments}
}

\details{
When applying the BIC in a longitudinal context, there is some debate as to whether the sample size
should be taken to mean the number of subjects or the total number of observations across all
subjects (see Section 7.3 of Hedeker and Gibbons, 2006).

Assuming the default \code{BIC} function accounts for the latter case, this generic function can be
implemented for longitudinal models where the number of subjects can be extracted in order to
calculate the BIC under the alternative definition.
}

\value{
A numeric value with the BIC of the longitudinal model, with the penalty taken as a function of the
number of subjects as described.
}

\references{
Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint

Hedeker, D. & Gibbons, D. R. (2006). \emph{Longitudinal Data Analysis}. Wiley
}
\author{Maurice Berk \email{maurice.berk01@imperial.ac.uk}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  BICn(fit)
}