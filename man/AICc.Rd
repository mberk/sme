% File sme/man/AICc.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{AICc}
\alias{AICc}
\title{The corrected AIC}
\description{Calculates the \emph{corrected} AIC (AICc) of Hurvich and Tsai (1989). The AICc
modifies the standard AIC with a correction for small sample sizes.}

\usage{AICc(object)}

\arguments{
  \item{object}{a fitted model object for which there exists a \code{logLik} method to extract the
  corresponding log-likelihood or an object inheriting from class \code{logLik}}
}

\details{
AIC is an asymptotic result and may be inappropriate when the sample size is small, the number of
estimated model parameters is large, or both. For the specific case of a linear model with
homogeneous errors then Hurvich and Tsai (1989) derived a corrected AIC (AICc) which includes a
correction for small sample sizes as \emph{AICc = AIC + (2*k*(k+1))/(n-k-1)} where \emph{AIC} is the
standard AIC, \emph{k} is the number of parameters in the model and \emph{n} is the number of
observations.

While this is an exact result, it only applies in the very specific circumstances in which it was
derived. However, as Burnham and Anderson (2002) point out, whenever the sample size is small some
form of correction to the standard AIC is necessary, to the extent that they argue the AICc of
Hurvich and Tsai (1989) should be used regardless of context unless a specific correction can be
derived. In fact Burnham and Anderson (2004) go so far as to argue that it should be used
irrespective of sample size as it tends to the standard AIC when \emph{n} is large.
}

\value{
A numeric value with the AICc of the model
}

\references{
Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint

Hurvich, C. M. & Tsai, C.-L. (1989). \emph{Regression and Time Series Model Selection in Small Samples}. Biometrika, 76, 297-307

Burnham, K. P. & Anderson, D. R. (2002). \emph{Model Selection and Multimodel Inference: a Practical Information-theoretic Approach}. Springer

Burnham, K. P. & Anderson, D. R. (2004). \emph{Multimodel Inference: Understanding AIC and BIC in Model Selection}. Sociological Methods Research, 33, 261-304
}
\author{Maurice Berk \email{maurice.berk01@imperial.ac.uk}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  AICc(fit)
}