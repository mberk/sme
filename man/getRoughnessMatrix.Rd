% File sme/man/getRoughnessMatrix.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{getRoughnessMatrix}
\alias{getRoughnessMatrix}
\title{Extract the roughness matrix used in an SME model fit}
\description{Given an \code{sme} object, returns the roughness matrix which can be used to quickly
calculate the integrated squared second derivative of the mean and individual level curves}


\arguments{
  \item{object}{a fitted SME model object returned from the \code{sme} function}
}

\details{
The parameters of the SME model are estimated using maximum \emph{penalized} likelihood, where the
penalized likelihood is the standard likelihood with the addition of penalty terms corresponding to
the roughness of the mean and individual level curves. Typically the roughness of a curve is
quantified as its integrated squared second derivative. Green and Silverman (1994) show that, for a
natural cubic spline interpolating the vector of points $f$, there exists a \emph{roughness matrix}
$G$ such that the integrated squared second derivate is $f'Gf$ where $f'$ denotes $f$ transposed.
For details on constructing the matrix $G$, either refer to the original source of
Green and Silverman (1994) or it may prove easier to access Berk and Montana (2009) where they can
be found in the appendix.
}

\value{
The roughness matrix corresponding to the SME model fit
}

\references{
Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint

Berk, M. & Montana, G. (2009). \emph{Functional modelling of microarray time series with covariate curves}. Statistica, 2-3, 158-187

Green, P. J. & Silverman, B. W. (1994). \emph{Nonparametric Regression and Generalized Linear Models}. Chapman and Hall
}
\author{Maurice Berk \email{maurice.berk01@imperial.ac.uk}}
\examples{
  data(MTB)
  fit <- sme(MTB[MTB$variable==6031,c("y","tme","ind")])
  getRoughnessMatrix(fit)
}