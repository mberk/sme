% File sme/man/inflammatory.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{inflammatory}
\alias{inflammatory}
\title{Example gene transcript from an experiment on children suffering from an inflammatory
condition}
\description{
    The \code{inflammatory} data frame has 318 rows and 3 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
        \item{y}{
            a numeric vector of observed gene expression levels
        }
        \item{tme}{
            a numeric vector of time points corresponding to the observations in \code{y}
        }
        \item{ind}{
            a factor indicating which subject the observations in \code{y} belong to
        }
    }
}
\source{
    Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
}
