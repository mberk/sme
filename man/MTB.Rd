% File sme/man/MTB.Rd
% Part of the sme package for R
% Distributed under GPL 3 or later: see sme/LICENCE

\name{MTB}
\alias{MTB}
\title{Ten example gene transcripts from an experiment investigating the \emph{in vitro}
transcriptional response to \emph{M.Tuberculosis}}
\description{
    The \code{MTB} data frame has 290 rows and 4 columns.
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
        \item{variable}{
            a factor indicating which gene transcript the observations in \code{y} belong to
        }        
    }
}
\source{
    Berk, M. (2012). \emph{Smoothing-splines Mixed-effects Models in R}. Preprint
}

