
\name{ardt}
\alias{ardt}
\description{Simulate Data from an AR Process}
\title{Simulate Data from an AR Process}

\usage{ardt(alpha,rvar,n,seed=0)}

\arguments{
\item{alpha}{Array of length $p$ containing AR coefficients }
\item{rvar}{Real scalar containing error variance \code{sigma^2(>0)}.}
\item{n}{Integer \code{> p} containing the length of the
realization.}
\item{seed}{Real scalar containing the seed for the random
number generator.}
}
\value{
\item{ier}{Integer variable indicating whether or not the AR
process is stationary (0 means yes, anything else means no).}
\item{x}{Array of length \code{n} containing the realization.}
}

