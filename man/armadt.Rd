
\name{armadt}
\alias{armadt}
\description{Simulate Data from an ARMA Process}
\title{Simulate Data from an ARMA Process}

\usage{armadt(alpha,beta,rvar,n,seed=0)}

\arguments{
\item{alpha}{Array of length $p$ containing AR coefficients 
$\Valpha$.}
\item{beta}{Array of length $q$ containing MA coefficients 
$\Vbeta$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{n}{Integer $(>\max(p,q))$ containing the length of
the realization.}
\item{seed}{Real scalar containing the seed for the random
number generator.}
}

\value{
\item{x}{Array of length {\code{n}} containing the realization.
\item{ier}{Integer variable indicating whether or not the ARMA
process is stationary (0 means yes, anything else means no).}
}
