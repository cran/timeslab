
\name{madt}
\alias{madt}
\description{Simulate Data from an MA Process}
\title{Simulate Data from an MA Process}

\usage{madt(beta,rvar,n,seed=0)}

\arguments{
\item{beta}{Array of length $q$ containing the MA coefficients $\Vbeta$.}n
\item{rvar}{Real scalar containing the error variance $\sigma^2(>0)$.}
\item{n}{Integer $(>q)$ containing the desired number of observations.}
}

\value{
\item{madt}{returns an array {\code{x}} of length {\tsl{n}} 
containing the realization.}
}
