
\name{corrdt}
\alias{corrdt}
\description{Simulate Data Having Specified Autocorrelations}
\title{Simulate Data Having Specified Autocorrelations}

\usage{corrdt(rho,r0,n,seed=0)}

\arguments{
\item{rho}{Array containing autocorrelations..}
\item{r0}{Real scalar containing the variance of the process $(>0)$.}
\item{n}{Length of the desired realization.}
\item{seed}{Real scalar containing the seed for the random number generator.}
}

\value{
\item{corrdt}{returns a realization of length {\code{n}} from a
Gaussion process having variance {\code{ro}} and autocorrelations {\code{rho}}}
}
