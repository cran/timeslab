
\name{arsp}
\alias{arsp}
\description{Calculate AR Spectral Density}
\title{Calculate AR Spectral Density}

\usage{arsp(alpha,rvar=1,Q=256)}

\arguments{
\item{alpha}{Array of length $p$ containing AR coefficients $\Valpha$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{Q}{Integer $(>{\code{p}})$ containing the number of
frequencies between 0 and 1 at which to calculate the spectral
density.}
}

\value{
\item{arsp}{returns the array {\code{f}} of length 
$m=[{\tt{Q}}/2]+1$ containing the values of the spectral density at the 
frequencies $(j-1)/{\tt{Q}},j=1,\ldots,m$.}
}

