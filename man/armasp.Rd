
\name{armasp}
\alias{armasp}

\description{Calculate ARMA Spectral Density Function}
\title{Calculate ARMA Spectral Density Function}

\usage{armasp(alpha,beta,rvar=1,Q=256)}

\arguments{
\item{alpha}{Array of length $p$ containing AR coefficients $\Valpha$.}
\item{beta}{Array of length $q$ containing MA coefficients $\Vbeta$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{Q}{Integer $(\ge\max(p,q))$ containing the number of
frequencies between 0 and 1 at which to calculate the spectral
density.}
}

\value{
\item{armasp}{returns an array {\code{f}} of length 
$m=[{\tt{Q}}/2]+1$ containing the values of the spectral density 
at the frequencies $(j-1)/{\tt{Q}},j=1,\ldots,m$.}
}
