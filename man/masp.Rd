

\name{masp}
\alias{masp}
\description{Calculate MA Spectral Density Function}
\title{Calculate MA Spectral Density Function}

\usage{masp(beta,rvar=1,Q=256)}

\arguments{
\item{beta}{Array of length $q$ containing coefficients $\Vbeta$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{Q}{Integer $(>q)$ containing the number of frequencies between 0
and 1 at which to calculate the spectral density. }
}

\value{
\item{masp}{returns the array {\code{f}} of length 
$m = [{\tt{Q}}/2]+1$ containing the MA spectra at frequencies 
$(j-1)/{\tt{Q}},j=1\ldots,m$. }
}
