
\name{filt}
\alias{filt}
\description{Apply a Linear Filter to a Vector}
\title{Apply a Linear Filter to a Vector}

\usage{filt(beta,beta0,x)}

\arguments{
\item{beta}{Array of length $p$ containing coefficients of lags 
$1,\ldots,p$.}
\item{beta0}{Real scalar containing the coefficient for lag 0.}
\item{x}{Array of length $n$ containing the data to be filtered.}
}

\value{
\item{filt}{returns an array of length $n-p$ 
containing the result of the filter. }
}
