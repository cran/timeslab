\name{acf}
\alias{acf}
\title{Calculate Sample Autocorrelation Function}
\description{Calculate Sample Autocorrelation Function}

\usage{acf(x,m=0)}

\arguments{
\item{x}{A time series.}
\item{m}{The number of lags at which to find the acf.}
}

\value{
\item{corr}{ A vector containing the autocorrelations.}
\item{var}{A scalar containing the sample variance.}
}
\note{acf and acf1 differ through the use of different internal
  algorithms but are otherwise equivalent.}
\seealso{
  \code{\link{acf1}}
  }
