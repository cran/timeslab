
\name{acf1}
\alias{acf1}
\description{Calculate Sample Autocorrelation Function}
\title{Calculate Sample Autocorrelation Function}

\usage{acf1(x,m=0)}
\arguments{
\item{x}{A time series.}
\item{m}{The number of lags at which to find the autocovariance 
function.}
}

\value{
\item{corr}{A vector containing the autocorrelations.}
\item{var}{A scalar containing the sample variance.}
}
\seealso{
  \code{\link{acf}}
  }
