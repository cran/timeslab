
\name{arcorr}
\alias{arcorr}
\description{Calculate AR Autocorrelation Function}
\title{Calculate AR Autocorrelation Function}

\usage{arcorr(alpha,rvar=1,m=0)}

\arguments{
\item{alpha}{Array containing AR coefficients $\Valpha$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{m}{Integer containing the number of autocorrelations to
calculate $(\ge0)$.}
}

\value{
\item{var}{Real scalar containing the variance of the process.}
\item{corr}{Array of length {\code{m}} containing the 
autocorrelations}
\item{ier}{Integer variable indicating whether or not the AR 
process is stationary (0 means yes, anything else means no).}
}
