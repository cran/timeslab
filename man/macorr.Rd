
\name{macorr}
\alias{macorr}
\descrition{Calculate MA Autocorrelation Function}
\title{Calculate MA Autocorrelation Function}


\usage{macorr(beta,rvar=1,m)}

\arguments{
\item{beta}{Array containing MA coefficients $\Vbeta$.}
\item{rvar}{Real scalar containing error variance $(>0)$.}
\item{m}{Integer containing the number of autocorrelations to 
calculate $(>0)$.}
}

\value{
\item{macorr}{returns a list containing the following 
four items:
\item{var}{Real scalar containing the variance of the process.}
\item{corr}{Array of length {\code{m}} containing the 
autocorrelations.}
\item{ier}{Integer variable indicating whether or not the MA 
process is stationary (0 means yes, anything else means no).}
}
