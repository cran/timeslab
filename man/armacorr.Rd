
\name{armacorr}
\alias{armacorr}
\description{Calculate ARMA Autocorrelation Function}
\title{Calculate ARMA Autocorrelation Function}

\usage{armacorr(alpha,beta,rvar=1,m)}

\arguments{
\item{alpha}{Array of length $p$ containing AR coefficients 
$\Valpha$.}
\item{beta}{Array of length $q$ containing MA coefficients 
$\Vbeta$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{m}{Integer $(\ge\max(p,q))$ containing the number
of autocorrelations to calculate.}
}
\value{
\item{var}{Real scalar containing variance of process.}
\item{corr}{Array of length {\code{m}} containing autocorrelations
$\rho(1), \ldots, \rho(m)$}
\item{ier}{Integer variable indicating whether or not the ARMA
process is stationary (0 means yes, anything else means no).}
}

