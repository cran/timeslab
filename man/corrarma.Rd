
\name{corrarma}
\alias{corrarma}
\description{Calculate ARMA Parameters from Autocorrelations}
\title{Calculate ARMA Parameters from Autocorrelations}

\usage{corrarma(rho,r0,p,q,maxit=100,del=1.e-5)}

\arguments{
\item{rho}{Array of length ${\tt{p}}+{\tt{q}}$ containing the 
autocorrelations of the process.}
\item{r0}{Real scalar containing the variance of the process $(>0)$.}
\item{p}{Integer containing AR order $p(>0)$.}
\item{q}{Integer containing MA order $q(>0)$.}
\item{maxit}{Integer containing the maximum number of iterations to allow
in Wilson's algorithm.}
\item{del}{Real scalar containing convergence criterion $(>0)$.}
}


\value{
\item{alpha}{Array of length {\code{p}} containing AR coefficients.}
\item{beta}{Array of length {\code{q}} containing MA coefficients.}
\item{rvar}{Real scalar containing the error variance $\sigma^2$.}
\item{ier}{Integer variable containing an error/convergence indicator.  The
following values are possible:}
\item{ }{0 ~~ {\code{CORRARMA}} successfully found the ARMA parameters.}
\item{ }{1 ~~ A singular matrix was encountered trying to find AR 
parameters.}
\item{ }{2 ~~ Wilson's algorithm for finding the MA parameters didn't 
converge.}
}

