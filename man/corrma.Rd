
\name{corrma}
\alias{corrma}
\description{Calculate MA Parameters from Autocorrelations}
\title{Calculate MA Parameters from Autocorrelations}

\usage{corrma(rho,r0,q,maxit=100,del=1.e-5)}

\arguments{
\item{rho}{Array of length {\code{q}} containing autocorrelations of 
lags ${\tt{1}},\ldots,{\tt{q}}$.}
\item{r0}{Real scalar containing the variance of the MA process.}
\item{q}{Integer containing order $q(>0)$.}
\item{maxit}{Integer containing the maximum number of iterations to use
in Wilson's algorithm $(>0)$.}
\item{del}{Real scalar containing the convergence criterion to use in
Wilson's algorithm $(>0)$.}
}


\value{
\item{corrma}{returns a list containing the following 
three items:}
\item{beta}{Array of length {\code{q}} containing MA coefficients.}
\item{rvar}{Real scalar containing the error variance $\sigma^2$ 
of the MA process.}
\item{ier}{Integer variable indicating whether or not Wilson's 
algorithm converged (0 means yes, 1 means no).}
}
