
\name{dtarma}
\alias{dtarma}
\description{Calculate Exact ARMA MLE's}
\title{Calculate Exact ARMA MLE's}

\usage{dtarma(x,alpha,beta,maxit=200,eps=.0001)}

\arguments{
\item{x}{Array containing the data to be used in the estimation 
of the procedure.}
\item{alpha}{Array containing the starting value for the AR 
coefficients.}
\item{beta}{Array containing the starting value for the MA 
coefficients.}
\item{maxit}{Integer $(1\le{\tt{maxit}}\le 500)$ containing 
the maximum number of iterations in the optimization procedure.}
\item{eps}{Real scalar containing the convergence criterion.}
}

\value{
\item{dtarma}{returns a list containing the following 
six items:
\item{ier}{Integer variable containing termination information 
(0 means convergence, $>1$ means some convergence error.)}
\item{alpha}{Array containing the values of the AR coefficients at 
termination.}
\item{beta}{Array containing the values of the MA coefficients at 
termination.}
\item{rvar}{Real scalar variable containing the value of the error 
variance at termination.}
\item{m2ll}{Real scalar variable containing the value of $-2$ times
the log likelihood evaluated at the output values of the parameters.}
\item{var}{Real scalar containing a measure of degree of convergence 
at termination.}
}
