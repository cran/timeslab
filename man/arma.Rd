
\name{arma}
\alias{arma}
\description{Form Plots Illustrating Patterns in ARMA Processes}
\title{Form Plots Illustrating Patterns in ARMA Processes}

\usage{arma(alpha,beta,x,iopt,p,q,rvar,n,m,seed=0)}

\arguments{
\item{alpha}{Array of length {\code{p}} containing AR coefficients 
$\Valpha$.}
\item{beta}{Array of length {\code{q}} containing MA coefficients 
$\Vbeta$.}
\item{x}{An ARMA process.}
\item{iopt}{Integer indicating which part of the ARMA
process is to be simulated.  ${\tt{iopt}=0}$ means to use input
alpha, beta and x,  ${\tt{iopt}=1}$ means to use input alpha
and beta and simulate x, and ${\tt{iopt}=2}$ means to simulate
alpha, beta and x.}
\item{p}{Integer containing order $p$ of the array $\Valpha$.}
\item{q}{Integer containing order $q$ of the array $\Vbeta$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{n}{Length of the realization.}
\item{m}{Number of autocorrelations to be calculated.}
\item{seed}{Real scalar containing the seed for the random
number generator.}
}

\value{
\item{ }{ \code{arma} returns plots illustrating patterns in ARMA 
processes and a list containing the following three items:}
\item{alpha}{Array containing the AR coefficients.}
\item{beta}{Array containing the MA coefficients.}
\item{x}{The realization.}
}
\example{
> arma(alpha,beta,x,2,5,4,0.2,500,20,0)         
}
