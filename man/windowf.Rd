
\name{windowf}
\alias{windowf}
\description{Calculate Nonparamteric Spectral Density Estimate}
\title{Calculate Nonparamteric Spectral Density Estimate}

\usage{windowf(rho,R0,Q,ioptw,M,n,alpha=0.05)}

\arguments{
\item{rho}{Array of length {\code{M}} (if {\code{ioptw}} is between 1 and 5) or 
length ${\tt{n}}-1$ if {\code{ioptw}} is between 6 and 8 containing 
autocorrelations.}
\item{R0}{Real scalar containing the sample variance $(>0)$.}
\item{Q}{Integer containing the number of frequencies between 0 
and 1 at which to calculate spectra.}
\item{ioptw}{Integer containing the number of the window to be used in the 
estimation procedure as indicated by the following:}
\item{ }{1 ~~ Truncated periodogram}
\item{ }{2 ~~ Bartlett}
\item{ }{3 ~~ Tukey}
\item{ }{4 ~~ Parzen}
\item{ }{5 ~~ Bohman}
\item{ }{6 ~~ Daniell}
\item{ }{7 ~~ Bartlett--Priestley}
\item{ }{8 ~~ Parzen--Cogburn--Davis}
\item{M}{Integer $(>0)$ containing scale parameter.}n
\item{n}{(If either {\code{ioptw}} is between 6 and 8 or the factor for 
determining confidence intervals is desired.) Integer containing the length 
of the data set being analyzed.}n
\item{alpha}{Real scalar ($0<${\code{alpha}}$<1$) indicating the level
of confidence.}
}

\value{
\item{f}{Array of length $[{\tt{Q}}/2]+1$ containing the spectral 
estimator at the frequencies $(j-1)/{\tt{Q}},j=1,\ldots,[{\tt{Q}}/2]+1$.}
\item{c}{Real scalar variable that can be used to find 95\% confidence 
intervals for the true spectral density.  The interval at the $i$th 
frequency would be from {\code{f(i)/c}} to {\code{f(i)*c}}.}
}
