
\name{diffeq}
\alias{diffeq}
\description{Find Future Values of a Difference Equation}
\title{Find Future Values of a Difference Equation}

\usage{diffeq(alpha,p,n,e)}

\arguments{
\item{alpha}{Array of length {\code{p}} containing AR coefficients $\Valpha$.}
\item{p}{Integer containing order $p(>0)$.}
\item{n}{Integer $(>{\tt{p}})$ containing the length of the realization.}
\item{e}{Array of length {\code{n}} containing the values for the
general difference equation.}
}

\value{
\item{x}{Array of length {\code{n}} containing the realization.}
}
