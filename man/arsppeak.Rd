
\name{arsppeak}
\alias{arsppeak}
\description{Find Peak Frequencies in AR Spectra}
\title{Find Peak Frequencies in AR Spectra}

\usage{arsppeak(alpha,rvar,n,start=0)}

\arguments{
\value{
\item{alpha}{Array of length \code{p} containing AR coefficients
$\Valpha$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{n}{If \code{ARSPPEAK} is being used for estimation purposes, 
{\code{n}} is an integer containing the length of the realization that 
was used to estimate the parameters of the process.  If the
parameters are the true values, let ${\code{n}}=1$.}
\item{start}{An optional argument that is a real scalar
containing a starting value $(0<{\tt{start}}<.5)$ for the
maximum finding procedure.}
}

\value{
\item{arsppeak}{returns a list containing the following 
three items:
\item{ier}{An integer indicating whether or not ARSPPEAK was
successful in finding a peak.  The possible values of
{\code{ier}} are:}
\item{ }{0 ~~ ARSPPEAK was successful in finding a peak.}
\item{ }{1 ~~ ARSPPEAK judged that the AR spectral density
has no relative maxima.}
\item{ }{2 ~~ A zero second derivative was encountered.}
\item{ }{3 ~~ The maximum finder converged to frequency 0
or 0.5.}
\item{ }{4 ~~ The maximum finder didn't converge.}
\item{peakf}{If ${\tt{ier}}=0$, the peak frequency.}
\item{se}{If ${\tt{ier}}=0$, the standard error of the peak
frequency estimator.}
}

