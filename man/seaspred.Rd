
\name{seaspred}
\alias{seaspred}
\description{Calculate Box-Jenkins Forecasts}
\title{Calculate Box-Jenkins Forecasts}

\usage{seaspred(x,ords,coeffs,lags,rvar,tf,tl,hl,conf)}

\arguments{
\item{x}{Array of length $n$ containing the data.}
\item{ords}{An array of length 8 containing the full and subset AR orders, 
followed by the full and subset MA orders, followed by a 1 if a constant 
term is in the model or a 0 if it is not, followed by the number of first 
differences in the model, the number of $S$th differences in the model, and 
finally the value of $S$.}
\item{coeffs}{Values for full AR, subset AR, full MA and subset MA 
coefficients, followed by the constant if there is one and the values of $m$ 
and $\lambda$ for the power transform.}
\item{lags}{Array containing the lags (if any) in the model.  If both the 
subset AR and MA orders are zero, no array called {\code{lags}} need be formed, 
but an argument must be included.}
\item{rvar}{Real scalar containing the value of noise variance.}
\item{tf,tl}{Integers containing the prediction origins to use.  The values 
must be at least ${\rm{maxp}}+{\rm{maxq}}+1$ (maxp and maxq are the largest 
AR and MA lags in the expanded version of the model) and at most $n$, and 
{\code{tf}} must be less than or equal to {\code{tl}}.}
\item{hl}{Integer containing the maximum number of steps ahead to forecast 
from each origin.}
\item{conf}{Real scalar containing the confidence level for the probability 
limits to be placed on the forecasts $(0<{\tt{conf}}<1)$.}
}

\value{
\item{seaspred}{returns a list containing the following 
six items:}
\item{xp}{Array of length $({\tt{tl}}-{\tt{tf}}+1){\tt{hl}}$ containing 
the forecasts.  The first {\code{hl}} elements are from origin {\tt{tf}}, 
the next {\code{hl}} are from origin {\code{tf+1}}, etc.}
\item{xpl}{Array containing the lower probability limits on the corresponding 
elements of {\code{xp}}.}
\item{xpu}{Array containing the upper probability limits on the corresponding
elements of {\code{xp}}.}
\item{ier}{Integer variable indicating whether {\code{SEASPRED}} finished
without error (0), or an illegal power transform was requested (1).}
\item{check}{?*?*?*?*?*?.}
\item{xx}{?*?*?*?*?*?.}
}
