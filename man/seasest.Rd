
\name{seasest}
\alias{seasest}
\description{Calculate Box-Jenkins Estimates for a Seasonal ARIMA Model}
\title{Calculate Box-Jenkins Estimates for a Seasonal ARIMA Model}

\usage{seasest(y,ords,coeffs,lags,back,maxit=50,eps=0.000001)}

\arguments{
\item{y}{Array of length $n$ containing the data.}
\item{ords}{An array of length 5 containing the full and subset AR 
orders, followed by the full and subset MA orders, followed by a 1 if a 
constant term is in the model or a 0 if it is not.}
\item{coeffs}{An array containing starting values for the coefficients 
that are included in the model in the order full AR, subset AR, full MA, 
subset MA and the mean of {\code{y}}.}
\item{lags}{An array containing the lags (if any) in the model.  If both 
the subset AR and MA orders are zero, no array called {\code{lags}} need 
be formed, but an argument must be included.}
\item{back}{An integer containing the number of back forecasts to used in 
determining initial values in the recursion used in evaluating the sum of 
squares of residuals functions $(\ge0)$.}
\item{maxit}{An integer containing the number of iterations to allow 
in the estimation procedure.  If {\code{maxit}} is negative, then 
--{\code{maxit}} iterations are allowed and the values of the coefficients 
for the successive iterations are displayed on the screen.  If {\code{maxit}} 
is 1 then {\code{SEASEST}} only evaluates {\code{rvar}} and {\tsl{sds}}.}
\item{eps}{Real scalar containing a convergence criterion.  If the 
maximum value of successive iterates differs by less than {\code{eps}}, 
then {\code{SEASEST}} judges that the algorithm has converged.
}}

\value{
\item{seasest}{returns a list containing the following 
five elements:}
\item{coeffs}{Array containing the final values reached for the parameters 
in the iterative process.  {\code{coeffs}} is not changed from input if 
{\code{maxit}}$=1$.}
\item{e}{Array of length $n$ containing the one step ahead prediction 
errors corresponding to the $n$ values of {\code{x}}.}
\item{ier}{An integer variable indicating whether or not convergence was
achieved (0 means yes, 1 means no), if a singular matrix was encountered 
(2), or whether the algorithm could not continue even though convergence 
was no reached (3 or 4).  If this final alternative happens, different 
starting values or convergence criteria may lead to convergence.}
\item{rv}{Real scalar containing an estimate of the error variance.}
\item{se}{An array containing the standard errors of the estimates. }
}
