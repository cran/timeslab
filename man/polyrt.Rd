
\name{polyrt}
\alias{polyrt}
\description{Find the Roots of a Polynomial Given its Coefficients}
\title{Find the Roots of a Polynomial Given its Coefficients}

\usage{polyrt(coeffs,m=100,eps=1.e-6)}

\arguments{
\item{coeffs}{Array of length $p+1$ containing the coefficients of the
polynomial of degree $p$, ({\code{coeffs}}(i) is the coefficient of power
$i-1$).}
\item{m}{Integer containing the maximum number of iterations in the
procedure.}
\item{eps}{Real scalar indicating the convergence criterion.}
}

\value{
\item{polyrt}{returns a list containing the following 
items:}
\item{real}{Array of length $p$ containing the real parts of the roots.}
\item{imag}{Array of length $p$ containing the imaginary parts of the
roots.}
\item{ier}{Integer variable indicating whether or not the procedure
converged (0 means yes, 1 means no).}
}
