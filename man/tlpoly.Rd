
\name{tlpoly}
\alias{tlpoly}
\description{Evaluate a Polynomial at a Vector of Points}
\title{Evaluate a Polynomial at a Vector of Points}

\usage{tlpoly(coeffs,x)}

\arguments{
\item{coeffs}{Array of length $p+1$ containing the coefficients of the 
polynomial of degree $p$, ({\code{coeffs}}(i) is the coefficient of power
$i-1$).}
\item{x}{Array containing the points at which to evaluate the polynomial.}
}

\value{
\item{tlpoly}{returns an array containing the values 
of the polynomial.}
}
