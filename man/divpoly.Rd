
\name{divpoly}
\alias{divpoly}
\description{Divide Two Polynomials}
\title{Divide Two Polynomials}

\usage{divpoly(num,den,n)}

\arguments{
\item{num}{Array containing coefficients of the numerator polynomial whose 
zeroth coefficient is one.}
\item{den}{Array containing coefficients of the denominator polynomial
with zeroth coefficient equal to one.}
\item{n}{Integer containing the order of resulting polynomial. }
}

\value{
\item{divpoly}{returns the array {\code{ratio}} of length {\tsl{n}} 
containing the coefficients of the polynomial that results from the 
division the original two polynomials.}
}
