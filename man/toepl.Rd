
\name{toepl}
\alias{toepl}
\description{Form Symmetric Toeplitz Matrix Given its First Row}
\title{Form Symmetric Toeplitz Matrix Given its First Row}

\usage{toepl(R,R0,M)}

\arguments{
\item{R}{Array of length ${\tt{n}}-1$ containing the second through 
{\code{n}}th elements of the first row of the Toeplitz matrix.}
\item{R0}{Real scalar containing the value for the diagonal of the Toeplitz 
matrix.}
\item{M}{Size of the resulting matrix.}
}

\value{
\item{G}{The desired Toeplitz matrix.}
}
