
\name{swp}
\alias{swq}
\description{Sweep a Matrix}
\title{Sweep a Matrix}

\usage{swp(a,k1,k2)}

\arguments{
\item{a}{The matrix to be swept.}
\item{k1,k2}{{\code{a}} is swept on the diagonals {\code{k1}} through 
{\code{k2}}.}
}

\value{
\item{ }{\code{swp} returns a list containing the following two items:}
\item{A}{Matrix that results from sweeping {\code{a}}.}
\item{ier}{Integer variable containing an indicator of whether or not a zero 
diagonal was encountered during the sweeping (0 means no, 1 means yes).}
}
