
\name{pacf}
\alias{pacf}
\description{Calculate Sample Partial Autocorrelation Function}
\title{Calculate Sample Partial Autocorrelation Function}

\usage{pacf(x,m)}

\arguments{
\item{x}{Array of length $n$ containing the data.}
\item{m}{Integer $(0<{\code{m}}<n)$ containing the number of partial 
autocorrelations to find.}
}

\value{
\item{pacf}{returns an array of length {\code{m}} 
containing the partial autocorrelations.}
}
