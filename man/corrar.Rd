
\name{corrar}
\alias{corrar}
\description{Calculate AR Parameters from Autocorrelations}
\title{Calculate AR Parameters from Autocorrelations}

\usage{corrar(rho,R0,p)}

\arguments{
\item{rho}{Array of length {\code{p}} containing autocorrelations.}
\item{R0}{Real scalar containing sample variance $(>0)$.}
\item{p}{Integer containing the AR order $(>0)$.}
}
\value{ 
\item{rvar}{Real scalar variable containing error variance.}
\item{alpha}{Array of length {\code{p}} containing AR coefficients.}
}
