
\name{armapart}
\alias{armapart}
\title{Calculate ARMA Partial Autocorrelation Function}

\usage{armapart(alpha,beta,rvar,m)}

\arguments{
\item{alpha}{ Array containing AR coefficients}
\item{beta}{Array containing MA coefficients.}
\item{rvar}{Real scalar containing error variance \code{sigma^2(>0)}}.
\item{m}{Integer indicating the number of partial autocorrelations 
to calculate.}
}

\value{
\item{theta}{ Array of length {\code{m}} containing the partial 
autocorrelations.}
\item{ier}{Integer variable indicating whether or not the ARMA
process is stationary (0 means yes, anything else means no)}
}
