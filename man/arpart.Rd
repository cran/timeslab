
\name{arpart}
\alias{arpart}
\description{Calculate AR Partial Autocorrelation Function}
\title{Calculate AR Partial Autocorrelation Function}

\usage{arpart(alpha)}

\arguments{
\item{alpha}{Array of length $p$ containing AR coefficients $\Valpha$.}
}
\value{
\item{theta}{Array of length $p$ containing partial autocorrelations.}
\item{ier}{Integer variable indicating whether the zeros of the
characteristic polynomial corresponding to {\code{alpha}} are all
outside the unit circle (0 means they are, anything else means
they are not.)}
}

