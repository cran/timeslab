
\name{arfilt}
\alise{arfilt}
\description{Apply an AR Filter to a Matrix}
\title{Apply an AR Filter to a Matrix}

\usage{arfilt(alpha,rvar,x)}

\argument{
\item{alpha}{Array containing AR coefficients}
\item{rvar}{Real scalar containing error variance \code{sigma^2}}
\item{x}{Array containing the matrix to be filtered.}
}

\value{
\item{w}{Matrix containing the filtered version of {\code{x}}.}
\item{ier}{Integer variable indicating whether or not the AR
process is stationary (0 means yes, j $(>0)$ means j$^{th}$ 
partial outside $(-1,1)$).}
}

