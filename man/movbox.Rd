
\name{movbox}
\alias{movbox}
\description{Form Quantities Needed for Moving Box Plot}
\title{Form Quantities Needed for Moving Box Plot}

\usage{movbox(x,k)}

\arguments{
\item{x}{Array of length $n$ containing the data.}
\item{k}{Integer indicating the length of the moving average smoother.}
}

\value{
\item{movbox}{returns a list containing the following 
three objects:}
\item{summ}{Matrix ($n\times 5$) with columns defined as follows:}
\item{ }{Column 1 is largest value in $[u4, u4 + 1.5*IQR]$,
          or $u4$ if none,}
\item{ }{Column 2 is the upper fourth (median of largest 
          ({\code{n}}+1)/2),}
\item{ }{Column 3 is the median,}
\item{ }{Column 4 is lower fourth (median of smallest
          ({\code{n}}+1)/2),}
\item{ }{Column 5 is smallest value in $[l4, l4 - 1,5*IQR]$,
          or $l4$ if none,}
\item{ }{where $IQR$ is the interquartile range.}
\item{inds}{Integer array of indices.}
\item{outs}{Real array containing values corresponding to {\code{inds}} 
that fall outside $[l4 - 1.5*IQR, u4 + 1.5*IQR]$. }
}
