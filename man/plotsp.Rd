
\name{plotsp}
\alias{plotsp}
\description{Form Plot of a (True or Sample) Spectral Density}
\title{Form Plot of a (True or Sample) Spectral Density}

\usage{plotsp(f,n,div,main=``Log Std Spectra'')}

\arguments{
\item{f}{Array of length $m=[{\tt{n}}/2]+1$ containing some spectral 
quantity the frequencies $(j-1)/{\tt{n}},j=1,\ldots,m$.}
\item{n}{Integer containing the length of the array {\code{f}}.}
\item{div}{Real scalar containing divisor of the array {\code{f}}.}
\item{main}{Main title of the resulting plot.}
}

\value{
\item{plotsp}{produces a plot of $\log{({\tt{f}}(j)/{\tt{div}})}$
versus $(j-1)/{\tt{n}}$ for $j=1,\ldots,[{\tt{n}}/2]+1$ with vertical 
scale running from $-6$ to 6.}
}
