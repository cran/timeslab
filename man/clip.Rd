
\name{clip}
\alias{clip}
\description{Clip a Vector Above and Below}
\title{Clip a Vector Above and Below}

\usage{clip(x,low=-1.e20,up=1e.20)}

\arguments{
\item{x}{A vector.}
\item{low}{The value that will replace any element of {\code{x}} that is
less than {\code{low}}.}
\item{up}{The value that will replace any element of {\code{x}} that is 
greater than {\code{up}}.}
}

\value{
\item{x}{The original {\code{x}} with any element less that {\code{low}} 
replaced by {\code{low}} and any element greater than {\code{up}} replaced 
by {\code{up}}.}
}

