
\name{wn}
\alias{wn}
\description{Simulate White Noise Data}
\title{Simulate White Noise Data}

\usage{wn(seed,n,dist=1)}

\arguments{
\item{seed}{Real scalar containing the seed for the random number generator.}
\item{n}{Integer containing the length of the desired realization.}
\item{dist}{Integer containing the number of the distribution to use based
on the following values:}
\item{ }{1 ~~ N(0,1)}
\item{ }{2 ~~ U(0,1)}
\item{ }{3 ~~ Unit exponential}
\item{ }{4 ~~ Logistic}
\item{ }{5 ~~ Standard Cauchy}
\item{ }{6 ~~ Extreme value}
\item{ }{7 ~~ Lognormal}
\item{ }{8 ~~ Double exponential}
}

\value{
\item{wn}{returns a realization of white noise from the specified
distribution.}
}
