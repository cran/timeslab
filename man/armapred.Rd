
\name{armapred}
\alias{armapred}
\description{Calculate Exact Predictions for an ARMA Process}
\title{Calculate Exact Predictions for an ARMA Process}

\usage{armapred(x,alpha,beta,rvar,t1,t2,h1,h2)}

\arguments{
\item{x}{Array of length $n$ containing the realization to be
used in the prediction.}
\item{alpha}{Array containing AR coefficients $\Valpha$.}
\item{beta}{Array containing MA coefficients $\Vbeta$.}
\item{rvar}{Real scalar containing error variance $\sigma^2(>0)$.}
\item{t1,t2}{Integers $(1\le{\tt{t1}}\le{\tt{t2}}\le n)$
specifying the range of memories to be used.}
\item{h1,h2}{Integers $(1\le{\tt{h1}}\le{\tt{h2}})$ specifying
the horizons to be used.}
}

\value{
\item{xp}{Array of length $({\tt{t2}}-{\tt{t1}}+1)({\tt{h2}}-
{\tt{h1}}+1)$ containing predictors.}
\item{se}{Real scalar containing the prediction standard errors
for the predictors in the array {\code{xp}}.}
}
