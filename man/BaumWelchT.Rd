\name{BaumWelchT}
\alias{BaumWelchT}
\title{BaumWelchT HMM}
\usage{
  BaumWelchT(x, series.length, m = 2, Q, mu, S, nu = TRUE,
    model, initial.prob, maxiter = 500, overflow = 1e-09,
    num.inst = 1e-09, tol = 1e-05, dig = 3)
}
\arguments{
  \item{x}{matrix of observations.}

  \item{series.length}{length of independent blocks in
  matrix x.}

  \item{m}{number of states.}

  \item{Q}{transition matrix.}

  \item{mu}{means of the emissions in different states.}

  \item{S}{correlation matrix of the emissions.}

  \item{nu}{degrees of freedom (in all states).}

  \item{model}{a model matrix for parameter estimation.}

  \item{initial.prob}{state probability at the start of
  series.}

  \item{maxiter}{maximum number of iterations before
  returning.}

  \item{overflow}{the probability ratio in case of numeric
  overflow.}

  \item{num.inst}{threshold for no update (see E.step()).}

  \item{tol}{threshold to difference in log-likelihood
  before returning.}

  \item{dig}{numeric precision for parameter estimation.}
}
\value{
  the estimate of the HMM with emissions distributed as
  multi dimensional correlated T variables.
}
\description{
  Estimate the HMM with emissions distributed as
  multi-dimensional correlated T variables
}
\examples{
library(HMMt)
x <- c(rt(1000, df=3), rt(1000, df=3)+1)

# x has a t distribution with a jump at position 1001.
plot(x, type = 'l')

# Check the output of BaumWelchT.
BaumWelchT(x)

# See that it usually finds the transition.
lines(BaumWelchT(x)$ViterbiPath-1, col=2)
}
\author{
  Guillaume Filion. date: June 17, 2011.
}

