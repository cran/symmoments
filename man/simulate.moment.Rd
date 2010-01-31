\name{simulate.moment}
\alias{ simulate.moment } 
\title{ Method to compute a multivariate moment using Monte Carlo integration  }
\description{Computes a multivariate normal moment by Monte Carlo integration}
\usage{\method{simulate}{moment}(object,nsim,seed,Mean,Sigma,...)}
\arguments{
  \item{object}{object of class 'moment' representing E[X1**k1,...,Xn**kn]}
  \item{nsim}{the number of samples to generate in computing the integral}
  \item{seed}{integer for random number generator (set.seed)}
  \item{Mean}{the mean of (X1,...,Xn)}
  \item{Sigma}{covariance of (X1**k1,...,Xn**kn), dimension nXn, expressed as a vector by row }
  \item{...}{Included only for consistency with generic function}
}
\value{Approximate value of the moment}
\references{Rizzo ML (2008). Statistical Computing with R. Chapman & Hall/CRC}
\author{Kem Phillips <kemphillips@comcast.net>}
\note{Non-central moments can be approximated by specifying Mean. For central moments, set Mean to a vector of 0s.

The mvtnorm package must be loaded for the function rmvnorm.}

\seealso{ callmultmoments and the methods toLatex and evaluate from symmoments }
\examples{
# Using 10000 samples, estimate the central moment for the moment c(2,4) at the covariance matrix
#  2 1
#  1 4

# and mean (0,0)
library(mvtnorm)
simulate(callmultmoments(c(2,4)),10000,NULL,c(0,0),c(2,1,1,4))
} 

