\name{ evaluate_expected.polynomial}
\alias{ evaluate_expected.polynomial} 
\title{ Evaluate the expected value of a multivariate polynomial }
\description{Evaluate the expected value of a multivariate polynomial 
assuming a specified non-central multivariate distribution.} 

\usage{ evaluate_expected.polynomial(poly,mu,sigma, envir='symmoments') }
\arguments{
  \item{ poly }{ either an object of class 'mpoly' or 'multipol', or a list with components for coefficients and powers.}
  \item{ mu }{ a vector of real numbers representing the mean of the multivariate distribution }
  \item{ sigma }{an vector giving an upper-triangular matrix representing the covariance matrix of the multivariate distribution }
  \item{ envir }{a character variable specifying the environment
containing the central moments needed for the calculation }

}
\details{ This function looks in the environment specified in the envir
argument for the central moments needed in the calculation.  The default
is the symmoments environment. The computation stops with an error 
message if a required moment is not found in envir.}

\value{expected value of the multivariate polynomial at the specified multivariate normal mean and covariance matrix}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, Journal of Statistical Software, 2010. 
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ See the evaluate_noncentral and make.all.moments functions.}

\examples{

# define a mpoly object for a multivariate polynomial and determine
# its expected value at specified mean and covariance matrix:
# note that all moments up to c(2,3,2) must exist in the symmoments
# environment.  Use make.all.moments(c(2,3,2)) if necessary.
# use library(mpoly) for first statement below.

# t0 <- mpoly(list(c(coef=3,x1=2),c(coef=2,x1=1,x2=3),c(coef=-4,z=2),c(coef=1,x1=1,x2=2,z=1)))
# evaluate_expected.polynomial(t0,c(1,2,3),c(1,0,0,1,0,1))



}

