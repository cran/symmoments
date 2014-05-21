\name{integrate.polynomial}
\alias{ integrate.polynomial}   
\title{Numerically integrate a multivariate polynomial }
\description{ Integrates a multivariate polynomial against a specified non-central multivariate distribution 
using 
ordinary integration by invoking the adaptIntegrate function from the cubature package.}
\usage{integrate.polynomial(poly,mu,sigma,lower=NULL,upper=NULL) }

\arguments{
  \item{poly}{ either an object of class 'mpoly' or 'multipol', or a list with two components for coefficients and powers. }
  \item{mu}{ a vector giving the mean of the multivariate distribution }
  \item{sigma}{ a square matrix giving the covariance matrix of the multivariate distribution }
  \item{lower}{ vectors of the lower limits of integration, one element for each dimension of the moment }
  \item{upper}{ vectors of the upper limits of integration, one element for each dimension of the moment }
}

\value{the expected value of the polynomial integrated against the multivariate normal distribution}

\details{Defaults for lower and upper are -/+ 6 times the standard deviations (square roots of diagonal elements of the covariance matrix).
If the polynomial is defined by a list, it has two components, coeff and powers.  
powers is a matrix.  Each row represents the powers for a term in the polynomial. 
coeff is a vector.  Each element is the coefficient of the corresponding power.
Example corresponding to example below:  list(coeff=c(3,2,-4,1),powers=matrix(c(2,0,0,1,3,0,0,0,2,1,2,1),ncol=3,byrow=TRUE))}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ evaluate.expected.polynomial, multmoments, evaluate, and simulate in symmoments}

\examples{
# define a mpoly object for a multivariate polynomial, and 
# determine its expected value at specified mean and covariance matrix:

# t0 <- mpoly(list(c(coef=3,x1=2),c(coef=2,x1=1,x2=3),c(coef=-4,z=2),c(coef=1,x1=1,x2=2,z=1)))


# integrate.polynomial(t0,c(1,2,3),matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,byrow=TRUE))


}

