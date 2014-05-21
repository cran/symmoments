\name{convert.mpoly}
\alias{ convert.mpoly }   
\title{ Convert between mpoly and list representations of multivariate polynomials  }
\description{Converts an mpoly object to a simple list representing a multivariate polynomial or a simple list to an mpoly object}
\usage{convert.mpoly(poly)}

\arguments{
  \item{poly}{ an mpoly object or a list giving powers and coefficients defining the polynomial}
 }

\value{if poly is of class 'mpoly', it is a list with two components shown below.  If poly is such a list, the value is the corresponding mpoly object}

\details{The list representation consists of 2 components: 'powers' is a matrix with each row representing the powers of X in one term of the multivariate polynomial. 'coeff' is a vector with each element being the coefficient of the corresponding term in powers}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ convert.multipol, evaluate.expected.polynomial, integrate.polynomial}

\examples{
# create an mpoly object here

# convert from mpoly to list representation
library(mpoly)
t0 <- mpoly(list(c(coef=3,x1=2),c(coef=2,x1=1,x2=3),
                   c(coef=-4,z=2),c(coef=1,x1=1,x2=2,z=1)))  
t1 <- convert.mpoly(t0)    
# convert from list representation back to an mpoly object
t2 <- convert.mpoly(t1) 

}

