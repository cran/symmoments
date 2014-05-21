\name{convert.multipol}
\alias{ convert.multipol }   
\title{ Convert between multipol and list representations of multivariate polynomials  }
\description{Converts an multipol object to a simple list representing a multivariate polynomial or a simple list to an multipol object}
\usage{convert.multipol(poly)}

\arguments{
  \item{poly}{ a multipol object or a list giving powers and coefficients defining the polynomial}
}

\value{if poly is of class 'multipol', it is a list with two components described below.  If poly is such a list, the value is the corresponding multipol object}


\details{The list representation consists of 2 components: 'powers' is a matrix with each row representing the powers of X in one term of the multivariate polynomial. 'coeff' is a vector with each element being the coefficient of the corresponding term in powers}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ convert.multipol, evaluate.expected.polynomial, integrate.polynomial}

\examples{
# create an mpoly object to work with
library(mpoly)
t0 <- mpoly(list(c(coef=3,x1=2),c(coef=2,x1=1,x2=3),
                   c(coef=-4,z=2),c(coef=1,x1=1,x2=2,z=1))) 

# convert from mpoly to list representation
t1 <- convert.mpoly(t0)    
# convert from list representation to a multipol object
t2 <- convert.multipol(t1) 
# convert back to a list representation
t3 <- convert.multipol(t2) 
}

