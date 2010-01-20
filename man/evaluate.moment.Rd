\name{ evaluate}
\alias{ evaluate.moment} 
\alias{ evaluate} 
\title{ Evaluate a multivariate moment  }
\description{Generic method for class moment to compute the numerical value of a moment at a specified covariance matrix from the output of callmultmoments} 
\usage{ \method{evaluate}{moment}(object,sigma) }
\arguments{
  \item{ object }{ an object of class 'moment' }
  \item{ sigma }{an upper-triangular matrix of covariance terms expressed as a vector at which the moment is to be evaluated}
}
\details{object is normally the output of a call to callmultmoment.
 This is a list with first component the moment itself,
 the second component the set of upper-triangular matrices representing the moment,
 and the third component containing their corresponding coefficients.
 This is an object of class 'moment'.
}
\value{numeric value of the moment at the specified covariance matrix}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, Journal of Statistical Software, 2010. 
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ callmultmoments and the simulate and toLatex methods from the symmoments package}

\examples{

evaluate(callmultmoments(c(1,2,3,4)),c(4,2,1,1,3,1,1,2,1,2))
# evaluates the moment at c(1,2,3,4) at the following covariance matrix
#    4 2 1 1
#    2 3 1 1
#    1 1 2 1
#    1 1 1 2

}

