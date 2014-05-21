\name{ evaluate_noncentral}
\alias{ evaluate_noncentral} 
\title{ Evaluate a noncentral multivariate moment  }
\description{Computes the numerical value of a non-central moment at a specfied mean and specified covariance matrix} 
\usage{evaluate_noncentral(moment,mu,sigma,envir='symmoments') }
\arguments{
  \item{ moment }{ a vector of non-negative integers representing the non-central moment to be evaluated}
  \item{ mu }{ a vector of real numbers representing the mean of the multivariate normal distribution }
  \item{ sigma }{an upper-triangular matrix of covariance terms for the 
         multivariate normal distribution expressed as a vector at which the moment is to be evaluated}
  \item{ envir }{a character variable specifying the environment
containing the central moments needed for the calculation }

}
\details{This function looks in the environment specified in the envir
argument for the central moments needed in the calculation.  The default
is the symmoments environment.  All even moments less than or equal to 
the moment argument are required.  The computation stops with an error 
message if a required moment is not found in envir.}
\value{numeric value of the moment at the specified mean and covariance matrix}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, Journal of Statistical Software, 2010. 
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ See the evaluate.moment and make.all.moments functions.}

\examples{

# evaluate_noncentral(c(3,1,2),c(3,4,1),c(4,2,1,3,1,2))
# evaluates the expected value of X1**3 X2 X3**2 at mean c(1,2,3) 
# and at the following covariance matrix
#    4 2 1 
#    2 3 1 
#    1 1 2 

#  requires all moments up to c(3,1,2) to exist in the symmoments environment.
#  use  make.all.moments(c(3,1,2)) if necessary.

# use moments in the global environment:
# evaluate_noncentral(c(3,1,2),c(3,4,1),c(4,2,1,3,1,2),'.GlobalEnv')

}

