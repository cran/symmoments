\name{callmultmoments}
\alias{ callmultmoments }   
\title{ Compute multivariate moment symbolically  }
\description{Computes a multivariate normal moment by initializing variables, calling multmoments, and constructing output}
\usage{callmultmoments(moment)}

\arguments{
  \item{moment}{ vector c(k1,... ,kn) specifying the moment X1**k1 *...* Xn**kn }
}

\value{A object of class 'moment', which is a list with three components:

\item{moment}{the input moment vector}

\item{representation}{a matrix containing the representation in terms of upper-triangular matrices}

\item{coefficients}{the coefficients corresponding to the rows of the representation}

If the sum of the exponents is odd, returns -1 and prints "Sum of powers is odd. Moment is 0."

If any exponent is negative, returns -2 and prints "All components of the moment must be non-negative."

If any exponent is not an integer, returns -3 and prints "All components of the moment must be integers."
}
\details{Each row of the representation gives the exponents for a single product of covariance terms.
For example, (1,2,0) represents S11**1 S12** S22**0, where the Sij are the covariances.  
The full moment is the sum of these terms multiplied by their respective coefficients.
If the sum of the exponents is odd, the moment is 0.  
}
\references{
Phillips K (2010). Symbolic Computation of the Central Moments of the
  Multivariate Normal Distribution. \emph{Journal of Statistical Software, Code Snippets},
  \bold{33}(1), 1--14. \url{http://www.jstatsoft.org/v33/c01/}.
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ multmoments and the methods toLatex, evaluate, and simulate in symmoments}
\examples{
# Compute the moment for the 4-dimensional moment c(1,2,3,4):

m.1234 <- callmultmoments(c(1,2,3,4))
}

