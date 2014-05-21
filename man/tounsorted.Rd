\name{tounsorted}
\alias{ tounsorted }   
\title{Compute an unsorted central moment object from a sorted object }
\description{Produces an unsorted central moment object from a sorted object of class "moment".
 
Unsorted moments are those with exponents not in numeric order, e.g., m312.}
\usage{tounsorted(moment,sorted.moment)}

\arguments{
  \item{moment}{ unsorted moment to obtain moment is in vector form, eg, c(3,1,2)}
  \item{sorted.moment}{ sorted moment to use in obtaining unsorted moment
 }
}

\value{A object of class 'moment', which is a list with three components:

\item{moment}{the input moment vector}

\item{representation}{a matrix containing the representation in terms of upper-triangular matrices}

\item{coefficients}{the coefficients corresponding to the rows of the representation}

}

\details{The unsorted moment is obtained by resorting the rows and columns of the sorted moment sucessively. 
}
\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ multmoments}
\examples{
# obtain moment m312 from m123
tounsorted(c(3,1,2),callmultmoments(c(1,2,3)))

}

