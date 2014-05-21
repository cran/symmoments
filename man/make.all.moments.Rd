\name{make.all.moments}
\alias{ make.all.moments }   
\title{ Create all moments up to specified size in environment symmoments  }
\description{Create all central moment objects of a specified or smaller size in environment symmoments}
\usage{make.all.moments(moment,verbose=TRUE)}

\arguments{
  \item{moment}{ vector c(k1,... ,kn) specifying the highest moment to compute }
  \item{verbose}{if TRUE (default), the names of the moments are shown as the algorithm progresses; if FALSE, progress is not shown}

}

\value{all objects of class 'moment' up to the value given in moment are created in environment symmoments}

\details{Unsorted moments, those with exponents are not in numeric order,
are created in the symmoments environment using the \code{tounsorted}
function to transform from the sorted moment. If symmoments does not 
exist, the user is prompted to create it using symmoments <- new.env().

If the sorted moment does not exist, it is created.

Moments of lower dimension are not created; for example, if c(2,4) is  input, m20 is created, but m2 is not.


Moments are named mij..l, e.g., m136. 
If any exponent is greater than 9,
lower case letters and then upper case letters are used.  
For example,
m3bA is the name of the moment c(3,11,36).

The largest exponent allowed by this scheme is 9+26+26=61,

If an object with a name of this form exists but is not an object of class "moment", 
it is replaced (overwritten) by the moment object.
}
\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ callmultmoments, tounsorted}
\examples{
# Compute all moments up to c(3,3)
# First create the symmoments environment if it does not exist 
# symmoments <- new.env()
# make.all.moments(c(3,3))
}
