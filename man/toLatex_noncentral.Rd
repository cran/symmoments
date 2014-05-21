
\name{toLatex_noncentral}
\alias{ toLatex_noncentral }   
\title{Compute a Latex expression for a noncentral moment   }
\description{Compute a Latex expression for a noncentral moment }
\usage{toLatex_noncentral(moment,envir='symmoments')}

\arguments{
  \item{moment}{ vector c(k1,... ,kn) specifying the moment X1**k1 *...* Xn**kn }
  \item{envir}{character variable specifying the environment that contains the required central moments}
}

\value{A text value giving the Latex representation of moment where X is multivariate normal}
\details{ All required moment objects must exist in the specified environment, with default 'symmoments'. 
However, if the sorted version of an unsorted moment exists, the tounsorted function is used to obtain it.}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ make.all.moments, tounsorted, callmultmoments and the method toLatex}
\examples{
# Compute the Latex representation of the 2-dimensional moment c(1,3) (not run).
# This requires that all moments up to c(1,3) exist in the symmoments environment.
# toLatex_noncentral(c(1,3))
}

