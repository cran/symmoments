\name{print.moment}
\alias{print.moment }   
\title{ Print the representation of a multivariate moment }
\description{Prints an object of class 'moment'} 
\usage{ \method{print}{moment}(x,...) }
\arguments{
  \item{ x }{ an object of class 'moment', usually the output of callmultmoments}
  \item{...}{Included only for consistency with generic function}
}
\details{ Prints the moment as E[X1**k1 X2**k2 ...]:  followed by the lines of the representation
with the corresponding coefficient attached}

\references{
Phillips K (2010). Symbolic Computation of the Central Moments of the
  Multivariate Normal Distribution. \emph{Journal of Statistical Software, Code Snippets},
  \bold{33}(1), 1--14. \url{http://www.jstatsoft.org/v33/c01/}.
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ callmultmoments (symmoments) }

\examples{print(callmultmoments(c(1,2,3)))}
