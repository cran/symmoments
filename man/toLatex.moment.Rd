\name{toLatex.moment}
\alias{toLatex.moment }   
\title{ LaTeX a multivariate moment }
\description{Computes a LaTeX representation sorted lexicographically of an object of class 'moment'} 
\usage{ \method{toLatex}{moment}(object,...) }
\arguments{
  \item{ object }{ an object of class 'moment', usually the output of callmultmoments}
  \item{...}{Included only for consistency with generic function}
}
\details{The first element of the result is the moment expressed as an expected value (E[...] =).
The remaining lines are the LaTex representation broken at appropriate intervals for printing.
(Individual terms for high dimensions will still overrun a printed line.)
Double backslashes are inserted where LaTeX requires a backslash.
These can be reset to single backslashes by writing the output to a file using the R function writeLines from the base package.
 }
\value{Character vector giving the LaTeX code for the symbolic moment}

\references{
Phillips K (2010). Symbolic Computation of the Central Moments of the
  Multivariate Normal Distribution. \emph{Journal of Statistical Software, Code Snippets},
  \bold{33}(1), 1--14. \url{http://www.jstatsoft.org/v33/c01/}.
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ callmultmoments and the evaluate method (symmoments) }

\examples{toLatex(callmultmoments(c(1,2,3)))}
