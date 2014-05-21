\name{toMatching}
\alias{toMatching}
\title{Convert representation of a phylogenetic tree as a moment L-matrix to matching form}
\description{Function converts a tree in moment format to matching format.

The input can be an \code{L-matrix} object, 
a square L matrix,
or an L matrix in reduced upper-triangular (vector) form.

The \code{toMatching} function sets its list output to class \code{L-matching},
which has 5 components, 
including the tree in matching format.
}
\usage{
toMatching(L, type = NULL, tip.label = NULL) 
}

\arguments{
  \item{L}{An \code{L-matrix} object, a square L matrix,
or an L matrix in reduced  upper-triangular (vector) form.}
  \item{type}{If object is not of class "L-matrix" and is a square L matrix,
then type should be "square'.  If it is an L matrix in upper triangular form,
type should be "ut". }
  \item{tip.label}{Character vector containing labels for tips. 
If null, labels default to "a"-"z" it at most 26; otherwise,
3-letter labels of the form "aaa", "aab",... }
}
\details{ An \code{L-matrix} object is a list with 5 components: "L" is the L-matrix in square form. "L.ut" is the L-matrix in upper-triangular form. "Newick" is the Newick representation of the tree.
"tip.label" is the character vector of tip labels. "tip.label.n" is the number of tips.}
\value{a matching representation of the phylogenetic tree corresponding to the input}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.

P.W. Diaconis and S. Holmes, Matchings and Phylogenetic Trees,
Proc. Natl. Acad. Sci., 1998,
95(25), 14600-14602}

\author{Kem Phillips <kemphillips@comcast.net}

\seealso{functions toMoment and toNewick}

\examples{
# create a Newick object
exam.Newick      <- "(((a,b),c),d);"
# convert to a moment L-matrix
exam.moment <- toMoment(exam.Newick)
# convert to matching format
exam.matching <- toMatching(exam.moment)
}
