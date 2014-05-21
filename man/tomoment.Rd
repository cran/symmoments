\name{toMoment}
\alias{toMoment}
\title{Converts a tree from Newick or matching to moment format}
\description{Converts a tree from Newick or matching to moment format}
\usage{
toMoment(inputobject, tip.label = NULL)
}

\arguments{
  \item{inputobject}{a tree in Newick format or a \code{matching} object defined in the \pkg{ape} package}
  \item{tip.label}{rearranged labels for tips; these must be the original labels}
}
\details{The L-matrix class consists of $5$ components:
"L" is the L-matrix in square form. "L.ut" is the L-matrix in upper-triangular form.
"Newick" is the Newick representation of the tree.
"tip.label" is the character vector of tip labels. "tip.label.n" is the number of tips.
}

\value{a moment L-matrix corresponding to the input phylogentic tree object}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.

J. Felsenstein, The {N}ewick tree format, 1990,
http://evolution.genetics.washington.edu/phylip/newicktree.html 

P.W. Diaconis and S. Holmes, Matchings and Phylogenetic Trees,
Proc. Natl. Acad. Sci., 1998,
95(25), 14600-14602
}

\author{Kem Phillips <kemphillips@comcast.net>}



\seealso{functions toNewick and toMatching}

\examples{
# create a Newick object
exam.Newick      <- "(((a,b),c),d);"
# convert to a moment L-matrix
exam.moment <- toMoment(exam.Newick)
# convert to matching object
exam.matching <- toMatching(exam.moment)
# convert back to moment object
backto.moment <- toMoment(exam.matching)
}

