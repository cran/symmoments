\name{toNewick}
\alias{toNewick}
\title{
convert representation of phylogenetic tree as a moment L-matrix to Newick form 
}
\description{function converts a tree in moment format to Newick format.

The input can be an \code{L-matrix} object, 
a square $L$ matrix,
or an $L$ matrix in reduced upper-triangular (vector) form.

The \code{toNewick} function sets its list output to class \code{L-Newick},
which has $5$ components, 
including the tree in Newick format.
}
\usage{
toNewick(L, type = NULL, tip.label = NULL)
}

\arguments{
  \item{L}{L can be an \code{L-matrix} object, a square $L$ matrix,
or an $L$ matrix in reduced upper-triangular (vector) form.}
  \item{type}{if L is not a \code{L-matrix} object, either 'square' or 'ut' as listed above}
  \item{tip.label}{Character vector containing labels for tips. 
If null, labels default to "a"-"z" it at most 26; otherwise,
3-letter labels of the form "aaa", "aab",... }
}
\details{ An \code{L-matrix} object is a list with 5 components: "L" is the L-matrix in square form. "L.ut" is the L-matrix in upper-triangular form. "Newick" is the Newick representation of the tree.
"tip.label" is the character vector of tip labels. "tip.label.n" is the number of tips.}

\value{a Newick representation of the phylogenetic tree corresponding to the input}

\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010.

 J. Felsenstein, The {N}ewick tree format,
 1990,
 http://evolution.genetics.washington.edu/phylip/newicktree.html 

P.W. Diaconis and S. Holmes, Matchings and Phylogenetic Trees,
Proc. Natl. Acad. Sci., 1998,
95(25), 14600-14602
}

\author{Kem Phillips <kemphillips@comcast.net>}


\seealso{functions toMoment and toMatching}

\examples{
# create a Newick object
exam.Newick      <- "(((a,b),c),d);"
# convert to a moment L-matrix
exam.moment <- toMoment(exam.Newick)
# convert back to Newick format
backto.Newick <- toNewick(exam.moment)
}
