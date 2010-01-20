\name{ multmoments}
\alias{ multmoments }   

\title{ Recursive function to compute a multivariate moment }
\description{Called by callmultmoments to compute representation of a multivariate normal moment using recursive algorithm} 
\usage{multmoments(moment,current.matrix,current.cell,moment.rep,row_col)}
\arguments{
  \item{ moment}{vector c(k1,...,kn) specifying the moment X1**k1 *...* Xn**kn}
  \item{ current.matrix }{ upper-triangular integer matrix under consideration in recursion}
  \item{ current.cell }{ cell in current matrix under consideration in recursion}
  \item{ moment.rep }{ current set of representations; mult.moments adds each satisfying matrix to moment.rep}
  \item{ row_col }{ matrix giving rows and columns for square matrix for each cell}
}
\details{Each row of the representation gives the exponents for a single product of covariance terms. For example, (1,2,0) represents S11**1 S12**2 S22**0, where the Sij are the covariances.

This function would normally only be called by callmultmoments.}

\value{ moment representation, moment.rep, augmented with additional representations}
\references{ K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010. 
}
\author{Kem Phillips <kemphillips@comcast.net>}
\seealso{ callmultmoments (symmoments) }

