\name{symmoments-package}
\alias{ symmoments-package}
\alias{ symmoments}
\docType{package}
\title{Symbolically compute multivariate normal central moments}

\description{
Symbolically computes multivariate normal central moments  
E[X1**2 *...* Xn**n], where (X1,...,Xn) ~ N(0,S), 
in terms of S elements.

Produces Latex code for the moment.

Computes numerical moments at specified covariance matrices.
}

\details{
\tabular{ll}{
Package: \tab symmoments\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-01-20\cr
License: \tab GPL 2\cr
LazyLoad: \tab yes\cr}

A representation of a central moment of the multivariate normal distribution, given by a 
positive integer vector c(k1,k2,...,kn), is obtained from the function callmultmoments.
This function initializes variables and calls the function multmoments which determines 
a representation of a multivariate moment using a recursive algorithm.   
The representation is given class 'moment'.

The print method prints the representation of a multivariate moment.

The toLatex method uses the output of callmultmoments to determine the LaTeX code 
for the moment sorted lexicographically.  

The generic evaluate method uses the output of callmultmoments to determine
the value of the moment for a specified covariance matrix.

The simulate method is used to approximate a (possibly non-central) moment 
using Monte Carlo integration.

}
\author{Maintainer: Kem Phillips <kemphillips@comcast.net>}
\references{
K Phillips, Symbolic Computation of the Central Moments of the Multivariate Normal Distribution, 
Journal of Statistical Software, 2010. 
}

\note{The mvtnorm package must be loaded for the simulate method.}

\examples{
# Compute the moment for the 4-dimensional moment c(1,2,3,4):
callmultmoments(c(1,2,3,4)) 

# Print the representation of the 4-dimensional moment c(1,2,3,4):
print(callmultmoments(c(1,2,3,4))) 

# Compute the LaTeX representation of the central moment c(1,2,3,4):
toLatex(callmultmoments(c(1,2,3,4)))

# Write the LaTeX representation to a file using the standard R function (not run):
# writeLines(callmultmoments(c(1,2,3,4))),con="yourfilename", sep = "\n")

# evaluate the moment c(1,2,3,4) at the following variance-covariance matrix
  #  4 2 1 1
  #  2 3 1 1
  #  1 1 2 1

evaluate(callmultmoments(c(1,2,3,4)),c(4,2,1,1,3,1,1,2,1,2))

# Using 10000 samples, estimate the central moment for c(2,4) at the covariance matrix (not run)
#  2 1
#  1 4

# and mean (0,0)
# library(mvtnorm)
# simulate(callmultmoments(c(2,4)),10000,NULL,c(0,0),c(2,1,1,4))
}
