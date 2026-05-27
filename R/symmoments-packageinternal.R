#' Symbolically compute and numerically evaluate multivariate central moments
#'
#' Symbolically computes and numerically evaluates multivariate normal moments
#'  \eqn{E[X_1^{k_1} \cdots X_n^{k_n}]}, where \eqn{(X_1,\ldots,X_n) \sim N(\mu, \Sigma)},
#' in terms of mu and S elements.
#'
#' Produces Latex code for the moment.
#'
#' Computes numerical moments at specified means and covariance matrices.
#'
#' Also converts between moment L-matrices, phylo objects, and matching objects.
#'
#' A representation of a central moment of the multivariate normal distribution,
#' given by a positive integer vector c(k1,k2,...,kn), is obtained from the
#' function [callmultmoments()]. This function initializes variables and calls
#' the function [multmoments()] which determines a representation of a
#' multivariate moment using a recursive algorithm. The representation is given
#' class 'moment'.
#'
#' The [print()] method prints the representation of a multivariate moment.
#'
#' The [toLatex()] method uses the output of [callmultmoments()] to determine
#' the LaTeX code for the moment sorted lexicographically.
#'
#' The generic [evaluate()] method uses the output of [callmultmoments()] to
#' determine the value of the moment for a specified covariance matrix.
#'
#' The [simulate()] method is used to approximate a (possibly non-central)
#' moment using Monte Carlo integration.
#'
#' The [toLatex_noncentral()] function computes the Latex representations of a
#' non-central moment.
#'
#' The [evaluate_noncentral()] function computes the value of a non-central
#' moment.
#'
#' The [evaluate_expected.polynomial()] function evaluates the expected value
#' of a multivariate polynomial defined by a list, multipol object, or mpoly
#' object.
#'
#' The [convert.multipol()] function converts between multipol objects and
#' multivariate polynomials defined by lists.
#'
#' The [convert.mpoly()] function converts between mpoly objects and
#' multivariate polynomials defined by lists.
#'
#' The [tounsorted()] function converts a sorted moment (e.g. m123) to an
#' unsorted moment (e.g. m312).
#'
#' The [make.all.moments()] function computes all moments up to a specified
#' size and places them in the symmoments environment.
#'
#' The [integrate.polynomial()] function integrates a multivariate polynomial
#' against the normal distribution using ordinary integration.
#'
#' The functions [toMoment()], [toNewick()], and [toMatching()] convert among
#' moment L-matrices, Newick trees, and `ape` matching objects.
#'
#' @note The mvtnorm package must be loaded for the simulate method.
#'   The cubature package must be loaded for the integrate.polynomial function.
#'   The combinat package must be loaded for the toMoment function.
#'
#' @author 
#' Kem Phillips \email{kemphillips@@comcast.net}
#' \cr\cr
#' Maintainer: Florian Schuberth \email{f.schuberth@utwente.nl}
#'
#' @references
#' \insertRef{Phillips2010}{symmoments}
#'
#' @examples
#' # Compute the moment for the 4-dimensional moment c(1,2,3,4):
#' callmultmoments(c(1,2,3,4))
#'
#' # Print the representation of the 4-dimensional moment c(1,2,3,4):
#' print(callmultmoments(c(1,2,3,4)))
#'
#' # Compute the LaTeX representation of the central moment c(1,2,3,4):
#' toLatex(callmultmoments(c(1,2,3,4)))
#'
#' # evaluate the moment c(1,2,3,4) at the following variance-covariance matrix
#' #  4 2 1 1
#' #  2 3 1 1
#' #  1 1 2 1
#' evaluate(callmultmoments(c(1,2,3,4)), c(4,2,1,1,3,1,1,2,1,2))
#'
#' # Using 10000 samples, estimate the central moment for c(2,4) (not run)
#' # at the covariance matrix:
#' #  2 1
#' #  1 4
#' # and mean (0,0)
#' \dontrun{
#' library(mvtnorm)
#' simulate(callmultmoments(c(2,4)), 10000, NULL, c(0,0), c(2,1,1,4))
#' }
#'
#' \dontrun{
#' # Compute Latex representation of a non-central moment
#' as.matrix(toLatex_noncentral(c(1,3)))
#'
#' # Create all 2-dimensional moment objects with exponents up to 3
#' symmoments <- new.env()
#' make.all.moments(c(3,3))
#'
#' # Evaluate a non-central moment (requires moments of order up to c(1,3)
#' # to exist in environment symmoments)
#' evaluate_noncentral(c(1,3), c(1,2), c(1,0,1))
#' }
#'
#' # Create an mpoly object
#' library(mpoly)
#' library(multipol)
#' t0 <- mpoly(list(c(coef=3, x1=2), c(coef=2, x1=1, x2=3),
#'                  c(coef=-4, z=2), c(coef=1, x1=1, x2=2, z=1)))
#'
#' # Convert an mpoly object to a moment object
#' t1 <<- convert.mpoly(t0)
#'
#' # Convert a moment object to a multipol object
#' t2 <<- convert.multipol(t1)
#'
#' # Convert from multipol back to mpoly through moment
#' mpoly(convert.mpoly(convert.multipol(t2)))
#'
#' \dontrun{
#' # Evaluate the expected value of a multivariate polynomial
#' # (required moments must exist in environment symmoments)
#' evaluate_expected.polynomial(t0, c(1,2,3), c(1,0,0,1,0,1))
#' }
#'
#' # Create a Newick representation of a tree
#' exam.Newick   <- "(((a,b),c),d);"
#'
#' # Convert to phylo format
#' library(ape)
#' exam.phylo    <- read.tree(text = exam.Newick)
#'
#' # Convert to matching format
#' exam.matching <- as.matching(exam.phylo)
#'
#' # Convert to L-matrix format
#' exam.L.matrix <- toMoment(exam.matching)
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats dmultinom
#' @importFrom multipol as.multipol
#'
#' @keywords internal
"_PACKAGE"