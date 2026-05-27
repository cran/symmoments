#' Convert a Phylogenetic Tree from a Moment L-Matrix to Matching Form
#'
#' Converts a tree structure represented in a moment format into an \code{ape} 
#' matching format structure.
#'
#' @param L The input tree structure. This can be an \code{L-matrix} object, a 
#'   square \eqn{L} matrix, or an \eqn{L} matrix in reduced upper-triangular (vector) form.
#' @param type A character string, either \code{'square'} or \code{'ut'}. This must be 
#'   specified if \code{L} is a raw matrix or vector rather than a formal \code{L-matrix} object. 
#'   Defaults to \code{NULL}.
#' @param tip.label A character vector containing custom labels for the tips. If \code{NULL} 
#'   (the default), labels fallback to \code{"a"} through \code{"z"} if there are at most 26 tips; 
#'   otherwise, 3-letter combinations of the form \code{"aaa"}, \code{"aab"}, etc., are generated.
#'
#' @return A matching representation of the phylogenetic tree corresponding to the input. 
#' The output list is assigned the class \code{'L-matching'}, which contains 5 components 
#' including the tree in matching format.
#'
#' @details An \code{L-matrix} object is a list containing the following 5 components:
#' \itemize{
#'   \item \code{L}: The L-matrix in full square form.
#'   \item \code{L.ut}: The L-matrix in reduced upper-triangular form.
#'   \item \code{Newick}: The Newick string representation of the tree structure.
#'   \item \code{tip.label}: A character vector of the tip labels.
#'   \item \code{tip.label.n}: An integer specifying the total number of tips.
#' }
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#' \cr\cr
#' \insertRef{Diaconis1998}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{toMoment}}, \code{\link{toNewick}}
#'
#' @examples
#' # Create a Newick character string
#' exam.Newick <- "(((a,b),c),d);"
#' 
#' # Convert to a moment L-matrix
#' exam.moment <- toMoment(exam.Newick)
#' 
#' # Convert to matching format
#' exam.matching <- toMatching(exam.moment)
#'
#' @export
`toMatching` <- 
  function (L, type = NULL, tip.label = NULL) 
  {
    toSquare <- function (L.ut) 
    {
      n <- (-1 + sqrt(1 + 8 * length(L.ut)))/2
      L <- 0 * diag(n)
      start <- 1
      for (irow in 1:n) {
        L[irow, ] <- c(rep(0, (irow - 1)), L.ut[start:(start + 
                                                         n - irow)])
        start <- start + n - irow + 1
      }
      return(L)
    }
    
    if (inherits(L,"L-matrix")) {
      L.matrix <- L$L
      if (!is.null(L$tip.label) & is.null(tip.label)) {
        tip.label <- L$tip.label
      }
    }
    if (!is.null(type)) {
      if (type == "square") {
        L.matrix <- L
      }
      if (type == "ut") {
        L.matrix <- toSquare(L)
      }
    }
    nL <- dim(L.matrix)[1]
    n <- 1 + nL/2
    temp.matching <- matrix(rep(0, 3 * (n - 1)), ncol = 3)
    if (!is.null(tip.label)) {
      temp.tip.label <- tip.label
    }
    if (is.null(tip.label)) {
      alphab <- gsub(" ", "", paste(rep(letters, times = rep(676, 
                                                             26)), paste(rep(letters, times = rep(26, 26)), letters)))
      if (!is.null(tip.label)) {
        temp.tip.label <- tip.label
      }
      if (is.null(tip.label)) {
        if (n <= 26) {
          temp.tip.label <- letters[1:n]
        }
        if (n > 26) {
          temp.tip.label <- alphab[1:n]
        }
      }
    }
    temp.node.label <- temp.tip.label
    nrows <- 0
    for (icol in 2:nL) {
      irow <- (1:icol)[L.matrix[1:icol, icol] == 1]
      if (length(irow) == 1) {
        ccol <- temp.node.label[icol]
        crow <- temp.node.label[irow]
        if (icol <= n) {
          couple <- gsub(" ", "", paste("(", crow, ",", 
                                        ccol, ")"), fixed = TRUE)
        }
        if (icol > n) {
          couple <- gsub(" ", "", paste("(", ccol, ",", 
                                        crow, ")"), fixed = TRUE)
        }
        nrows <- nrows + 1
        temp.node.label <- c(temp.node.label, couple)
        temp.matching[nrows, ] <- c(irow, icol, length(temp.node.label))
      }
    }
    temp.node.label <- gsub("!", "", temp.node.label, fixed = TRUE)
    match.obj <- list(matching = temp.matching, tip.label = temp.tip.label, 
                      node.label = temp.node.label)
    class(match.obj) <- "matching"
    return(match.obj)
  }