#' Convert a Tree from Newick or Matching to Moment Format
#'
#' Converts a phylogenetic tree from a Newick character string or an \code{ape} 
#' matching matrix into a moment L-matrix object.
#'
#' @param inputobject A tree structure represented as a Newick format character 
#'   string, or a \code{matching} object as defined in the \pkg{ape} package.
#' @param tip.label A character vector specifying rearranged labels for the 
#'   tips. If provided, these must be the original tip labels. Defaults to \code{NULL}.
#'
#' @return A moment L-matrix object corresponding to the input phylogenetic 
#' tree object.
#'
#' @details The returned L-matrix class object consists of 5 internal components:
#' \itemize{
#'   \item \code{L}: The L-matrix represented in full square form.
#'   \item \code{L.ut}: The L-matrix represented in upper-triangular form.
#'   \item \code{Newick}: The Newick string representation of the tree structure.
#'   \item \code{tip.label}: A character vector holding the labels of the tips.
#'   \item \code{tip.label.n}: An integer specifying the total number of tips.
#' }
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#' \cr\cr
#' Felsenstein, J. (1990). The Newick tree format. 
#' \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
#' \cr\cr
#' \insertRef{Diaconis1998}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{toNewick}}, \code{\link{toMatching}}
#'
#' @examples
#' # Create a Newick character string
#' exam.Newick <- "(((a,b),c),d);"
#' 
#' # Convert to a moment L-matrix
#' exam.moment <- toMoment(exam.Newick)
#' 
#' # Convert to a matching object
#' exam.matching <- toMatching(exam.moment)
#' 
#' # Convert back to a moment object
#' backto.moment <- toMoment(exam.matching)
#'
#' @export
`toMoment` <- 
  function (inputobject, tip.label = NULL) 
  {
    if (!inherits(inputobject,"matching") & !is.matrix(inputobject)) {
      if (!inherits(inputobject,"L-Newick")) {
        temp <- inputobject
        Newick.out <- inputobject
      }
      if (inherits(inputobject,"L-Newick")) {
        temp <- inputobject$Newick
        Newick.out <- inputobject$Newick
      }
      temp <- gsub(";", " ", temp, fixed = TRUE)
      temp <- gsub("(", " ", temp, fixed = TRUE)
      temp <- gsub(")", " ", temp, fixed = TRUE)
      temp <- gsub(",", " ", temp, fixed = TRUE)
      temp <- gsub("  ", " ", temp, fixed = TRUE)
      tips <- strsplit(temp, " +")[[1]]
      tips <- tips[tips != ""]
      n <- length(tips)
      neworder <- 1:n
      if (!is.null(tip.label)) {
        for (itips in 1:n) {
          neworder[itips] <- grep(tip.label[itips], tips)
        }
      }
      tips <- tips[neworder]
      np <- n
      temp.tips <- rep(" ", n)
      temp.Newick <- gsub(";", " ", inputobject, fixed = TRUE)
      for (itips in 1:n) {
        temp.tips[itips] <- gsub(" ", "", paste(";", as.character(itips), 
                                                ";"), fixed = TRUE)
        temp.Newick <- gsub(tips[itips], temp.tips[itips], 
                            temp.Newick, fixed = TRUE)
      }
      n.temp.tips <- n
      temp.tip.names <- temp.tips
      for (itips in 1:n) {
        temp.tips[itips] <- gsub(" ", "", temp.tips[itips], 
                                 fixed = TRUE)
      }
      couples.n <- 1
      while (couples.n > 0) {
        couples.n <- 0
        for (cspec in 1:np) {
          for (rspec in 1:np) {
            couple <- gsub(" ", "", paste("(", temp.tips[cspec], 
                                          ",", temp.tips[rspec], ")"), fixed = TRUE)
            if (length(grep(couple, temp.Newick)) > 0) {
              couples.n <- couples.n + 1
              n.temp.tips <- n.temp.tips + 1
              temp.tip.names <- c(temp.tip.names, couple)
              temp.tips <- c(temp.tips, gsub(" ", "", paste(";", 
                                                            as.character(n.temp.tips), ";"), fixed = TRUE))
              temp.Newick <- gsub(couple, gsub(" ", "", 
                                               paste(";", as.character(n.temp.tips), ";"), 
                                               fixed = TRUE), temp.Newick, fixed = TRUE)
            }
          }
        }
        np <- np + couples.n
      }
      temp.tip.names.n <- length(temp.tip.names)
      L <- 0 * diag(2 * (n - 1))
      for (icouple in (n + 1):temp.tip.names.n) {
        temp <- gsub("(", " ", temp.tip.names[icouple], fixed = TRUE)
        temp <- gsub(")", " ", temp, fixed = TRUE)
        temp <- gsub(",", " ", temp, fixed = TRUE)
        temp <- gsub(";", " ", temp, fixed = TRUE)
        temp <- gsub("  ", " ", temp, fixed = TRUE)
        ind <- strsplit(temp, " +")[[1]]
        ind <- ind[ind != ""]
        indn <- as.integer(ind)
        indn <- ind[order(indn)]
        irow <- as.integer(indn[1])
        icol <- as.integer(indn[2])
        L[irow, icol] <- 1
      }
      L.ut <- NULL
      for (irow in 1:(2 * (n - 1))) {
        L.ut <- c(L.ut, L[irow, irow:(2 * (n - 1))])
      }
      SMNM <- list(L = L, L.ut = L.ut, Newick = Newick.out, 
                   tip.label = tips, tip.label.n = n)
      class(SMNM) <- "L-matrix"
      return(SMNM)
    }
    if (inherits(inputobject,"matching") | is.matrix(inputobject)) {
      temp.tip.label <- NULL
      if (inherits(inputobject,"matching")) {
        temp.matching <- inputobject$matching[, c(1, 2)]
        if (!is.null(inputobject$tip.label)) {
          temp.tip.label <- inputobject$tip.label
        }
      }
      if (is.matrix(inputobject)) {
        temp.matching <- inputobject[, c(1, 2)]
      }
      nL <- max(temp.matching)
      n <- 1 + nL/2
      L <- 0 * diag(nL)
      if (is.matrix(temp.matching)) {
        L[temp.matching] <- 1
      }
      if (!is.matrix(temp.matching)) {
        L[matrix(temp.matching, nrow = 1)] <- 1
      }
      L.ut <- NULL
      for (irow in 1:(2 * (n - 1))) {
        L.ut <- c(L.ut, L[irow, irow:(2 * (n - 1))])
      }
      if (!is.null(temp.tip.label)) {
        Newick <- toNewick(L, type = "square", tip.label = temp.tip.label)
      }
      if (is.null(temp.tip.label)) {
        Newick <- toNewick(L, type = "square")
      }
      if (is.null(tip.label)) {
        alphab <- gsub(" ", "", paste(rep(letters, times = rep(676, 
                                                               26)), paste(rep(letters, times = rep(26, 26)), 
                                                                           letters)))
        if (n <= 26) {
          tips <- letters[1:n]
        }
        if (n > 26) {
          tips <- alphab[1:n]
        }
      }
      if (!is.null(tip.label)) {
        tips <- tip.label
      }
      if (is.null(tip.label) & inherits(inputobject,"matching")){
        if (!is.null(inputobject$tip.label)) {
          tips <- inputobject$tip.label
        }
      }
      SMNM <- list(L = L, L.ut = L.ut, Newick = Newick$Newick, 
                   tip.label = tips, tip.label.n = n)
      class(SMNM) <- "L-matrix"
      return(SMNM)
    }
  }