#' Recursive Function to Compute a Multivariate Moment
#'
#' Called by \code{\link{callmultmoments}} to compute the representation of a 
#' multivariate normal moment using a recursive algorithm.
#'
#' @param moment A vector \code{c(k1, ..., kn)} specifying the moment 
#'   \eqn{X_1^{k_1} \cdots X_n^{k_n}}.
#' @param current.matrix Upper-triangular integer matrix under consideration 
#'   in the recursion.
#' @param current.cell Cell in the current matrix under consideration in 
#'   the recursion.
#' @param moment.rep Current set of representations; \code{multmoments} adds 
#'   each satisfying matrix to \code{moment.rep}.
#' @param row_col Matrix giving rows and columns for a square matrix for 
#'   each cell.
#'
#' @return The moment representation, \code{moment.rep}, augmented with 
#' additional representations.
#'
#' @details Each row of the representation gives the exponents for a single 
#' product of covariance terms. For example, \code{(1, 2, 0)} represents 
#' \eqn{S_{11}^1 S_{12}^2 S_{22}^0}, where the \eqn{S_{ij}} are the covariances.
#' \cr\cr
#' This function would normally only be called by \code{\link{callmultmoments}}.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#' 
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{callmultmoments}}
#'
#' @keywords internal
#' @export

`multmoments` <- 
  function (moment, current.matrix, current.cell, moment.rep, row_col) 
  {
    # A recursive function to compute the representation of a multivariate moment
    #    using upper - triangular matrices
    #
    # moment: a vector of integers representing the moment
    #     eg, c(3, 2, 4) for a 3-dimensional normal vector
    #     corresponding to the moment E[(X1^3)(X2^2)(X3^4)]
    
    # current.matrix: input/output matrix under consideration in recursion
    #     this is an upper-triangular integer matrix (l(i, j))
    #     organized by row, with   lmom * (lmom + 1) / 2   total elements  
    #
    # current.cell: input/output cell in current matrix under consideration in recursion
    
    # moment.rep: input/output current set of representations
    #     function adds each satisfying matrix to moment.rep
    
    # row_col: 2x(lmom * (lmom + 1) / 2 matrix giving rows and columns for square matrix
    #     for each cell so that they don't have to be calculated each time
    
    #  algorithm:
    #  loop through cell values from 0 to min(moments[row] - rowsum, moments[col] - colsum)
    #    if the this new matrix satisfies moment criterion, 
    #       add to moment.rep and return
    #    if the current matrix is too great in any dimension, 
    #       return
    #    if the current matrix is < moment in any dimension, and
    #       at most moment(i) for all other indexs i, continue
    
    summomentmatrix <- function (moment.matrix) 
    {
      #  compute the row/col sums of moment.matrix
      #  uses Matrix package
      #  construct a square upper diagonal matrix from moment.matrix
      
      #  length of moment based on representation:
      lmom <- (sqrt(8 * length(moment.matrix) + 1) - 1) / 2 
      #  make a square matrix for use with row and column sum functions
      tempmatrix <- matrix(c(rep(0, lmom^2)), nrow=lmom)
      tempmatrix[1, ] <- moment.matrix[1:lmom]
      endrow <- lmom
      if (lmom > 1)
      {
        for (irow in 2:lmom)
        {startrow <- endrow + 1
        endrow <- endrow + (lmom - irow + 1)
        tempmatrix[irow, ] <- cbind(c(rep(0, (irow - 1)), moment.matrix[startrow:endrow])) }
        
        summomentmatrix <- colSums(tempmatrix) + rowSums(tempmatrix) 
      }
      if (lmom == 1){summomentmatrix <- 2 * sum(moment.matrix)}     
      
      return(summomentmatrix)}
    
    
    lmom <- length(moment)
    totcells <- lmom * (lmom + 1) / 2   # total cells in a moment representation
    
    thisrow <- row_col[1, current.cell]
    thiscol <- row_col[2, current.cell]
    
    moment.row <- moment[thisrow]
    moment.col <- moment[thiscol]
    rowcells <- row_col[1, row_col[1, ] == thisrow]
    colcells <- row_col[2, row_col[2, ] == thiscol]
    
    #  determine sums for row and columm
    #  maxvalue is the largest value that can be added to a
    #    a row/column sum and still be no more than criterion
    #  maxvalue will be the minimum of the moments minus these sums
    
    rowsum <- summomentmatrix(current.matrix)[thisrow]
    colsum <- summomentmatrix(current.matrix)[thiscol] 
    
    maxvalue <- min(moment[thisrow] - rowsum, moment[thiscol] - colsum)
    
    for (ivalue in 0:maxvalue)
    {         # for:loop through all possible values for cell
      
      current.matrix[current.cell] <- ivalue
      current.sum <- summomentmatrix(current.matrix)
      # determine current sum for this matrix
      
      #  if matrix fulfills criterion, add to reps and return
      if (sum(moment == current.sum) == lmom)
      {moment.rep <- rbind(moment.rep, current.matrix)
      return(moment.rep)}
      
      #   if the sum is too large, return because any other sum will also be too large
      if (sum(current.sum > moment) > 0){return(moment.rep)}
      
      #   if at least one term in current sum is smaller than moment, 
      #       go down one cell unless this is the last cell
      
      if ((sum(current.sum == moment) < lmom) & (current.cell != totcells))     
      {cc <- current.cell + 1
      # recursive step:
      moment.rep <- multmoments(moment, current.matrix, cc, moment.rep, row_col)}
      
    }        # end of for 
    
    #   note: the three conditions above are exhaustive except for last cell 
    #      In that case there is nothing to return
    return(moment.rep)
  }