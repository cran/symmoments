#' Compute an Unsorted Central Moment Object from a Sorted Object
#'
#' Produces an unsorted central moment object from a sorted object of class 
#' \code{'moment'}.
#'
#' @param moment The unsorted target moment to obtain, specified in vector 
#'   form (e.g., \code{c(3, 1, 2)}).
#' @param sorted.moment A sorted object of class \code{'moment'} to use as the 
#'   base for creating the unsorted moment.
#'
#' @return An object of class \code{'moment'}, which is a list containing the 
#' following three components:
#' \item{moment}{The input unsorted moment vector.}
#' \item{representation}{A matrix containing the representation in terms of 
#'   upper-triangular matrices, rearranged to match the target unsorted order.}
#' \item{coefficients}{A numeric vector of the coefficients corresponding to 
#'   the rows of the representation matrix.}
#'
#' @details Unsorted moments are those whose exponents are not in sorted 
#' numerical order (e.g., \code{m312} vs \code{m123}). The unsorted moment's 
#' representation is calculated by rearranging the rows and columns of the 
#' sorted moment's matrices successively.
#'
#' @references 
#' #' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{multmoments}}, \code{\link{callmultmoments}}
#'
#' @examples
#' # Obtain unsorted moment m312 from sorted base m123
#' tounsorted(c(3, 1, 2), callmultmoments(c(1, 2, 3)))
#'
#' @export
`tounsorted` <- 
  function (moment,sorted.moment) 
  {
    # converts a sorted moment to a specified unsorted moment
    # eg,   m(2,3,5) ->  m(5,2,3)
    # each row is sorted separately
    # the rows may be ordered differently
    
    # moment: noncanonical moment to obtain 
    #         moment is in vector form, eg, c(3,1,2)
    
    # sorted.moment: canonical moment of class "moment"
    #         this moment is monotone in its powers;
    
    
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
    
    
    
    nestedreps <- function (input.vector, inner.rep, outer.rep) 
    {
      # replicates input.vector, first by inner.rep, then by outer.rep
      
      # input.vector: vector to replicate
      # inner.rep:    count of replicates of elements
      # outer.rep:    count of replicates of resulting vector
      
      temp.vector <- NULL
      for (ivec in 1:length(input.vector))
      {temp.vector <- c(temp.vector, rep(input.vector[ivec], inner.rep))}
      
      total.vector <- rep(temp.vector, outer.rep) 
      
      return(total.vector)}
    
    
    output.moment <- sorted.moment
    output.moment$moment <- moment
    lmom <- length(moment)
    r <- order(moment)
    lrep <- lmom*(lmom+1)/2
    representation <- sorted.moment$representation
    
    nrep <- dim(representation)[1]
    overallcoeff <- ((1/2)^(sum(moment)/2))*prod(combinat::fact(moment))/combinat::fact(sum(moment)/2)
    
    m <- matrix(rep(1, lmom^2), nrow=lmom) 
    
    limits <- c((lmom:1)%*%(m * !(lower.tri(m))))
    # sum of row lengths: lmom, lmom + lmom - 1, ... , for use in row_col
    row_col <- matrix(rep(0, (2 * lrep)), nrow=2)
    for (icell in (1:lrep))
    {
      row_col[1, icell] <- min((1:lmom)[icell<=limits])
      if (row_col[1, icell] == 1){row_col[2, icell] <- icell}
      if (row_col[1, icell]>1){row_col[2, icell] <- icell - limits[row_col[1, icell] - 1] + 
        row_col[1, icell] - 1 }
    }
    #  2x(nm * (nm + 1) / 2 matrix giving rows and columns for each cell
    #    compute here so that they don't have to be calculated each time 
    
    
    for (irep in 1:nrep)
    {
      
      utri <- toSquare(representation[irep,])
      noncanonical.matrix <- 0*utri
      
      for (irow in 1:lmom)
      {
        jrow <- r[irow]
        for (icol in irow:lmom)
        {
          jcol <- r[icol]
          if (jrow <= jcol)
          {noncanonical.matrix[jrow,jcol] <- utri[irow,icol]}
          if (jrow > jcol)
          {noncanonical.matrix[jcol,jrow] <- utri[irow,icol]}
        }
      }
      thisrep <- t(noncanonical.matrix)[t(!lower.tri(noncanonical.matrix))]
      
      
      output.moment$representation[irep,] <- thisrep
      
      
      #  determine the coefficient for each term based on switching equivalent terms
      #          this is taken from callmultmoments
      
      #  "base" gives the number of switches that can be made to each element of the l-matrix
      #  diagonal elements are not switchable, but are included to allow subtraction below
      
      base <- c(rep(1, lmom * (lmom + 1) / 2))
      base[1] <- 1   # first diagonal element  -  not switchable
      
      totreps <- 1   # total number of transpostions
      # if there is only one element, it must be the diagonal, so is not switchable  -  skip
      if (lmom > 1){
        base[1] <- 1   # first diagonal element
        for (cell in 2:length(base))
        {
          icol = row_col[1, cell]  #  determine if diagonal element
          irow = row_col[2, cell]
          if (irow == icol){base[cell] <- 1}  # diagonal  -  not switchable
          if (icol != irow)
          {totreps <- totreps * (1 + thisrep[cell])
          base[cell] <- 1 + thisrep[cell]} 
        }  #  done with computing base and total transpositions (totreps)
      }
      
      mcoeff <- 1  #  sum of multinomial coefficients
      if (totreps > 1){ 
        
        #  baserep will represent the lower diagonal (including diagonal)
        #  of the augmented matrices
        baserep <- matrix(rep(0, totreps * length(base)), nrow=totreps)
        basegt1 <- base[base>1]
        nbase <- 0
        for (ibase in 1:length(base))
        {if (base[ibase] > 1)
        {nbase = nbase + 1
        if (nbase == 1){baserep[, ibase] <- nestedreps(c(0:(basegt1[nbase] - 1)), 1, totreps / prod(basegt1[1:nbase])) }    
        if (nbase > 1) {baserep[, ibase] <- nestedreps(c(0:(basegt1[nbase] - 1)), prod(basegt1[1:(nbase - 1)]), totreps / prod(basegt1[1:nbase])) }
        }
        }
        
        #  now go through each transposition
        if ( !is.na(totreps) & totreps != 1)
        {
          mcoeff <- 0
          for (jrep in (1:totreps)) # check each transposition
          {newrep <- baserep[jrep, ]       #  added lower diagonal elements
          addrep <- sort(newrep, decreasing=TRUE)[1:(lmom * (lmom - 1) / 2)] 
          fulnrep <- c((thisrep - newrep), addrep)
          thiscoeff <- ((length(fulnrep))^sum(fulnrep)) * dmultinom(x=fulnrep, prob=rep(1.0, length(fulnrep)))
          mcoeff <- mcoeff + thiscoeff
          #  the multinomial coefficient is obtained from the multinomial distribution
          #  multiply by an appropriate power to get rid of probability
          }  
        }
      }  
      if (is.na(totreps)){mcoeff <- 1}
      if (totreps == 1)
      {mcoeff <- (length(thisrep))^sum(thisrep) * dmultinom(x=thisrep, prob=rep(1.0, length(thisrep)))}
      
      #  determine full coefficient  -  round because all coefficients should be integers
      #                                (Note - this statement has not been proved)
      
      output.moment$coefficients[irep] <- round(overallcoeff * mcoeff)
      
    }  # end of representations
    
    
    return(output.moment)} 