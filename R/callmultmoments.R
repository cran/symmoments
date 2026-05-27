#' Compute Multivariate Moment Symbolically
#'
#' Computes a multivariate normal moment by initializing variables, calling 
#' \code{multmoments}, and constructing output.
#'
#' @param moment A vector \code{c(k1, ..., kn)} specifying the moment 
#'   \eqn{X_1^{k_1} \cdots X_n^{k_n}}.
#'
#' @return An object of class \code{'moment'}, which is a list with three 
#' components:
#' \itemize{
#'   \item \code{moment}: The input moment vector.
#'   \item \code{representation}: A matrix containing the representation in 
#'     terms of upper-triangular matrices.
#'   \item \code{coefficients}: The coefficients corresponding to the rows of 
#'     the representation.
#' }
#' If the sum of the exponents is odd, returns \code{-1} and prints 
#' "Sum of powers is odd. Moment is 0."
#' \cr\cr
#' If any exponent is negative, returns \code{-2} and prints 
#' "All components of the moment must be non-negative."
#' \cr\cr
#' If any exponent is not an integer, returns \code{-3} and prints 
#' "All components of the moment must be integers."
#'
#' @details Each row of the representation gives the exponents for a single 
#' product of covariance terms. For example, \code{(1, 2, 0)} represents 
#' \eqn{S_{11}^1 S_{12}^2 S_{22}^0}, where the \eqn{S_{ij}} are the covariances.
#' \cr\cr
#' The full moment is the sum of these terms multiplied by their respective 
#' coefficients. If the sum of the exponents is odd, the moment is 0.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{multmoments}}, and the methods \code{\link{toLatex}}, 
#' \code{\link{evaluate}}, and \code{\link{simulate}} in the \code{symmoments} 
#' package.
#'
#' @examples
#' # Compute the moment for the 4-dimensional moment c(1,2,3,4):
#' m.1234 <- callmultmoments(c(1, 2, 3, 4))
#' 
#' @export

`callmultmoments` <- 
  function (moment) 
  {
    #  function to compute the representation of a multivariate moment
    # 
    #  moment: a vector of integers representing the moment
    #          eg, c(3, 2, 4) for a 3-dimensional normal vector
    #          corresponding to the moment (E[X1^3)(X2^2)(X3^4)]
    #  sum of exponents must be even; otherwise moment is 0
    
    #  returns a list of 3 components:
    #    1: $moment  -  the input moment vector 
    #    2: $representation  -  a matrix containing the representation 
    #       in terms of upper-triangular matrices
    #    3: $coefficients  -  the coefficients corresponding to the rows of the representation
    #  if sum odd, returns  -1 and prints "Sum of powers is odd. Moment is 0."
    #  if any component is negative, returns  -2 and prints "All components of the moment must be non-negative."
    #  if any component is not an integer, returns  -3 and prints "All components of the moment must be integers."
    
    subblank <- function (inputstring) 
    {
      # remove blanks from a character string
      outputstring <- sub(" ", "", inputstring)
      if (outputstring == inputstring) {
        return(outputstring)
      }
      if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
      }
      return(outputstring)
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
    
    
    mrepnames <- function (ndim) 
    {
      # get colnames for a representation
      combs <- expand.grid(1:ndim, 1:ndim)
      char.combs <- paste("S(", combs[, 2], ", ", combs[, 1], ")")
      for (elem in 1:(ndim^2))
      {char.combs[elem] <- subblank(char.combs[elem])}
      m <- sort(matrix((1:ndim^2), nrow=ndim, byrow=TRUE)[ !lower.tri(matrix((1:ndim^2), nrow=ndim, byrow=TRUE))])
      return(char.combs[m])
    }
    
    
    lmom <- length(moment)
    lrep <- (lmom^2 + lmom) / 2  # length of representation using upper-triangular matrices
    moment.rep <- matrix(rep(0, lrep), nrow=1)   # initial null representation to be augmented
    if (sum(trunc(moment) == moment) < length(moment)) 
    {print("All components of the moment must be integers.") 
      return( -3)}
    if (sum(moment<0) > 0)
    {print("All components of the moment must be non-negative.")
      return( -2)}
    if (trunc(sum(moment) / 2) != sum(moment) / 2)
    {print("Sum of powers is odd. Moment is 0.")
      return( -1)}
    
    icells <- matrix(rep(0, lmom^2), nrow=lmom)
    m <- matrix(rep(1, lmom^2), nrow=lmom) 
    limits <- c((lmom:1)%*%(m * !(lower.tri(m))))
    # sum of row lengths: lmom, lmom + lmom - 1, ... , for use in row_col
    
    row_col <- matrix(rep(0, (2 * lrep)), nrow=2)
    #  2x(nm * (nm + 1) / 2 matrix giving rows and columns for each cell
    #       so that they don't have to be calculated each time in multmoment
    
    for (icell in (1:lrep))
    {
      row_col[1, icell] <- min((1:lmom)[icell<=limits])
      if (row_col[1, icell] == 1){row_col[2, icell] <- icell}
      if (row_col[1, icell]>1){row_col[2, icell] <- icell - limits[row_col[1, icell] - 1] + 
        row_col[1, icell] - 1 }
    }
    
    # initial current.matrix and current.cell
    current.matrix <- c(rep(0, lrep))  
    current.cell <- 1
    
    # call recursive function to determine upper-triangular representations
    
    moment.rep <- multmoments(moment, current.matrix, current.cell, moment.rep, row_col)
    if (dim(moment.rep)[1] == 2)    # get rid of initial 0 matrix representation
    { moment.rep <- matrix(moment.rep[2, ], nrow=1) }
    if (dim(moment.rep)[1] > 2)
    {moment.rep <- moment.rep[2:dim(moment.rep)[1], ]}
    rownames(moment.rep) <- 1:(dim(moment.rep)[1])
    
    ##################################################################
    # now determine coefficients for upper-triangular representations
    
    l.representation <- moment.rep
    lmom <- length(moment)
    nrep = dim(l.representation)[1]  
    totlength <- lmom^2
    rep.coefficients <- c(1:nrep)  # coefficients corresponding to nrep representations
    
    #  multiplier for all terms
    overallcoeff <- ((1 / 2)^(sum(moment) / 2)) * prod(factorial(moment)) / factorial((sum(moment) / 2)) 
    
    for (irep in 1:(dim(l.representation)[1]))
    {
      #  loop through all matrices
      thisrep <- l.representation[irep, ] 
      
      #  determine the coefficient for each term based on switching equivalent terms
      
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
      rep.coefficients[irep] <- round(overallcoeff * mcoeff)
      
      cell <- 0
      for (irow in (1:length(moment)))
      {
        for (icol in (irow:length(moment)))
        {cell <- cell + 1
        #   exponent of term
        expo <- l.representation[irep, cell]
        
        }
      }
    }
    output <- list(moment, moment.rep, rep.coefficients)
    names(output) <- c('moment', 'representation', 'coefficients')
    colnames(output$representation) <- mrepnames(length(moment)) 
    names(output$coefficients) <- paste("rep", (1:length(output$coefficients)))
    class(output) <- "moment" 
    return(output)}