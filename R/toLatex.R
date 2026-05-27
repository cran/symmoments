#' LaTeX a Multivariate Moment
#'
#' Computes a LaTeX representation sorted lexicographically of an object of 
#' class \code{'moment'}.
#'
#' @param object An object of class \code{'moment'}, usually the output of 
#'   \code{\link{callmultmoments}}.
#' @param ... Included only for consistency with the generic function.
#'
#' @return A character vector giving the LaTeX code for the symbolic moment.
#'
#' @details The first element of the result is the moment expressed as an 
#' expected value (\code{E[...] =}). The remaining lines are the LaTeX 
#' representation broken at appropriate intervals for printing. (Individual 
#' terms for high dimensions will still overrun a printed line.)
#' \cr\cr
#' Double backslashes (\code{\\\\}) are inserted where LaTeX requires a 
#' backslash. These can be reset to single backslashes by writing the output 
#' to a file using the standard R function \code{\link[base]{writeLines}}.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{callmultmoments}}, and the \code{\link{evaluate}} method.
#'
#' @examples
#' toLatex(callmultmoments(c(1, 2, 3)))
#'
#' @importFrom utils toLatex
#' @method toLatex moment
#' @export
`toLatex.moment` <- 
  function (object, ...) 
  {
    #  build latex code for the l-matrix representation of the moment
    #  object is the representation of the l-matrices for moment
    #  each row is such an l-matrix
    
    # object: list from callmultmoment (class moment)
    #     with first component the moment itself, 
    #     the second component the set of upper-triangular
    #          representations of the moment, 
    #     and third component, their correpsonding coefficients
    
    #  note that Latex backslashes are doubled to allow writing to file 
    #      with writeLines
    
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
    sortmatrix <- function (xmat) 
    {
      # sort matrix xmat by successive columns
      # output is the index vector, not the sorted matrix
      
      ncol <- dim(xmat)[2]
      
      sortstring <- "order("
      for (icol in 1:ncol)
      {
        sortstring <- paste(sortstring, "xmat[, ", icol, "]")
        if (icol != ncol){sortstring <- paste(sortstring, ", ")}
      }
      sortstring <- paste(sortstring, ", decreasing=TRUE)")
      matindex <- eval(parse(text=sortstring))
      
      return(matindex)}
    
    
    moment.fullrep <- object
    
    maxchars.latex <- 300
    # maximum latex characters to print on one line
    # this includes formatting characters (subscript notation, etc)
    
    #  extract the components for convenience
    moment <- moment.fullrep[[1]]
    moment.rep <- moment.fullrep[[2]][sortmatrix(moment.fullrep[[2]]), ]
    #       latex representation will be in sorted order
    coefficients.rep <- moment.fullrep[[3]][sortmatrix(moment.fullrep[[2]])]
    
    
    doubquote <- subblank(paste("\\", "\\"))
    
    if ( !is.matrix(moment.rep)){moment.rep <- matrix(moment.rep, nrow=1)}
    numrep <- dim(moment.rep)[1]
    
    latex.moment <- rep(" ", numrep + 1)  
    
    #   write the left hand side, ie, E[X1 ... Xn] =
    
    latex.moment[1] <- "E["
    
    for (imoment in(1:(length(moment))) )
    {latex.moment[1] <- paste(latex.moment[1], 
                              "X_{", imoment, "}^{",  moment[imoment],  "}")  }
    
    latex.moment[1] <- paste(latex.moment[1], "] =", doubquote)
    latex.moment[1] <- subblank(latex.moment[1])
    
    #  write the right hand side, ie, the set of terms
    
    if (sum(moment==0) == length(moment))
    {latex.moment[2] <- "1"
    return(latex.moment)}
    
    totchars <- 0        # used with totchars.latex
    for (irep in (1:numrep))
    {
      mcoeff <- coefficients.rep[irep]
      thisrep <- moment.rep[irep, ] 
      
      if (mcoeff != 1){latex.moment[irep + 1] <- as.character(mcoeff)}
      #              omit coefficient if it is 1
      cell <- 0
      for (irow in (1:length(moment)))
      {
        for (icol in (irow:length(moment)))
        {cell <- cell + 1
        #   exponent of term, that is, the "l" value
        exponent <- moment.rep[irep, cell]
        
        #   if exponent is 1, omit it as obvious
        if (exponent == 1){latex.moment[irep + 1] <- paste(latex.moment[irep + 1], 
                                                           "\\sigma_{ ",  irow,  ", ",  icol,  "}")}
        
        if (exponent > 1){latex.moment[irep + 1] <- paste(latex.moment[irep + 1], 
                                                          "\\sigma_{ ",  irow,  ", ",  icol, "}^{",  exponent, "} ")}
        }  # end of for
      }  # end of for
      totchars <- totchars + nchar(latex.moment[irep + 1])
      if (irep < numrep){latex.moment[irep + 1] <- paste(latex.moment[irep + 1], " + ")}
      if (totchars > maxchars.latex) 
      {totchars <- 0
      latex.moment[irep + 1] <- paste(latex.moment[irep + 1], doubquote)}
      latex.moment[irep + 1] <-  subblank(latex.moment[irep + 1])
    }  # end of for
    
    return(latex.moment)}




#' Compute a LaTeX Expression for a Non-Central Moment
#'
#' Computes a LaTeX expression for a non-central multivariate normal moment.
#'
#' @param moment A vector \code{c(k1, ..., kn)} specifying the moment 
#'   \eqn{X_1^{k_1} \cdots X_n^{k_n}}.
#' @param envir A character string specifying the environment that contains 
#'   the required central moments. Defaults to \code{'symmoments'}.
#'
#' @return A character string giving the LaTeX representation of the non-central 
#' moment where \eqn{X} follows a multivariate normal distribution.
#'
#' @details All required central moment objects must already exist in the 
#' specified environment (the default is \code{'symmoments'}). However, if only 
#' the sorted version of an unsorted moment exists in that environment, the 
#' \code{\link{tounsorted}} function will automatically be called to transform 
#' and obtain it.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{make.all.moments}}, \code{\link{tounsorted}}, 
#'   \code{\link{callmultmoments}}, and the method \code{\link{toLatex}}.
#'
#' @examples
#' \dontrun{
#' # Compute the LaTeX representation of the 2-dimensional non-central moment c(1,3).
#' # Note: This requires that all central moments up to c(1,3) have already been 
#' # generated in the symmoments environment using make.all.moments.
#' toLatex_noncentral(c(1, 3))
#' }
#'
#' @export
`toLatex_noncentral` <-
  function (moment,envir='symmoments') 
  {
    # Compute the Latex representation of a noncentral moment 
    
    
    exists.envir <- FALSE
    if (exists(envir))
    {
      if (inherits(eval(parse(text=envir)),'environment'))
      {exists.envir <- TRUE}
    }
    if (!exists.envir)
    {return(paste('There is no environment named ', envir)) }
    
    
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
    
    make.moment.name <- function (moment) 
    {
      # returns character name of moment
      # moment is the moment as a vector, eg, c(1,2,3)
      
      allchar <- c(0:9,letters,toupper(letters))
      ndim <- length(moment)
      moment.name <- subblank(paste("m",allchar[moment[1]+1]))
      if (ndim > 1)
      {for (idim in 2:ndim)
      {Ak <- allchar[moment[idim]+1]
      moment.name <- subblank(paste(moment.name,Ak))
      }
      }
      
      return(moment.name)
    }
    muproduct.Latex <- function (exponent) 
    {
      # create Latex expression for a mu-product with this exponent
      ndim <- length(exponent)
      term <- " "
      if (exponent[1] > 0)
      {
        if (exponent[1] == 1)
        {term <- "\\mu_{1}"}     
        if (exponent[1] > 1)  
        {term <- paste("\\mu_{1}^{",exponent[1],"}")}
      }
      if (ndim > 1)
      {
        for (idim in 2:ndim)
        {
          if (exponent[idim] == 1)
          {term <- paste(term,"\\mu_{",idim,"}")}
          if (exponent[idim] > 1)
          {term <- paste(term,"\\mu_{",idim,"}^{",exponent[idim],"}")}
        }
      }
      if (sum(exponent) == 0)
      {term <- " "}
      term <- subblank(term)
      
      return(term)
    }
    
    doubquote <- subblank(paste("\\", "\\"))
    fullmoment <- "E["
    for (imoment in(1:(length(moment))) )
    {fullmoment <- paste(fullmoment, 
                         "X_{", imoment, "}^{",  moment[imoment],  "}")  }
    fullmoment <- paste(fullmoment, "\\mid \\mu, \\Sigma] =", doubquote)
    fullmoment <- subblank(fullmoment)
    
    product = prod(moment+1)
    cumproduct = cumprod(moment+1)
    mdim = length(moment)
    thisone <- rep(0,mdim)
    
    # special case: X**0
    if (mdim == 1 & sum(moment) == 0)
      fullmoment <- paste(fullmoment," 1 \\\\")
    
    if (mdim > 1 | sum(moment) > 0)
    { 
      
      for (mcount in 0:(product-1))
      {
        if (mdim > 1)
        {remain <- mcount
        for (mindex in mdim:2)
        {
          denom <- cumproduct[mindex-1]
          thisone[mindex] <- trunc(remain/denom)
          remain <- remain - denom*thisone[mindex]
        }
        thisone[1] <- remain
        }
        if (mdim == 1)
        {thisone[1] <- mcount}
        
        if (sum(thisone)%%2 ==0)
        {
          sortmoment <- sort(thisone)
          Tmoment <- make.moment.name(sortmoment)
          if (eval(parse(text=subblank(paste("!exists('",Tmoment,"',envir=",envir,",inherits=FALSE)")))))   
          {return(paste(Tmoment,' does not exist in environment ',envir))}
          
          Tmoment <- subblank(paste(envir,'$',Tmoment))
          
          eval(parse(text=paste("thismoment <- ",Tmoment)))       
          
          if (sum(thismoment$moment == thisone) < length(thismoment$moment))
          {thismoment <- tounsorted(thisone,thismoment)   
          }
          
          muproduct <- muproduct.Latex(moment-thismoment$moment)
          
          combinproduct <- prod(choose(moment,thismoment$moment))
          if (combinproduct == 1)
          {Acombinproduct <- " "}
          if (combinproduct > 1)
          {Acombinproduct <- subblank(paste(combinproduct))}
          thisLatex <- toLatex(thismoment)[-1] # get rid of initial expression
          
          if (length(thisLatex) > 1)
            
          {thisLatex[1] <- paste("(",thisLatex[1])
          thisLatex[length(thisLatex)] <- subblank(paste(thisLatex[length(thisLatex)],")"))
          thisLatex <- c(thisLatex," \\\\")
          }
          
          if (length(thisLatex) == 1)
          {#print(paste(thisLatex,thisLatex=="1"))
            if (subblank(thisLatex) != "1")
            {#print(paste("length",length(thisLatex),thisLatex))
              thisLatex[1] <- paste("(",thisLatex[1])
              thisLatex[length(thisLatex)] <- subblank(paste(thisLatex[length(thisLatex)],")"))
              thisLatex <- c(thisLatex," \\\\")
            }
            
            else {thisLatex <- c("\\\\")}
          }  
          if (mcount == 0)
          {thisLatex <- c(paste(Acombinproduct,"\\;",muproduct),thisLatex)}
          if (mcount > 0)
          {thisLatex <- c(paste("+",Acombinproduct,"\\;",muproduct),thisLatex)}
          
          fullmoment <-c(fullmoment,thisLatex)
        }
      }
    }
    return(fullmoment)
  }