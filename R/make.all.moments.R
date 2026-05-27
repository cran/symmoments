#' Create All Moments Up to Specified Size in Environment symmoments
#'
#' Create all central moment objects of a specified or smaller size in the 
#' \code{symmoments} environment.
#'
#' @param moment A vector \code{c(k1, ..., kn)} specifying the highest moment 
#'   to compute.
#' @param verbose If \code{TRUE} (the default), the names of the moments are 
#'   shown as the algorithm progresses; if \code{FALSE}, progress is hidden.
#'
#' @return All objects of class \code{'moment'} up to the value given in 
#' \code{moment} are created in the environment \code{symmoments}.
#'
#' @details Unsorted moments (those whose exponents are not in numeric order) 
#' are created in the \code{symmoments} environment using the 
#' \code{\link{tounsorted}} function to transform them from their sorted counterpart. 
#' If the \code{symmoments} environment does not exist, the user is prompted to 
#' create it using \code{symmoments <- new.env()}.
#' \cr\cr
#' If a sorted moment does not exist, it is automatically created. 
#' Moments of lower dimension are not created; for example, if \code{c(2, 4)} is 
#' input, \code{m20} is created, but \code{m2} is not.
#' \cr\cr
#' **Naming Conventions:**
#' \itemize{
#'   \item Moments are named using the structure \code{mij..l}, e.g., \code{m136}.
#'   \item If any exponent is greater than 9, lowercase letters, and then 
#'     uppercase letters are used. For example, \code{m3bA} represents the 
#'     moment \code{c(3, 11, 36)}.
#'   \item The largest single exponent allowed by this alphanumeric encoding 
#'     scheme is \eqn{9 + 26 + 26 = 61}.
#' }
#' If an object with a name of this form already exists in the target environment 
#' but is not an object of class \code{"moment"}, it will be silently overwritten 
#' by the new moment object.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{callmultmoments}}, \code{\link{tounsorted}}
#'
#' @examples
#' \dontrun{
#' # Create the symmoments environment if it does not exist 
#' symmoments <- new.env()
#' 
#' # Compute all moments up to c(3,3)
#' make.all.moments(c(3, 3))
#' }
#'
#' @export

`make.all.moments` <- 
  function (moment,verbose=TRUE) 
  {
    # make all moments less than or equal to input moment
    # put these into the symmoments namespace
    
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
    
    make.moment.vector <- function (moment) 
    {
      # returns character vector of moment
      # moment is the moment as a vector, eg, c(1,2,3)
      
      ndim <- length(moment)
      moment.vector <- subblank(paste("c(",moment[1]))
      if (ndim > 1)
      {
        moment.vector <- subblank(paste(moment.vector,","))
        for (idim in 2:ndim)
        {Ak <- moment[idim]
        moment.vector <- subblank(paste(moment.vector,Ak))
        if (idim < ndim)
        {moment.vector <- subblank(paste(moment.vector,","))}
        }
      }
      moment.vector <- subblank(paste(moment.vector,")"))
      return(moment.vector)
    }
    toSquare <- function(L.ut)
    {
      n <- (-1 + sqrt(1+8*length(L.ut)))/2
      L <- 0*diag(n)
      start <- 1
      for (irow in 1:n)
      {
        L[irow,] <- c(rep(0,(irow-1)),L.ut[start:(start+n-irow)])
        start <- start+n-irow + 1
      }
      return(L)}
    
    
    create.envir <- TRUE
    if (!exists('symmoments'))
    {symmoments <- NULL}
    if (exists('symmoments'))
    {
      if (inherits(symmoments,'environment'))
      {create.envir <- FALSE}
    }
    if (create.envir)
    {
      print('Environment symmoments must exist to receive the moment objects.')
      return('Please create it using   symmoments <- new.env()')
    }
    
    product = prod(moment+1)
    cumproduct = cumprod(moment+1)
    mdim = length(moment)
    thisone <- rep(0,mdim)
    for (mcount in 0:(product-1))
    {  #1
      remain <- mcount
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
      { #2
        notmoment <- 0
        Tmoment <- make.moment.name(thisone)
        if (exists(eval(parse(text=subblank(paste("'",Tmoment,"'")))),envir=symmoments,inherits=FALSE))
        { #3
          if (inherits(eval(parse(text=subblank(paste('symmoments$',Tmoment)))),'moment'))
          {
            notmoment <- 0
            if (verbose)
            {print(paste(subblank(subblank(paste("symmoments$",Tmoment)))," exists"))}
          }  # done with this one
          if (!inherits(eval(parse(text=subblank(paste("symmoments$",Tmoment)))),'moment'))
          {notmoment <- 1}
          
        } #3
        if (!exists(eval(parse(text=subblank(paste("'",Tmoment,"'")))),envir=symmoments,inherits=FALSE) | notmoment == 1)  
        { #3
          sortmoment <- sort(thisone)
          if (sum(sortmoment==thisone) == length(thisone))  
          {  #4
            thisvec <- make.moment.vector(thisone)
            if (verbose)
            {print(paste("Starting ",subblank(paste('symmoments$',Tmoment))))}
            eval(parse(text=subblank(paste("symmoments$",Tmoment," <- callmultmoments(",thisvec,")"))))
          }  #4  canonical, so just make it
          if (sum(sortmoment==thisone) != length(thisone))  #  unsorted
          { #4
            Smoment <- make.moment.name(sortmoment)
            if (exists(eval(parse(text=subblank(paste("'",Smoment,"'")))),envir=symmoments,inherits=FALSE))
            { #5
              if (inherits(eval(parse(text=subblank(paste("symmoments$",Smoment)))),'moment'))  # canonical moment exists
              { #6
                thisvec <- make.moment.vector(thisone)
                if (verbose)
                {print(paste("Starting ",subblank(paste("symmoments$",Tmoment))))}
                eval(parse(text=subblank(paste("symmoments$",Tmoment," <- tounsorted(",thisvec,",symmoments$",Smoment,")"))))
              } #6
              if (!inherits(eval(parse(text=subblank(paste("symmoments$",Smoment)))),'moment'))
              {notmoment <- 1}
            } #5
            
            if (!exists(eval(parse(text=subblank(paste("'",Tmoment,"'")))),envir=symmoments,inherits=FALSE) | notmoment == 1)
            { #5
              thisvec <- make.moment.vector(sortmoment)
              if (verbose)
              {print(paste("Starting ",subblank(paste("symmoments$",Smoment))," to create",subblank(paste("symmoments$",Tmoment))))}
              eval(parse(text=subblank(paste("symmoments$",Smoment," <- callmultmoments(",thisvec,")"))))
              thisvec <- make.moment.vector(thisone)
              if (verbose)
              {print(paste("Starting ",subblank(paste('symmoments$',Tmoment))))}
              eval(parse(text=subblank(paste("symmoments$",Tmoment," <- tounsorted(",thisvec,",symmoments$",Smoment,")"))))
            } #5
          } #4
        } #3
      } #2
    } #1
    return(NULL)
  }