#' Evaluate a Moment or Polynomial
#'
#' @param object An object of class `symmoment` or a multi-dimensional polynomial.
#' @param sigma A numeric matrix representing the covariance matrix.
#'
#' @return The evaluated result.
#' @export
evaluate <- function(object, sigma) {
  UseMethod("evaluate", object)
}


#' Evaluate a Multivariate Moment
#'
#' Generic method for class \code{moment} to compute the numerical value of a 
#' moment at a specified covariance matrix from the output of 
#' \code{\link{callmultmoments}}.
#'
#' @param object An object of class \code{'moment'}.
#' @param sigma An upper-triangular matrix of covariance terms expressed as a 
#'   vector at which the moment is to be evaluated.
#'
#' @return The numeric value of the moment evaluated at the specified 
#' covariance matrix.
#'
#' @details \code{object} is normally the output of a call to 
#' \code{\link{callmultmoments}}. This is a list with the first component being 
#' the moment itself, the second component being the set of upper-triangular 
#' matrices representing the moment, and the third component containing their 
#' corresponding coefficients. This is an object of class \code{'moment'}.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{callmultmoments}}, and the methods \code{\link{simulate}} 
#'   and \code{\link{toLatex}} from the \code{symmoments} package.
#'
#' @examples
#' # Evaluates the moment at c(1,2,3,4) at the following covariance matrix:
#' #   4 2 1 1
#' #   2 3 1 1
#' #   1 1 2 1
#' #   1 1 1 2
#' evaluate(callmultmoments(c(1, 2, 3, 4)), c(4, 2, 1, 1, 3, 1, 1, 2, 1, 2))
#'
#' @method evaluate moment
#' @export

`evaluate.moment` <- 
  function (object, sigma) 
  {
    #  evaluate the moment using the representation from callmultmoment
    #      at the upper-triangular value of sigma
    
    #  object from callmultmoment
    #      list with first component the moment itself
    #      the second component the set of upper-triangular 
    #      matrices representing the moment
    #      and third component, their corresponding coefficients
    #
    #  sigma is the upper-triangular matrix of covariance terms
    #      at which the moment is to be evaluated
    #
    #  returns the value of the moment at this sigma
    
    moment <- object[[1]]
    moment.rep <- object[[2]]
    coefficients.rep <- object[[3]]
    
    #  evaluate the moment by adding the value at each representation
    #  this is the product of all sigma[i, j]^l[i, j] 
    #     if sigma and l are thought of as square matrices and l is the representation
    
    moment.value <- 0
    for (irep in 1:(dim(moment.rep)[1]))
    {moment.value <- moment.value + coefficients.rep[irep] * prod(sigma^moment.rep[irep, ])}
    
    return(as.vector(moment.value))}



#' Evaluate a Non-Central Multivariate Moment
#'
#' Computes the numerical value of a non-central moment at a specified mean 
#' and specified covariance matrix.
#'
#' @param moment A vector of non-negative integers representing the 
#'   non-central moment to be evaluated: \eqn{X_1^{k_1} \cdots X_n^{k_n}}.
#' @param mu A vector of real numbers representing the mean vector \eqn{\mu} 
#'   of the multivariate normal distribution.
#' @param sigma An upper-triangular matrix of covariance terms for the 
#'   multivariate normal distribution, expressed as a vector stacked by row, 
#'   at which the moment is to be evaluated.
#' @param envir A character string specifying the environment containing the 
#'   central moments needed for the calculation. Defaults to \code{'symmoments'}.
#'
#' @return The numeric value of the non-central moment evaluated at the 
#' specified mean and covariance matrix.
#'
#' @details This function searches the environment specified in the \code{envir} 
#' argument for the central moments required to complete the non-central expansion. 
#' The default is the \code{symmoments} environment. All even central moments 
#' less than or equal to the requested \code{moment} vector must be present. 
#' The computation will stop with an error message if any required central moment 
#' is missing from \code{envir}.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{evaluate.moment}}, \code{\link{make.all.moments}}
#'
#' @examples
#' \dontrun{
#' # Evaluates the expected value of X1^3 X2 X3^2 at mean c(3,4,1) 
#' # and at the following covariance matrix:
#' #    4 2 1 
#' #    2 3 1 
#' #    1 1 2 
#' # Note: requires all central moments up to c(3,1,2) to exist in 'symmoments'.
#' # If needed, run: make.all.moments(c(3,1,2))
#' evaluate_noncentral(c(3, 1, 2), c(3, 4, 1), c(4, 2, 1, 3, 1, 2))
#' 
#' # Using central moments stored instead in the global environment:
#' evaluate_noncentral(c(3, 1, 2), c(3, 4, 1), c(4, 2, 1, 3, 1, 2), '.GlobalEnv')
#' }
#'
#' @export
`evaluate_noncentral` <- 
  function (moment,mu,sigma,envir='symmoments') 
  {
    # Evaluate noncentral moment with mean mu and covariance matrix sigma
    # moment and mu are vectors
    # sigma is a vector of the upper diagonal
    # envir is a character variable containing the name of the environment containing the required central moments
    
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
    
    product = prod(moment+1)
    cumproduct = cumprod(moment+1)
    mdim = length(moment)
    thisone <- rep(0,mdim)
    value <- 0
    
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
        if (mdim > 1)
        {
          if (sum(thismoment$moment == thisone) < length(thismoment$moment))
          {thismoment <- tounsorted(thisone,thismoment)} 
        }
        muproduct <- prod(mu^(moment-thismoment$moment))
        combinproduct <- prod(choose(moment,thismoment$moment))
        thisvalue <- combinproduct*evaluate.moment(thismoment,sigma)*muproduct
        value <- value + thisvalue
      }
    }
    return(value)
  }

#' Evaluate the Expected Value of a Multivariate Polynomial
#'
#' Evaluates the expected value of a multivariate polynomial assuming a 
#' specified non-central multivariate normal distribution.
#'
#' @param poly An object of class \code{'mpoly'} or \code{'multipol'}, or a simple 
#'   list containing coefficients and powers defining the polynomial.
#' @param mu A vector of real numbers representing the mean vector \eqn{\mu} 
#'   of the multivariate normal distribution.
#' @param sigma A vector giving an upper-triangular matrix, stacked by row, 
#'   representing the covariance matrix of the multivariate distribution.
#' @param envir A character string specifying the environment containing the 
#'   central moments needed for the calculation. Defaults to \code{'symmoments'}.
#'
#' @return The expected value of the multivariate polynomial evaluated at the 
#' specified multivariate normal mean and covariance matrix.
#'
#' @details This function searches the environment specified in the \code{envir} 
#' argument for the central moments required to complete the expected value expansion. 
#' The default is the \code{symmoments} environment. The computation will stop 
#' with an error message if any required central moment is missing from \code{envir}.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{evaluate_noncentral}}, \code{\link{make.all.moments}}
#'
#' @examples
#' \dontrun{
#' library(mpoly)
#' 
#' # Define an mpoly object for a multivariate polynomial and determine
#' # its expected value at a specified mean and covariance matrix:
#' # Note: All moments up to c(2,3,2) must exist in the symmoments environment. 
#' # Run make.all.moments(c(2,3,2)) beforehand if necessary.
#' 
#' t0 <- mpoly(list(
#'   c(coef = 3, x1 = 2),
#'   c(coef = 2, x1 = 1, x2 = 3),
#'   c(coef = -4, z = 2),
#'   c(coef = 1, x1 = 1, x2 = 2, z = 1)
#' ))
#' 
#' evaluate_expected.polynomial(t0, c(1, 2, 3), c(1, 0, 0, 1, 0, 1))
#' }
#'
#' @export
#' 
`evaluate_expected.polynomial` <- 
  function (poly,mu,sigma,envir='symmoments') 
  {
    # compute expected value of a multidimensional polynomial
    # poly is either a multipol objects or
    # a list with components "powers", which lists the moments involved,
    # and "coeff, their coefficients
    # mu is the mean as a vector
    # sigma is the variance-covariance matrix as the vector of the upper diagonal part, including diagonal
    # envir is a character variable containing the name of the environment containing the required central moments
    
    exists.envir <- FALSE
    if (exists(envir))
    {
      if (inherits(eval(parse(text=envir)),'environment'))
      {exists.envir <- TRUE}
    }
    if (!exists.envir)
    {return(paste('There is no environment named ', envir)) }
    
    
    temp.poly <- poly
    if (inherits(poly,"multipol"))
    {temp.poly <- convert.multipol(poly)}
    if (inherits(poly,"mpoly"))
    {temp.poly <- convert.mpoly(poly)}
    
    
    ndim <- dim(temp.poly$powers)[2]
    npowers <- dim(temp.poly$powers)[1]
    powers <- temp.poly$powers
    coeff <- temp.poly$coeff
    
    value <- 0
    for (imom in 1:npowers)  
    {
      if (coeff[imom] != 0)
      {
        thisvalue <- coeff[imom]*evaluate_noncentral(powers[imom,],mu,sigma,envir)
        value <- value + thisvalue
      }
    }
    
    return(value)
  }
