#' Print the Representation of a Multivariate Moment
#'
#' Prints an object of class \code{'moment'}.
#'
#' @param x An object of class \code{'moment'}, usually the output of 
#'   \code{\link{callmultmoments}}.
#' @param ... Included only for consistency with the generic function.
#'
#' @details Prints the moment as \code{E[X1**k1 X2**k2 ...]}: followed by 
#' the lines of the representation matrix with the corresponding coefficient 
#' attached to each row.
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{callmultmoments}}
#'
#' @examples
#' print(callmultmoments(c(1, 2, 3)))
#'
#' @method print moment
#' @export

`print.moment` <- 
  function(x, ...){
    
    # function: method to print a moment of the multivariate normal distribution
    # x is an object of class moment
    
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
    
    
    
    moment <- x$moment
    coef <- as.numeric(x$coefficient)
    representation  <- x$representation
    express <- "E["
    for (imom in 1:length(moment))
    {term <- subblank(paste("X",imom,"^",moment[imom])) 
    express <- paste(express,term)}
    express <- paste(express,"]:") 
    cat(express, "\n")
    print(cbind(coef,representation))
    
    invisible(x)}