#' Convert Between mpoly and List Representations of Multivariate Polynomials
#'
#' Converts an \code{mpoly} object to a simple list representation, or converts 
#' a simple list representation back to an \code{mpoly} object.
#'
#' @param poly An \code{mpoly} object, or a list containing powers and 
#'   coefficients that define a multivariate polynomial.
#'
#' @return If \code{poly} is of class \code{'mpoly'}, returns a list with two 
#' components (\code{powers} and \code{coeff}). If \code{poly} is a list of this 
#' form, returns the corresponding \code{mpoly} object.
#'
#' @details The list representation consists of 2 components:
#' \itemize{
#'   \item \code{powers}: A matrix where each row represents the exponents/powers 
#'     of \eqn{X} for a single term in the multivariate polynomial.
#'   \item \code{coeff}: A numeric vector where each element is the coefficient 
#'     for the corresponding row/term in \code{powers}.
#' }
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{convert.multipol}}, 
#'   \code{\link{evaluate_expected.polynomial}}, 
#'   \code{\link{integrate.polynomial}}
#'
#' @examples
#' \dontrun{
#' library(mpoly)
#' 
#' # Create an mpoly object
#' t0 <- mpoly::mpoly(list(
#'   c(coef = 3, x1 = 2),
#'   c(coef = 2, x1 = 1, x2 = 3),
#'   c(coef = -4, z = 2),
#'   c(coef = 1, x1 = 1, x2 = 2, z = 1)
#' ))  
#' 
#' # Convert from mpoly to list representation
#' t1 <- convert.mpoly(t0)    
#' 
#' # Convert from list representation back to an mpoly object
#' t2 <- convert.mpoly(t1) 
#' }
#'
#' @importFrom mpoly mpoly
#' @export
#' 
`convert.mpoly` <- 
  function (poly) 
  {
    # convert between a mpoly object and a list giving the corresponding moments and coefficients
    
    if (inherits(poly,"mpoly")) 
    {
      mpoly.list <- unclass(poly)
      matrix.size <- length(mpoly.list)
      
      poly.size <- prod(matrix.size)
      coeff <- rep(0,poly.size)
      variables <- "coef"
      vcount <- 0
      for (mcount in 1:poly.size)
      {   
        coeff[mcount] <- unlist(mpoly.list[[mcount]])['coef']
        variables <- c(variables,names(unlist(mpoly.list[mcount])))
      } 
      variables <- unique(variables)
      if (length(variables) == 1 & variables[1] == 'coef')
      {
        powers <-  matrix(0,ncol=1,nrow=1)
        return(list(coeff=coeff,powers=powers))
      }
      
      if (length(variables) > 1) 
      {
        variables <- variables[variables !='coef']
        ndim <- length(variables)
        powers <- matrix(rep(0,ndim*poly.size),ncol=ndim)
        
        for (mcount in 1:poly.size)
        {   
          thisone <- rep(0,ndim)
          thisexp <- unlist(unlist(mpoly.list[mcount]))
          thesevars <- c(names(unlist(mpoly.list[mcount])))
          nvars <- length(thesevars) - 1
          thesevars <- thesevars[thesevars != 'coef']
          nvars <- length(thesevars)
          if (nvars > 0)
          {for (idim in 1:ndim)
          {
            for (ivar in 1:nvars)
            {
              if (thesevars[ivar] == variables[idim])
              {
                thisone[idim] <- thisexp[ivar]
              }
            }
          }
          }
          powers[mcount,] <- thisone
        } 
        return(list(coeff=coeff,powers=powers))
      }
    }
    
    
    if (!inherits(poly,"mpoly"))  # assume conversion from matrix (ie, list) to mpoly
    { 
      n.powers <- dim(poly$powers)[1]
      if (!is.null(n.powers) | 0==0)
      {
        mpoly.list <- vector("list",length=n.powers)
        moment.size <- dim(poly$powers)[2]
        
        variables <- gsub(" ","",paste("X",1:moment.size))
        for (ipower in 1:n.powers)
        { 
          toeval <- " "
          thisterm <- poly$powers[ipower,]
          if (sum(thisterm) == 0)
          {mpoly.list[1] <- list(c(coef=poly$coeff[ipower]))} 
          
          if (sum(thisterm) > 0)
          {  
            thesevars <- variables[thisterm > 0]
            thesevalues <- thisterm[thisterm > 0] 
            toeval <- paste("c(",thesevars[1],"=",thesevalues[1])
            
            if (length(thesevars) > 1)
            { 
              for (ivar in 2:length(thesevars))
              {toeval <- paste(toeval,",",thesevars[ivar],"=",thesevalues[ivar])} 
            } 
            toeval    <- paste(toeval,", coef = ",poly$coeff[ipower],")")
            eval(parse(text=paste("mpoly.list[[ipower]] <- ",toeval)))  
          } 
        } 
        output <- mpoly.list
        return(output)
      }
    } 
    
  }


#' Convert Between multipol and List Representations of Multivariate Polynomials
#'
#' Converts a \code{multipol} object to a simple list representation, or converts 
#' a simple list representation back to a \code{multipol} object.
#'
#' @param poly A \code{multipol} object, or a list containing powers and 
#'   coefficients that define a multivariate polynomial.
#'
#' @return If \code{poly} is of class \code{'multipol'}, returns a list with two 
#' components (\code{powers} and \code{coeff}). If \code{poly} is a list of this 
#' form, returns the corresponding \code{multipol} object.
#'
#' @details The list representation consists of 2 components:
#' \itemize{
#'   \item \code{powers}: A matrix where each row represents the exponents/powers 
#'     of \eqn{X} for a single term in the multivariate polynomial.
#'   \item \code{coeff}: A numeric vector where each element is the coefficient 
#'     for the corresponding row/term in \code{powers}.
#' }
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{convert.mpoly}}, 
#'   \code{\link{evaluate_expected.polynomial}}, 
#'   \code{\link{integrate.polynomial}}
#'
#' @examples
#' \dontrun{
#' library(mpoly)
#' library(multipol)
#' 
#' # Create an mpoly object to work with
#' t0 <- mpoly::mpoly(list(
#'   c(coef = 3, x1 = 2),
#'   c(coef = 2, x1 = 1, x2 = 3),
#'   c(coef = -4, z = 2),
#'   c(coef = 1, x1 = 1, x2 = 2, z = 1)
#' )) 
#' 
#' # Convert from mpoly to list representation
#' t1 <- convert.mpoly(t0)    
#' 
#' # Convert from list representation to a multipol object
#' t2 <- convert.multipol(t1) 
#' 
#' # Convert back to a list representation
#' t3 <- convert.multipol(t2) 
#' }
#'
#' @export
`convert.multipol` <- function (poly) 
{
  # convert between a multipol object and a list giving the corresponding moments and coefficients
  
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
  if (inherits(poly,"multipol")) 
  {
    multipol.array <- as.array(poly)
    
    array.size <- dim(multipol.array)
    poly.size <- prod(array.size)
    moment.size <- length(dim(multipol.array)[dim(multipol.array)>1]) # problem of multipol array component x**0
    if (is.null(moment.size) | moment.size==0)
    {moment.size <- 1}
    
    coeff <- rep(0,poly.size)
    powers <- matrix(nrow=poly.size,ncol=moment.size)
    cumproduct = cumprod(dim(multipol.array))
    thisone <- rep(1,moment.size)
    
    nonzero <- 0
    for (mcount in 1:poly.size)
    {  
      thisone <- rep(1,moment.size)
      remain <- mcount
      if (moment.size == 1)
      {thisone <- c(remain)}
      if (moment.size > 1)
      {
        for (mindex in moment.size:2)
        {  
          thisone[mindex] <- trunc((remain-1)/cumproduct[mindex-1])
          remain <- remain - cumproduct[mindex-1]*thisone[mindex]
          thisone[mindex] <- thisone[mindex] + 1
        }  
      }
      thisone[1] <- remain
      
      powers[mcount,] <- thisone - rep(1,moment.size)  ###  account for indexing in multipol
      thistemp = paste(thisone[1])
      if (moment.size > 1)
      {
        for (mindex in 2:moment.size)
        {thistemp = paste(thistemp,",",thisone[mindex])}
      }
      thistemp = subblank(thistemp)
      coefftemp <- NULL
      eval(parse(text=paste("coefftemp <- multipol.array[",thistemp,"]")))
      if (coefftemp != 0)
      {
        nonzero <- nonzero + 1
        eval(parse(text=paste("coeff[nonzero] <- multipol.array[",thistemp,"]")))
        powers[nonzero,] <- thisone - rep(1,moment.size)  ###  account for indexing in multipol
      }
    } 
    coeff <- coeff[1:nonzero]
    powers <- as.matrix(powers[1:nonzero,],ncol=dim(powers)[2])
    output <- list(coeff=coeff,powers=powers)
    return(output)
  }
  
  
  if (!inherits(poly,"multipol"))  # assume conversion from matrix (ie, list) to multipol
  {
    n.powers <- dim(poly$powers)[1]
    moment.size <- dim(poly$powers)[2]
    if (dim(poly$powers)[1] == 1)
    {dimmult <- c(1,poly$powers + rep(1,moment.size))} 
    if (dim(poly$powers)[1]> 1)
    {dimmult <- apply(poly$powers, 2, max) + rep(1,moment.size)}
    n.elements <- prod(dimmult)
    multipol.array <- array(rep(0,n.elements),dim=dimmult)
    output <- as.multipol(multipol.array)
    for (ipower in 1:n.powers)
    {
      thistemp = paste(1+poly$powers[ipower,1])
      if (moment.size > 1)
      {
        for (mindex in 2:moment.size)
        {thistemp = paste(thistemp,",",(1+poly$powers[ipower,mindex]))}
      }
      thistemp = subblank(thistemp)
      eval(parse(text=paste("multipol.array[",thistemp,"] <- poly$coeff[",ipower,"]")))
    }
    output <- as.multipol(multipol.array) 
    
    return(output)
  }
  
  
}