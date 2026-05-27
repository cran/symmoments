#' Numerically Integrate a Multivariate Polynomial
#'
#' Integrates a multivariate polynomial against a specified non-central 
#' multivariate normal distribution using ordinary numerical integration via 
#' the \code{\link[cubature]{adaptIntegrate}} function from the \code{cubature} package.
#'
#' @param poly An object of class \code{'mpoly'} or \code{'multipol'}, or a 
#'   simple list containing two components (\code{coeff} and \code{powers}) 
#'   defining the polynomial.
#' @param mu A numeric vector giving the mean vector \eqn{\mu} of the 
#'   multivariate normal distribution.
#' @param sigma A square matrix specifying the covariance matrix of the 
#'   multivariate normal distribution.
#' @param lower A numeric vector of the lower limits of integration, containing 
#'   one element for each dimension. If \code{NULL} (the default), it defaults 
#'   to \eqn{-6} times the standard deviations from the mean.
#' @param upper A numeric vector of the upper limits of integration, containing 
#'   one element for each dimension. If \code{NULL} (the default), it defaults 
#'   to \eqn{+6} times the standard deviations from the mean.
#'
#' @return The expected value of the polynomial numerically integrated against 
#' the specified multivariate normal distribution.
#'
#' @details Defaults for \code{lower} and \code{upper} boundaries are set to 
#' \eqn{\pm 6} times the standard deviations (the square roots of the diagonal 
#' elements of the covariance matrix \code{sigma}).
#' \cr\cr
#' If the polynomial is defined by a simple list, it must contain two components:
#' \itemize{
#'   \item \code{powers}: A matrix where each row represents the exponents/powers 
#'     for a single term in the polynomial.
#'   \item \code{coeff}: A numeric vector where each element is the coefficient 
#'     of the corresponding row in \code{powers}.
#' }
#' For example, the list structure equivalent to the polynomial in the examples section is:
#' \cr
#' \code{list(coeff = c(3, 2, -4, 1), powers = matrix(c(2,0,0, 1,3,0, 0,0,2, 1,2,1), ncol = 3, byrow = TRUE))}
#'
#' @references 
#' \insertRef{Phillips2010}{symmoments}
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{evaluate_expected.polynomial}}, \code{\link{multmoments}}, 
#'   \code{\link{evaluate}}, \code{\link{simulate}}
#'
#' @examples
#' \dontrun{
#' library(mpoly)
#' 
#' # Define an mpoly object for a multivariate polynomial
#' t0 <- mpoly(list(
#'   c(coef = 3, x1 = 2),
#'   c(coef = 2, x1 = 1, x2 = 3),
#'   c(coef = -4, z = 2),
#'   c(coef = 1, x1 = 1, x2 = 2, z = 1)
#' ))
#' 
#' # Numerically integrate against a specified mean and covariance identity matrix
#' integrate.polynomial(t0, c(1, 2, 3), matrix(c(1,0,0, 0,1,0, 0,0,1), nrow = 3, byrow = TRUE))
#' }
#'
#' @importFrom cubature adaptIntegrate
#' @export
`integrate.polynomial` <- 
  function (poly,mu,sigma,lower=NULL,upper=NULL) 
  {
    # integrate polynomial moment against MVN
    
    # poly: either a multipol objects or 
    #       a multipol defined by a list with moment powers and coefficients
    # mu: mean of multivariate normal as vector
    # sigma: variance-covariance matrix of multivariate normal
    # lower, upper: vectors giving limits of integration
    #    if one is NULL, then make it the mean +/- 6 * SD
    
    if (is.null(lower)) 
    {lower <- mu - 6*sqrt(diag(sigma))}
    if (is.null(upper))
    {upper <- mu + 6*sqrt(diag(sigma))}
    
    thispoly <- poly
    if (inherits(poly,"multipol"))
    {thispoly <- convert.multipol(poly)}
    if (inherits(poly,"mpoly"))
    {thispoly <- convert.mpoly(poly)}
    
    ndim <- dim(thispoly$powers)[2]
    npowers <- dim(thispoly$powers)[1]
    powers <- thispoly$powers
    coeff <- thispoly$coeff
    value <- 0
    
    f <- function(x)
    {
      y <- x[1]^powers[imom,1]
      for (idim in 2:ndim)
      {
        y <- y*x[idim]^powers[imom,idim]
      }
      y <- y*mvtnorm::dmvnorm(x,mean=mu,sigma=sigma, log=FALSE)
      return(y)
    }
    
    for (imom in 1:npowers)   
    {thisvalue <- adaptIntegrate(f,lower,upper)$integral
    value <- value + coeff[imom]*thisvalue
    }
    return(value)}