#' Compute a Multivariate Moment Using Monte Carlo Integration
#'
#' Computes a multivariate normal moment by Monte Carlo integration.
#'
#' @param object An object of class \code{'moment'} representing 
#'   \eqn{E[X_1^{k_1} \cdots X_n^{k_n}]}.
#' @param nsim The number of samples to generate in computing the integral.
#' @param seed An integer for the random number generator (\code{\link[base]{set.seed}}).
#' @param Mean The mean vector of \eqn{(X_1, \dots, X_n)}.
#' @param Sigma Covariance matrix of \eqn{(X_1, \dots, X_n)}, dimension \eqn{n \times n}, 
#'   expressed as a vector stacked by row.
#' @param ... Included only for consistency with the generic function.
#'
#' @return An approximate numerical value of the specified moment.
#'
#' @note Non-central moments can be approximated by specifying \code{Mean}. 
#' For central moments, set \code{Mean} to a vector of zeros.
#' \cr\cr
#' The \code{mvtnorm} package must be installed for this function to utilize 
#' \code{\link[mvtnorm]{rmvnorm}}.
#'
#' @references 
#' \insertRef{Rizzo2008}{symmoments}
#' 
#'
#' @author Kem Phillips \email{kemphillips@@comcast.net}
#'
#' @seealso \code{\link{callmultmoments}}, and the methods \code{\link{toLatex}} 
#'   and \code{\link{evaluate}}.
#'
#' @examples
#' # Using 10000 samples, estimate the central moment for the moment c(2,4) 
#' # at the specified covariance matrix and mean (0,0):
#' library(mvtnorm)
#' simulate(callmultmoments(c(2, 4)), 10000, NULL, c(0, 0), c(2, 1, 1, 4))
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats simulate
#' @method simulate moment
#' @export

`simulate.moment` <- 
  function(object, nsim, seed=NULL, Mean, Sigma, ...){
    
    # function: method to calculate moment of the multivariate normal distribution
    #           using Monte-Carlo integration (Rizzo, 2008)
    # object is an object of class moment
    # nsim is the number of samples to generate
    # seed is the seed for the random number generator
    # Mean is the mean of the (X1, ..., Xn)
    # Sigma is the variance-covariance of (X1^k1, ..., Xn^kn), dimension nXn
    
    
    # requires package mvtnorm for function rmvnorm
    
    moment.fullrep <- object
    if (is.numeric(seed)){set.seed(seed)}
    
    if (inherits(moment.fullrep,"moment")){thismoment <- moment.fullrep$moment}
    if (!inherits(moment.fullrep,"moment"))
    {print("moment must be of class 'moment'")
      return(-1)}   
    
    ndim <- length(thismoment)                                                                                        
    sample <- mvtnorm::rmvnorm(n=nsim, mean=Mean, sigma=matrix(Sigma,nrow=length(Mean)))
    exponents <- matrix(rep(thismoment, nsim), nrow=nsim, byrow=TRUE)
    powers <- sample^exponents
    prods <- rep(1, nsim)        #  calculate product of powers of Xs
    for (icol in (1:ndim))
    {prods <- prods * powers[, icol]}
    moment.value <- mean(prods)
    
    return(moment.value)}