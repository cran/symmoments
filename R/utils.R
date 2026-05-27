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
