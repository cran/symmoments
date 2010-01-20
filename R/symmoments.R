
`multmoments` <- 
function (moment, current.matrix, current.cell, moment.rep, row_col) 
{
# A recursive function to compute the representation of a multivariate moment
#    using upper - triangular matrices
#
# moment: a vector of integers representing the moment
#     eg, c(3, 2, 4) for a 3-dimensional normal vector
#     corresponding to the moment E[(X1^3)(X2^2)(X3^4)]

# current.matrix: input/output matrix under consideration in recursion
#     this is an upper-triangular integer matrix (l(i, j))
#     organized by row, with   lmom * (lmom + 1) / 2   total elements  
#
# current.cell: input/output cell in current matrix under consideration in recursion

# moment.rep: input/output current set of representations
#     function adds each satisfying matrix to moment.rep

# row_col: 2x(lmom * (lmom + 1) / 2 matrix giving rows and columns for square matrix
#     for each cell so that they don't have to be calculated each time

#  algorithm:
#  loop through cell values from 0 to min(moments[row] - rowsum, moments[col] - colsum)
#    if the this new matrix satisfies moment criterion, 
#       add to moment.rep and return
#    if the current matrix is too great in any dimension, 
#       return
#    if the current matrix is < moment in any dimension, and
#       at most moment(i) for all other indexs i, continue

summomentmatrix <- function (moment.matrix) 
{
#  compute the row/col sums of moment.matrix
#  uses Matrix package
#  construct a square upper diagonal matrix from moment.matrix

#  length of moment based on representation:
lmom <- (sqrt(8 * length(moment.matrix) + 1) - 1) / 2 
#  make a square matrix for use with row and column sum functions
tempmatrix <- matrix(c(rep(0, lmom^2)), nrow=lmom)
tempmatrix[1, ] <- moment.matrix[1:lmom]
endrow <- lmom
if (lmom > 1)
  {
  for (irow in 2:lmom)
     {startrow <- endrow + 1
      endrow <- endrow + (lmom - irow + 1)
      tempmatrix[irow, ] <- cbind(c(rep(0, (irow - 1)), moment.matrix[startrow:endrow])) }
     
  summomentmatrix <- colSums(tempmatrix) + rowSums(tempmatrix) 
  }
if (lmom == 1){summomentmatrix <- 2 * sum(moment.matrix)}     

return(summomentmatrix)}


lmom <- length(moment)
totcells <- lmom * (lmom + 1) / 2   # total cells in a moment representation
 
thisrow <- row_col[1, current.cell]
thiscol <- row_col[2, current.cell]
 
moment.row <- moment[thisrow]
moment.col <- moment[thiscol]
rowcells <- row_col[1, row_col[1, ] == thisrow]
colcells <- row_col[2, row_col[2, ] == thiscol]
 
#  determine sums for row and columm
#  maxvalue is the largest value that can be added to a
#    a row/column sum and still be no more than criterion
#  maxvalue will be the minimum of the moments minus these sums

rowsum <- summomentmatrix(current.matrix)[thisrow]
colsum <- summomentmatrix(current.matrix)[thiscol] 
 
maxvalue <- min(moment[thisrow] - rowsum, moment[thiscol] - colsum)
 
for (ivalue in 0:maxvalue)
   {         # for:loop through all possible values for cell
 
   current.matrix[current.cell] <- ivalue
   current.sum <- summomentmatrix(current.matrix)
             # determine current sum for this matrix

#  if matrix fulfills criterion, add to reps and return
   if (sum(moment == current.sum) == lmom)
      {moment.rep <- rbind(moment.rep, current.matrix)
       return(moment.rep)}

#   if the sum is too large, return because any other sum will also be too large
   if (sum(current.sum > moment) > 0){return(moment.rep)}

#   if at least one term in current sum is smaller than moment, 
#       go down one cell unless this is the last cell

   if ((sum(current.sum == moment) < lmom) & (current.cell != totcells))     
      {cc <- current.cell + 1
# recursive step:
       moment.rep <- multmoments(moment, current.matrix, cc, moment.rep, row_col)}
                        
   }        # end of for 
 
#   note: the three conditions above are exhaustive except for last cell 
#      In that case there is nothing to return
return(moment.rep)
}

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

if (class(moment.fullrep) == "moment"){thismoment <- moment.fullrep$moment}
if (class(moment.fullrep) != "moment")
   {print("moment must be of class 'moment'")
    return(-1)}   
    
ndim <- length(thismoment)                                                                                        
sample <- rmvnorm(n=nsim, mean=Mean, sigma=matrix(Sigma,nrow=length(Mean)))
exponents <- matrix(rep(thismoment, nsim), nrow=nsim, byrow=TRUE)
powers <- sample^exponents
prods <- rep(1, nsim)        #  calculate product of powers of Xs
for (icol in (1:ndim))
   {prods <- prods * powers[, icol]}
moment.value <- mean(prods)
 
return(moment.value)}



`evaluate` <- function(object, sigma) 
   {UseMethod("evaluate",object)}

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
  

