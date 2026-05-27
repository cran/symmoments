
<!-- README.md is generated from README.Rmd. Please edit that file -->

# symmoments

<!-- badges: start -->

[![R-CMD-check](https://github.com/FloSchuberth/symmoments/actions/workflows/rhub.yaml/badge.svg)](https://github.com/FloSchuberth/symmoments/actions)
<!-- badges: end -->

Symbolically computes and numerically evaluates multivariate normal
moments , where , in terms of mu and S elements.

## Installation

You can install the development version of symmoments from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("FloSchuberth/symmoments")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(symmoments)
 # Compute the moment for the 4-dimensional moment c(1,2,3,4):
callmultmoments(c(1,2,3,4))
#> E[ X1^1 X2^2 X3^3 X4^4 ]: 
#>    coef S(1,1) S(1,2) S(1,3) S(1,4) S(2,2) S(2,3) S(2,4) S(3,3) S(3,4) S(4,4)
#> 1    72      0      0      0      1      0      0      2      1      1      0
#> 2   144      0      0      0      1      0      1      1      0      2      0
#> 3    72      0      0      0      1      0      1      1      1      0      1
#> 4    72      0      0      0      1      0      2      0      0      1      1
#> 5    24      0      0      0      1      1      0      0      0      3      0
#> 6    36      0      0      0      1      1      0      0      1      1      1
#> 7    72      0      0      1      0      0      0      2      0      2      0
#> 8    36      0      0      1      0      0      0      2      1      0      1
#> 9   144      0      0      1      0      0      1      1      0      1      1
#> 10   18      0      0      1      0      0      2      0      0      0      2
#> 11   36      0      0      1      0      1      0      0      0      2      1
#> 12    9      0      0      1      0      1      0      0      1      0      2
#> 13   48      0      1      0      0      0      0      1      0      3      0
#> 14   72      0      1      0      0      0      0      1      1      1      1
#> 15   72      0      1      0      0      0      1      0      0      2      1
#> 16   18      0      1      0      0      0      1      0      1      0      2

# Print the representation of the 4-dimensional moment c(1,2,3,4):
print(callmultmoments(c(1,2,3,4)))
#> E[ X1^1 X2^2 X3^3 X4^4 ]: 
#>    coef S(1,1) S(1,2) S(1,3) S(1,4) S(2,2) S(2,3) S(2,4) S(3,3) S(3,4) S(4,4)
#> 1    72      0      0      0      1      0      0      2      1      1      0
#> 2   144      0      0      0      1      0      1      1      0      2      0
#> 3    72      0      0      0      1      0      1      1      1      0      1
#> 4    72      0      0      0      1      0      2      0      0      1      1
#> 5    24      0      0      0      1      1      0      0      0      3      0
#> 6    36      0      0      0      1      1      0      0      1      1      1
#> 7    72      0      0      1      0      0      0      2      0      2      0
#> 8    36      0      0      1      0      0      0      2      1      0      1
#> 9   144      0      0      1      0      0      1      1      0      1      1
#> 10   18      0      0      1      0      0      2      0      0      0      2
#> 11   36      0      0      1      0      1      0      0      0      2      1
#> 12    9      0      0      1      0      1      0      0      1      0      2
#> 13   48      0      1      0      0      0      0      1      0      3      0
#> 14   72      0      1      0      0      0      0      1      1      1      1
#> 15   72      0      1      0      0      0      1      0      0      2      1
#> 16   18      0      1      0      0      0      1      0      1      0      2

# Compute the LaTeX representation of the central moment c(1,2,3,4):
toLatex(callmultmoments(c(1,2,3,4)))
#>  [1] "E[X_{1}^{1}X_{2}^{2}X_{3}^{3}X_{4}^{4}]=\\\\"                             
#>  [2] "18\\sigma_{1,2}\\sigma_{2,3}\\sigma_{3,3}\\sigma_{4,4}^{2}+"              
#>  [3] "72\\sigma_{1,2}\\sigma_{2,3}\\sigma_{3,4}^{2}\\sigma_{4,4}+"              
#>  [4] "72\\sigma_{1,2}\\sigma_{2,4}\\sigma_{3,3}\\sigma_{3,4}\\sigma_{4,4}+"     
#>  [5] "48\\sigma_{1,2}\\sigma_{2,4}\\sigma_{3,4}^{3}+\\\\"                       
#>  [6] "9\\sigma_{1,3}\\sigma_{2,2}\\sigma_{3,3}\\sigma_{4,4}^{2}+"               
#>  [7] "36\\sigma_{1,3}\\sigma_{2,2}\\sigma_{3,4}^{2}\\sigma_{4,4}+"              
#>  [8] "18\\sigma_{1,3}\\sigma_{2,3}^{2}\\sigma_{4,4}^{2}+"                       
#>  [9] "144\\sigma_{1,3}\\sigma_{2,3}\\sigma_{2,4}\\sigma_{3,4}\\sigma_{4,4}+\\\\"
#> [10] "36\\sigma_{1,3}\\sigma_{2,4}^{2}\\sigma_{3,3}\\sigma_{4,4}+"              
#> [11] "72\\sigma_{1,3}\\sigma_{2,4}^{2}\\sigma_{3,4}^{2}+"                       
#> [12] "36\\sigma_{1,4}\\sigma_{2,2}\\sigma_{3,3}\\sigma_{3,4}\\sigma_{4,4}+"     
#> [13] "24\\sigma_{1,4}\\sigma_{2,2}\\sigma_{3,4}^{3}+\\\\"                       
#> [14] "72\\sigma_{1,4}\\sigma_{2,3}^{2}\\sigma_{3,4}\\sigma_{4,4}+"              
#> [15] "72\\sigma_{1,4}\\sigma_{2,3}\\sigma_{2,4}\\sigma_{3,3}\\sigma_{4,4}+"     
#> [16] "144\\sigma_{1,4}\\sigma_{2,3}\\sigma_{2,4}\\sigma_{3,4}^{2}+"             
#> [17] "72\\sigma_{1,4}\\sigma_{2,4}^{2}\\sigma_{3,3}\\sigma_{3,4}\\\\"

# evaluate the moment c(1,2,3,4) at the following variance-covariance matrix
#  4 2 1 1
#  2 3 1 1
#  1 1 2 1
evaluate(callmultmoments(c(1,2,3,4)), c(4,2,1,1,3,1,1,2,1,2))
#> [1] 3480
```
