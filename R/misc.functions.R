#' Misc. Functions
#' Error Function and Complimentary Error Function
#' Logit and Inverse Logit Function
#' @author Daniel D Sjoberg \email{sjobergd@@mskcc.org}
#'
#' @param x function input
#'
#' @return function output
#' @examples
#' erf(0)
#' @export
# error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

## (see Abramowitz and Stegun 29.2.29)
## and the 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

## and the inverses
inverf <- function (x) qnorm((1 + x)/2)/sqrt(2)
inverfc <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)


## logistic function
logit=function(x){
  log(x/(1-x))
}

invlogit=function(x){
  exp(x)/(1 + exp(x))
}
