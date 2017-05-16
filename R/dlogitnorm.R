#' Logit-normal distribution functions
#'
#' The logit-normal distribution is a probability distribution of a random variable
#' whose logit has a normal distribution. If Y is a random variable with a normal
#' distribution, and P is the logistic function, then X = P(Y) has a logit-normal
#' distribution; likewise, if X is logit-normally distributed, then
#' Y = logit(X)= log (X/(1-X)) is normally distributed.
#'
#' @author Daniel D Sjoberg \email{sjobergd@@mskcc.org}
#'
#' @param x input in [0,1]
#' @param n number of random variables to generate
#' @param location is the location parameter
#' @param scale is the scale parameter
#'
#' @return Returns a list of results from analysis.
#'
#' @examples
#' x=seq(0,1,0.1)
#' y=dlogitnorm(x)
#' plot(x,y)
#'
#'
#' @export
#'
# density function
dlogitnorm=function(x, location = 0, scale = 1){
  #function doesn't exist at 0 and 1, os making it zero
  ifelse(0<x & x<1,
         dnorm(qlogis(x), mean = location , sd = scale)*((x*(1-x))^-1),
         0)
}
# CDF
plogitnorm=function(x, location = 0, scale = 1){
  pnorm(qlogis(x), mean = location, sd = scale)
}
# random var generator
rlogitnorm=function(n, location = 0, scale = 1){
  exp(rnorm(qlogis(x), mean = location, sd = scale))
}
