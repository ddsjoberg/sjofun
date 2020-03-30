#' One-diminsional Search Algorithm
#'
#' Find the minimum of a function using the Golded-section search algorithm
#'
#' @param x function with one argument
#' @param a min value of search values
#' @param b max value of search values
#' @param tol tolerance level. Default is `1e-04`
#'
#' @export
#' @examples
#' contfunc = function(x){
#'   -exp(-x) * sin(x)
#' }
#' golden_search(contfunc, a = 0, b = 1.5)
golden_search <- function(x, a, b, tol = 1e-04) {
  phi <- (-1 + sqrt(5)) / 2

  c <- b - phi * (b - a)
  d <- a + phi * (b - a)
  i <- 0
  iter <- NULL
  while (abs(c - d) > tol) {
    i <- i + 1

    fc <- x(c)
    fd <- x(d)

    # to find min, use < ; > is for max
    if (fc < fd) {
      b <- d
      d <- c
      # fd=fc;fc=f(c)
      c <- b - phi * (b - a)
    } else {
      a <- c
      c <- d
      # fc=fd;fd=f(d)
      d <- a + phi * (b - a)
    }

  }
  return((b + a) / 2)
}
