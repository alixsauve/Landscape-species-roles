# feature_scaling() returns x after being scaled up between a and b
feature_scaling <- function(x, a, b){
  # x is a vector of numbers to be normalised
  # [a, b] is the range of x', a < b
  if (a < b){
    if (min(x) != max(x)){
      x2 <- a + (x-min(x))*(b-a)/(max(x)-min(x))
      return(x2)
    }
    else {
      stop("The extrema of x must be different values.")
    }
  }
  else {
    stop("a must be inferior to b.")
  }
}
