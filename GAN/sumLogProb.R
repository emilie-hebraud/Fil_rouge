sumLogProb <- function(a, b) {
  #method that sums numbers in log scale without transforming back to non-log scale
  if (a == -Inf) {
    return(b)
  } else if (b == -Inf) {
    return(a)
  } else if (b < a) {
    return (a + log(1 + exp(b - a)))
  }
  return (b + log(1 + exp(a - b)))
}