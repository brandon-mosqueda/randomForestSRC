zap.rfsrc <- function() {
  g <- function(x, y) x / (1 - exp(-x)) - y
}

brandon <- function() {
  print("Brandon function from zap.rfsrc")
}