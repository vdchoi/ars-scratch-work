library(numDeriv)

# We can use a package to do that.

# Example:

f <- function(x) {
  return(x^2)
}

x <- c(2.00, 3.00, 4.00, 5.00)

dx <- grad(f, x)
