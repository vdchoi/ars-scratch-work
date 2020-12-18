library(numDeriv) # Accurate numerical derivatives
library(ramify) # Additional matrix functionality
library(assertthat) # Pre and post assertion
library(testthat) # Unit testing
library(docstring)


verify_bounded_integral <- function(func, a, b) {
  #' Verify that a given function has a finite positive integral
  #' 
  #' @param func The function that we want to check
  #' @param a The lower bound of the domain of the function. Can be -Inf.
  #' @param b The upper bound of the domain of the function. Can be +Inf.
  #' @return 0 if the function has a finite positive integral, 1 otherwise
  
  # Assertions
  # Assert that func is a function taking 1 argument
  assert_that(class(func) == "function")
  assert_that(length(args(func)) == 1)
  assert_that(is.numeric(a))
  assert_that(length(a) == 1)
  assert_that(is.numeric(b))
  assert_that(length(b) == 1)
  # Verify that a < b
  assert_that(a < b)
  
  
  # If the domain is not bounded, check that the limit value
  # of the PDF is zero
  
  if (a == -Inf) {
    assert_that(func(a) == 0,
                msg = 'The domain of the density function is not bounded from below, however the limit value is nonzero or indeterminate. Please check your density function.'
    )
  }
  if (b == Inf) {
    assert_that(func(b) == 0,
                msg = 'The domain of the density function is not bounded from above, however the limit value is nonzero. Please check your density function.'
    )
  }
  
  result = tryCatch({
    # See if it can be done!
    c <- integrate(func, lower = a, upper = b)[[1]]
    # If it can, then check if the integral is positive
    if (c > 0) {
      return(0)
    } else {
      # Negative integral?!
      return(1)
    }
  }, error = function(error_condition) {
    print(error_condition)
    return(1)
  })
  
}



verify_log_concavity <- function(func, a, b, npoints = 10000) {
  #' Verify that a given function is log-concave and always positive
  #' 
  #' @param func The function that we want to check whether it's log-concave. It must take only one argument.
  #' @param a The lower bound of the domain of the function. Can be -Inf.
  #' @param b The upper bound of the domain of the function. Can be +Inf.
  #' @param npoints: The number of samples that will be taken in the domain
  #' @return 0 if the function is log-concave, 1 otherwise
  
  # Assertions
  assert_that(is.numeric(npoints))
  assert_that(length(npoints) == 1)
  # Verify that npoints > 100
  assert_that(npoints > 100)

  # Determine sampling points to check log-concavity
  # Approximate the mean and standard deviation
  h <- function(x) log(func(x))
  f_mean <- function(x) x*func(x)
  mean <- integrate(f_mean, lower = a, upper = b)[[1]]
  f_var <- function(x) (x-mean)**2*func(x)
  variance <- integrate(f_var, lower = a, upper = b)[[1]]
  if (variance < 0.00) {
    return(1) # The density can't be negative
  }
  

  # Sample values to check log-concavity
  x_vec <- rnorm(npoints, mean=mean, sd = sqrt(variance))  # Sample
  # Restrict to lie in the domain of f
  # # keep as many samples as possible
  x_vec[x_vec < a] <- 2*a - x_vec[x_vec < a]
  x_vec[x_vec > b] <- 2*b - x_vec[x_vec > b]
  # and clip the rest, if any
  x_vec <- x_vec[x_vec > a]
  x_vec <- x_vec[x_vec < b]
  npoints <- length(x_vec)  # update npoints
  # sort
  x_vec <- sort(x_vec)
  
  # Evaluate the derivatives on those x values
  dy_vec <- grad(h, x_vec)

  # Verify that all first derivatives are decreasing as x increases
  
  # Calculate the difference between the next and the previous elements
  # of the vector of derivatives
  differences <- diff(dy_vec)

  # Evaluate the function at those x values
  f_vec <- func(x_vec)
  
  # If all differences are negative (log-concave) AND its always positive
  if (prod(differences < 0.00) & prod(f_vec>0.00)) {
    # log-concave -> 0
    return(0)
  # Otherwise
  } else {
    # Not log-concave -> 1
    return(1)
  }
}
