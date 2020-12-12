library(numDeriv) # Accurate numerical derivatives
library(ramify) # Additional matrix functionality
library(assertthat) # Pre and post assertion
library(testthat) # Unit testing

verify_concavity <- function(logf, xmin, xmax, npoints) {
  #' Verify that a given function is concave
  #' (Input)
  #' logf: The function that we want to check whether 
  #' it's concave. It must take only one argument.
  #' The (log should have already been applied to the function).
  #' xmin, xmax: The domain of the function. Can be +-Inf.
  #' npoints: The number of samples that will be taken in the domain
  #' (Output)
  #' 0 if the function is log-concave, 1 otherwise
  
  # Assertions
  # Assert that logf is a function taking 1 argument
  assert_that(class(logf) == "function")
  assert_that(length(args(logf)) == 1)
  # Assert that xmin, xmax, npoints are scalars
  assert_that(is.numeric(xmin))
  assert_that(length(xmin) == 1)
  assert_that(is.numeric(xmax))
  assert_that(length(xmax) == 1)
  assert_that(is.numeric(npoints))
  assert_that(length(npoints) == 1)
  # Verify that xmin < xmax
  assert_that(xmin < xmax)
  # Verify that npoints > 100
  assert_that(npoints > 1)

  # If the domain is not bounded, check that the limit value
  # of the PDF is zero
  f <- function(x) exp(logf(x))
  if (xmin == -Inf) {
    assert_that(f(xmin) == 0,
    msg = 'The domain of the density function is not bounded from below, however the limit value is nonzero. Please check your density function.'
    )
  }
  if (xmax == Inf) {
    assert_that(f(xmax) == 0,
    msg = 'The domain of the density function is not bounded from above, however the limit value is nonzero. Please check your density function.'
    )
  }
  # Determine the mean and standard deviation
  f_mean <- function(x) x*exp(logf(x))
  mean <- integrate(f_mean, lower = xmin, upper = xmax)[[1]]
  f_var <- function(x) (x-mean)**2*exp(logf(x))
  variance <- integrate(f_var, lower = xmin, upper = xmax)[[1]]

  # Sample values to check log-concavity
  x_vec <- rnorm(npoints, mean=mean, sd = sqrt(variance))  # Sample
  # Restrict to lie in the domain of f
  # # keep as many samples as possible
  x_vec[x_vec < xmin] <- 2*xmin - x_vec[x_vec < xmin]
  x_vec[x_vec > xmax] <- 2*xmax - x_vec[x_vec > xmax]
  # and clip the rest, if any
  x_vec <- x_vec[x_vec > xmin]
  x_vec <- x_vec[x_vec < xmax]
  npoints <- length(x_vec)  # update npoints
  # sort
  x_vec <- sort(x_vec)
  
  # Evaluate the derivatives on those x values
  dy_vec <- grad(logf, x_vec)
  
  # Verify that all first derivatives are decreasing as x increases
  
  # Calculate the difference between the next and the previous elements
  # of the vector of derivatives
  dy_vec_shifted <- c(dy_vec[2:npoints], dy_vec[1])
  differences <- (dy_vec_shifted - dy_vec)[1:npoints-1]
  
  # If all differences are negative
  if (prod(differences < 0.00)) {
    # log-concave -> 0
    return(0)
  # Otherwise
  } else {
    # Not log-concave -> 1
    return(1)
  }
}


test_file("tests/tests-log-concavity.R")
