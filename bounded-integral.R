library(numDeriv) # Accurate numerical derivatives
library(ramify) # Additional matrix functionality
library(assertthat) # Pre and post assertion
library(testthat) # Unit testing

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
    msg = 'The domain of the density function is not bounded from below, however the limit value is nonzero. Please check your density function.'
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


test_file("tests/tests-bounded-integral.R")
