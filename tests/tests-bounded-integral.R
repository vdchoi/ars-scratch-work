# Test file for the function "verify_concavity"

test_that("Test that density with bounded integral passes", {

  # Example function that has a positive finite integral

  # Normal distribution with mean 20.00 and stdev of 5.00
  f <- function(x) {
    return(1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x-20.0)/5.00)^2))
  }
  
  # Get output
  res <- verify_bounded_integral(f, -Inf, Inf)
  
  # The result should be 0
  expect_true(res == 0)

})


test_that("Test that density with unbounded integral fails", {
  
  # Abs(tan(x)) from 0 to pi theoretically integrates to Inf!

  f <- function(x) {
    return(abs(tan(x)))
  }
  
  # Get output
  res <- verify_bounded_integral(f, 0.00, pi)
  
  # The result should be 0
  expect_true(res == 1)
  
})


test_that("Test that density with unbounded integral fails (2)", {
  
  # Abs(tan(x)) from 0 to pi theoretically integrates to Inf!
  
  f <- function(x) {
    return(1/x)
  }
  
  # Get output
  res <- verify_bounded_integral(f, 1.00, Inf)
  
  # The result should be 0
  expect_true(res == 1)
  
})