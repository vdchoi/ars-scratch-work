# Test file for the function "verify_log_concavity"

test_that("Test that concave function passes", {

  # Example function that is log-concave

  # Normal distribution with mean 20.00 and stdev of 5.00
  f <- function(x) {
    return(1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x-20.0)/5.00)^2))
  }
  
  # Get output
  res <- verify_log_concavity(f, -Inf, Inf)
  
  # The result should be 0
  expect_true(res == 0)

})

test_that("Test that concave function passes", {
  
  # Example function that is log-concave
  
  # Linear
  f <- function(x) {
    return(0.20*x+5.00)
  }
  
  # Get output
  res <- verify_log_concavity(f, -10, 10)
  
  # The result should be 0
  expect_true(res == 0)
  
})

test_that("Test that not log-concave function fails", {
  
  # Example function that is not log-concave
  g <- function(x) {
    return(4.00^x)
  }
  
  # Get output
  res <- verify_log_concavity(g, 0.00, 1.00)
  
  # The result should be 1
  expect_true(res == 1)
  
})

test_that("Test that not log-concave function fails (2)", {
  
  # Example function that is not log-concave
  # Sum of two Gaussians
  f <- function(x) {
    return(
      1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x-20.0)/5.00)^2) + 
      1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x+20.0)/5.00)^2)
      )
  }
  
  # Get output
  res <- verify_log_concavity(f, -Inf, Inf)
  
  # The result should be 0
  expect_true(res == 1)
  
})

test_that("Test that not log-concave function fails (3)", {
  
  # Example function that is not log-concave
  f <- function(x) {
    return(
      tan(x)
    )
  }
  
  # Get output
  res <- verify_log_concavity(f, 0.10, pi/2.-0.10)
  
  # The result should be 0
  expect_true(res == 1)
  
})

test_that("Test that non-positive density fails", {
  
  # Example function that evaluates to non-positive values
  
  f <- function(x) {
    return(0.20*x-5.00)
  }
  
  # Get output
  res <- verify_log_concavity(f, -10, 10)
  
  # The result should be 0
  expect_true(res == 1)
  
})