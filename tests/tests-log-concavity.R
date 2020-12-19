# Test file for the function "verify_concavity"
source('../function_checks.R')

test_that("Test that density with bounded integral passes", {
  
  # Example function that has a positive finite integral
  
  # Normal distribution with mean 20.00 and stdev of 5.00
  f <- function(x) {
    return(1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x-20.0)/5.00)^2))
  }
  
  # Get output
  res <- verify_bounded_integral(f, -Inf, Inf)
  
  expect_true(res)
  
})


test_that("Test that density with unbounded integral fails", {
  
  # Abs(tan(x)) from 0 to pi theoretically integrates to Inf!
  
  f <- function(x) {
    return(abs(tan(x)))
  }
  
  # Get output
  res <- verify_bounded_integral(f, 0.00, pi)
  
  expect_false(res)
  
})


test_that("Test that density with unbounded integral fails (2)", {
  
  # Abs(tan(x)) from 0 to pi theoretically integrates to Inf!
  
  f <- function(x) {
    return(1/x)
  }
  
  # Get output
  res <- verify_bounded_integral(f, 1.00, Inf)
  
  expect_false(res)
  
})



# Test file for the function "verify_log_concavity"
test_that("Test that log-concave function passes (1)", {

  # Example function that is log-concave

  # Normal distribution with mean 20.00 and stdev of 5.00
  f <- function(x) {
    return(1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x-20.0)/5.00)^2))
  }
  
  # Get output
  res <- verify_log_concavity(f, -Inf, Inf)
  
  expect_true(res)

})

test_that("Test that log-concave function passes (2)", {
  
  # Example function that is log-concave
  
  # Linear
  f <- function(x) {
    return(0.20*x+5.00)
  }
  
  # Get output
  res <- verify_log_concavity(f, -10, 10)
  
  expect_true(res)
  
})


test_that("Test that log-concave function passes (3)", {
  
  # Sin
  f <- function(x) {
    return(0.50*sin(x))
  }
  
  # Get output
  res <- verify_log_concavity(f, 0.0, pi)
  
  expect_true(res)
  
})



test_that("Test that not log-concave function fails", {
  
  # Example function that is not log-concave
  g <- function(x) {
    return(4.00^x)
  }
  
  # Get output
  res <- verify_log_concavity(g, 0.00, 1.00)
  
  expect_false(res)
  
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
  
  expect_false(res)
  
})

test_that("Test that not log-concave function fails (3)", {
  
  # Example function that is not log-concave
  f <- function(x) {
    return(
      tan(x)
    )
  }
  
  # Get output
  res <- verify_log_concavity(f, 0.00+0.10, pi/2.-0.10)
  
  expect_false(res)
  
})

test_that("Test that non-positive density fails", {
  
  # Example function that evaluates to non-positive values
  
  f <- function(x) {
    return(0.20*x-5.00)
  }
  
  # Get output
  res <- verify_log_concavity(f, -10, 10)
  
  expect_false(res)
  
})

