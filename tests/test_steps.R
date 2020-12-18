# Set working directory to be in the tests folder
#setwd('/Users/vicchoi/Documents/Classes/github_berkeley/stat-243-proj/tests')

library(testthat)
source('../utils.R')
source('../steps.R')

test_that("Initialization gives x values at correctly spaced points", {
  g1 <- function (x) {
    return(x^2)
  }
  h1 <- function (x) {
    return(log(g1(x)))
  }
  
  data1 <- initialization_step(h1, 1, 2)
  
  test_x <- seq(1, 2, length.out=22)[2:21]
  
  # Test that all values in the two arrays are near identical
  expect_true(all.equal(data1$x, test_x))
})


test_that("A sample from a normal gives mean of roughly 0", {
  set.seed(1)
  g2 <- function (x) {
    return(dnorm(x))
  }
  
  h2 <- function(x) {
    return(log(g2(x)))
  }
  
  data2 <- initialization_step(h2, -Inf, Inf)
  
  sample_x <- exp_sampling(10000, data2, h2)
  
  final_sample <- sampling_step(sample_x, data2, h2)
  
  # Test that the mean is roughly 0
  expect_true(abs(mean(final_sample$sample)) < 0.01)
  
})

test_that("Test that the update function correctly adds values and sorts", {
  set.seed(1)
  g2 <- function (x) {
    return(dnorm(x))
  }
  
  h2 <- function(x) {
    return(log(g2(x)))
  }
  
  data2 <- initialization_step(h2, -Inf, Inf)
  
  sample_x <- exp_sampling(50, data2, h2)
  
  final_sample <- sampling_step(sample_x, data2, h2)
  
  updated_data <- update_step(data2, final_sample)
  
  # Test that the updated data is of the proper length
  expect_true(length(updated_data$x) == length(final_sample$x) + length(data2$x))
  
  # Test that the data is sorted
  expect_true(!is.unsorted(updated_data$x))
  
})












