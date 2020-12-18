library(assertthat)
library(testthat)
source('../utils.R')
source('../steps.R')
source('../ars.R')

test_that("Test ARS on standard normal distribution, 
          comparing mean/std with theoretical values 0/1", {
  set.seed(1)
  g <- function (x) {
    return(dnorm(x))
  }
  
  ars_results <- ars(g, a = -Inf, b= Inf, N = 100000)
  
  expect_lt(abs(mean(ars_results) - 0), 0.01)
  expect_lt(abs(sd(ars_results) - 1), 0.01)
})

test_that("Test ARS on the (unnormalized) square function distribution (D = [0,1]), 
          comparing mean, std with theoretical values 3/4, 1/4 * sqrt(3/5)", {
  set.seed(1)
  g <- function (x) {
    return(x^2)
  }
  
  ars_results <- ars(g, a = 0, b = 1, N = 100000)
  
  expect_lt(abs(mean(ars_results) - 3/4), 0.01)
  expect_lt(abs(sd(ars_results) - 1/4 * sqrt(3/5)), 0.01)
})

test_that("Test ARS on the sine function distribution (D = [0,pi]), 
          comparing mean, std with theoretical values pi/2, 1/2 * sqrt(-8 + pi^2)", {
  set.seed(1)
  g <- function (x) {
    return(sin(x))
  }
  
  ars_results <- ars(g, a = 0, b = pi, N = 100000)
  
  expect_lt(abs(mean(ars_results) - pi/2), 0.01)
  expect_lt(abs(sd(ars_results) - 1/2 * sqrt(-8 + pi^2)), 0.01)
})






