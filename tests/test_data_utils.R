setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(assertthat)
library(testthat)
source('../utils/data_utils.R')

test_that("The upper bound and lower bound hold for a finite domain pdf function (x squared)", {
  g1 <- function (x) {
    return(x^2)
  }
  h1 <- function (x) {
    return(log(g(x)))
  }
  
  init_abscissa(h1, 1, 2)
  data1 <- init_data(g1, 1, 2)
  
  test_x <- seq(1, 2, length.out=50)
  test_hx <- sapply(test_x, h1)
  test_uk <- sapply(test_x, uk, data=data1)
  test_lk <- sapply(test_x, lk, data=data1)
  
  expect_true(all(test_hx <= test_uk))
  expect_true(all(test_hx >= test_lk))
})

test_that("The upper bound and lower bound hold in a infinite domain pdf function (normal)", {
  g2 <- function (x) {
    return(dnorm(x))
  }
  
  h2 <- function(x) {
    return(log(g2(x)))
  }
  
  init_abscissa(h2, -Inf, Inf)
  data2 <- init_data(g2, -Inf, Inf)
  
  test_x <- seq(-10, 10, length.out=50)
  test_hx <- sapply(test_x, h2)
  test_uk <- sapply(test_x, uk, data=data2)
  test_lk <- sapply(test_x, lk, data=data2)
  test_hx <= test_uk
  test_hx >= test_lk
  
  expect_true(all(test_hx <= test_uk))
  expect_true(all(test_hx >= test_lk))
})

# uk(1.1, data)
# lk(1.1, data)

# a small test on x^2, finite bound
# test_x <- seq(1, 2, length.out=5)
# test_hx <- sapply(test_x, h1)
# test_uk <- sapply(test_x, uk, data=data1)
# test_lk <- sapply(test_x, lk, data=data1)
# test_hx <= test_uk
# test_hx >= test_lk

# plot_x <- seq(1, 2, length.out=10000)
# plot_hx <- sapply(plot_x, h1)
# plot(plot_x, plot_hx, type='l')
# lines(test_x,  test_uk, type='l', col='blue')
# lines(test_x,  test_lk, type='l', col='green')
# plot_x <- seq(-10, 10, length.out=10000)
# plot_hx <- sapply(plot_x, h2)
# plot(plot_x, plot_hx, type='l')
# lines(test_x,  test_uk, type='l', col='blue')
# lines(test_x,  test_lk, type='l', col='green')