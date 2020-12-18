library(testthat)
source('../utils.R')
# 
# test_that("Dummpy", {
#   expect_true(TRUE)
# })

test_that("The upper bound and lower bound hold for a finite domain pdf function (x squared)", {
  g1 <- function (x) {
    return(x^2)
  }
  h1 <- function (x) {
    return(log(g1(x)))
  }

  init_abscissa(h1, 1, 2)
  data1 <- initialization_step(h1, 1, 2)

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
  data2 <- initialization_step(h2, -Inf, Inf)

  test_x <- seq(-10, 10, length.out=50)
  test_hx <- sapply(test_x, h2)
  test_uk <- sapply(test_x, uk, data=data2)
  test_lk <- sapply(test_x, lk, data=data2)
  test_hx <= test_uk
  test_hx >= test_lk

  expect_true(all(test_hx <= test_uk))
  expect_true(all(test_hx >= test_lk))
})