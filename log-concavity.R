library(numDeriv)
library(ramify)
library(assertthat)

verify_log_concavity <- function(logf, xmin, xmax, npoints) {
  #' Verify that a given function is concave
  #' (Input)
  #' logf: The function. It must take only one argument.
  #' The log should have already been applied to the function.
  #' xmin, xmax: The domain of the function
  #' npoints: The number of samples that will be taken in the domain
  #' (Output)
  #' 0 if the function is log-concave, 1 otherwise
  
  # Assertions
  # Assert that logf is a function taking 1 argument
  assert_that(class(f) == "function")
  assert_that(length(args(f)) == 1)
  # Assert that xmin, xmax, npoints are scalars
  assert_that(is.numeric(xmin))
  assert_that(length(xmin) == 1)
  assert_that(is.numeric(xmax))
  assert_that(length(xmax) == 1)
  assert_that(is.numeric(npoints))
  assert_that(length(npoints) == 1)
  # Verify that xmin < xmax
  assert_that(xmin < xmax)
  # Verify that npoints > 1
  assert_that(npoints > 1)
  
  # Define the appropriate range of x values
  x_vec <- linspace(xmin, xmax, npoints)
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


# Example function that is log-concave
# Normal distribution with mean 20.00 and stdev of 5.00
f <- function(x) {
  return(log(1/(sqrt(2.00 * pi)*5.00)*exp(-0.50*((x-20.0)/5.00)^2)))
}

verify_log_concavity(f, -1.00e2, 1.00e2, 1000)

# Example function that is not log-concave
g <- function(x) {
  return(log(4.00^x))
}

verify_log_concavity(g, 0.00, 1.00, 100)


