library(numDeriv)
library(docstring)

init_abscissa <- function(h, a, b, init_k = 20) {
  #' Initialize Abscissas
  #'
  #' Initialize the abscissas (x vector), based on the function domain D. 
  #' If D is bounded, select xs with equal distances, excluding the bounds
  #' If D is left/right unbounded, select first/last x based on h's derivative values
  #'
  #' @param h log density function
  #' @param a lower bound of function domain D (could be -Inf)
  #' @param b upper bound of function domain D (could be Inf)
  #' @param init_k number of xs in the initial abscissa, default is 20
  #' @return a vector consisting of all xs in the abscissa

  # both left and right finite
  if (a != -Inf && b != Inf){
    return(seq(a, b, length.out = init_k+2)[2:(init_k+1)])
  }
  
  # if infinite on left side, find x_0, s.t. dhx_0 > 0
  if (a == -Inf){
    x_0 <- -10
    while (grad(h, x_0) <= 0){
      x_0 <- x_0 * 2
    }
  }
  # else set a finite left bound
  else {
    x_0 = a + 1
  }
  
  # similarly if infinte on right side
  if (b == Inf){
    x_k <- 10
    while (grad(h, x_k) >= 0){
      x_k <- x_k * 2
    }
  }
  else {
    x_k = b - 1
  }
  
  x <- seq(x_0, x_k, length.out = init_k)
  
  return(x)
}

# x_value
uk <- function (x_value, data) {
  #' Upper envelope function for h (u_k in formula (2))
  #'
  #' Calculating the upper envelope function at point x_value, 
  #' given the overall maintained data structure
  #'
  #' @param x_value the position to evaluate the u_k function
  #' @param data the maintained data structure for the overall process
  #' @return the value of the upper envelope for u_k(x_value)
  
  # get index j x_value in [z_{j-1}, z_{j}]
  j <- sum(x_value > data$z) 
  # boundary cases
  if (x_value == data$a){
    j = 1
  }
  if (x_value == data$b){
    j = data$k
  }
  
  # calculating with equation (2)
  result <- data$hx[j] + (x_value - data$x[j]) * data$dhx[j]
  
  return(result)
}

lk <- function (x_value, data) {
  #' Lower envelope function for h (l_k in formula (4))
  #'
  #' Calculating the lower envelope function at point x_value, 
  #' given the overall maintained data structure
  #'
  #' @param x_value the position to evaluate the u_k function
  #' @param data the maintained data structure for the overall process 
  #' @return the value of the lower envelope for l_k(x_value)
  
  # get index j for x_value in [x_{j}, x{j+1}]
  j <- sum(x_value > data$x)
  
  if (j == 0 || j == data$k) {
    result <- -Inf
  }
  else {
    # calculating with equation (4)
    result <- ((data$x[j+1] - x_value) * data$hx[j] + (x_value - data$x[j]) * data$hx[j+1] )/
      (data$x[j+1] - data$x[j])
  }
  return(result)
}

exp_sampling <- function(n, data, h) {
  #' Sample from a Piece-wise Exponential function
  #'
  #' Generate a sample of length n from the current data and the function h
  #'
  #' @param n number of samples
  #' @param data data structure output from the initialization_step function
  #' @param h log density function
  #' @return a vector of length n of samples from a piece-wise exponential corresponding to function h
  
  # Save intermediate values for k, the intersection points, and the upper hull
  k <- data$k
  intersections <- data$z
  upper <- sapply(data$x, uk, data=data)
  
  # Create a vector of probabilities relative to the area between intersections
  probabilities <- abs(diff(exp(h(intersections)))/(data$dhx))
  sum_prob <- sum(probabilities)
  prob_vector <- probabilities / sum_prob
  
  # Sample and return a vector of indexes
  sample_multinom <- rmultinom(n, 1, prob_vector)
  i_vector <- colSums(sample_multinom*(1:k))
  
  # Draw n uniform samples
  unif_draws <- runif(n)
  
  # Use inverse transform sampling to transform uniform samples
  slope_j <- data$dhx[i_vector]
  z_j <- intersections[i_vector]
  z_j_1 <- intersections[i_vector+1]
  constant_j <- slope_j*z_j + upper[i_vector]
  ac_pi <- (exp(z_j_1*slope_j) - exp(z_j*slope_j)) * exp(constant_j)
  
  sample <- log(exp(slope_j*z_j) + 
                  ac_pi*exp(-constant_j)*unif_draws) / slope_j
  
  return(sample)
}



