library(docstring)
library(assertthat)

initialization_step <- function (h, a, b) {
  #' Initialize Step
  #'
  #' Initialize the data structure we maintain given the log density function h as input, 
  #' as described in 2.2.1. 
  #' We calculate the initial abscissa, h(x) vector, h'(x) vector, 
  #' and z vector needed for further computation
  #'
  #' @param h log density function
  #' @param a lower bound of function domain D (could be -Inf)
  #' @param b upper bound of function domain D (could be Inf)
  #' @return a list, the main data structure we will work with, containing
  #' x, h(x), d'(x), z, and other information
  
  # get the initial abscissa vector 
  x <- init_abscissa(h, a, b)
  
  # calculated needed information, h(x) and h'(x)
  hx <- sapply(x, h)
  dhx <- grad(h, x)
  
  k <- length(x)
  
  # calculate z vector using formula (1)
  # 2:k -- indexing with (j+1),  1:k-1 -- indexing with j
  z <- (hx[2:k] - hx[1:k-1] - x[2:k] * dhx[2:k] + x[1:k-1] * dhx[1:k-1]) / (dhx[1:k-1] - dhx[2:k])
  
  # z's index start from 0, ends with k, i.e., z[1] represents z_0, z[k+1] represents z_k
  z <- c(a, z, b)
  
  # construct the main data structure we work with in this project
  data <- list(x = x, hx = hx, dhx = dhx, z = z, k = length(x), a = a, b = b)
  return(data) 
}


sampling_step <- function(sample_vector, data, h) {
  #' Rejection Sampling Step
  #'
  #' Given a sample output from the exp_sampling function, perform rejection
  #' squeezing and rejection sampling steps. Output will be a structured list
  #' of data to update
  #'
  #' @param sample_vector vector of x-values. Output of exp_sampling
  #' @param data data structure output from the initialization_step function
  #' @param h log density function
  #' @return output of sample data in our data structure
  
  # Find upper and lower shell functions
  lower_shell <- sapply(sample_vector, lk, data=data)
  upper_shell <- sapply(sample_vector, uk, data=data)
  
  # Uniform Sample
  unif_w <- runif(n = length(sample_vector))
  
  # Squeezing step
  log_vec <- unif_w <= exp(lower_shell - upper_shell)
  
  # Rejection step
  eval_x <- sample_vector[!log_vec]
  eval_hx <- h(eval_x)
  eval_dhx <- grad(h, eval_x)

  rejection_step_log_vec <- unif_w[!log_vec] <= 
    exp(eval_hx - upper_shell[!log_vec])
  
  log_vec[which(!log_vec)] <- rejection_step_log_vec
  
  # Return accepted values
  return(list(sample = sample_vector[log_vec], 
              x = eval_x,
              hx = eval_hx,
              dhx = eval_dhx))
}

update_step <- function (data, new_data) {
  #' Update Step
  #'
  #' Combine old data with newly computed data, 
  #' giving a new set of abscissa, and relevant information,
  #' as described in 2.2.3. 
  #' We combine the x, h(x), h'(x) values, and sort them based on x values,
  #' then calculate the new z vector, and return updated data structure
  #'
  #' @param data a list, the main data structure
  #' @param new_data a list, the data structure with similar members 
  #' but with newly calculated x, h(x) and h'(x)
  #' @return a list, the updated main data structure
  
  # checking inputs
  assert_that(is.list(data))
  assert_that(is.list(new_data))

  assert_that(is.vector(data$x))
  assert_that(is.vector(data$hx))
  assert_that(is.vector(data$dhx))

  assert_that(is.vector(new_data$x))
  assert_that(is.vector(new_data$hx))
  assert_that(is.vector(new_data$dhx))
  
  # combining data and sorting
  combined_x <- c(data$x, new_data$x)
  combined_hx <- c(data$hx, new_data$hx)
  combined_dhx <- c(data$dhx, new_data$dhx)
  
  sort_order <- order(combined_x)
  
  # new arrays
  x <- combined_x[sort_order]
  hx <- combined_hx[sort_order]
  dhx <- combined_dhx[sort_order]
  
  k <- length(x)
  data$k <- k
  
  # recalculate z after sorting
  z <- (hx[2:k] - hx[1:k-1] - x[2:k] * dhx[2:k] + x[1:k-1] * dhx[1:k-1]) / (dhx[1:k-1] - dhx[2:k])
  z <- c(data$a, z, data$b)
  
  data <- list(x = x, hx = hx, dhx = dhx, z = z, k = k, a = data$a, b = data$b)
  return(data)
}

