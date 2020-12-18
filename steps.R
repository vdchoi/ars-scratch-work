
# f is the function, a and b are its range
# assuming a, b are real numbers for now
initialization_step <- function (h, a, b) {
  # abscissa, arrays of data, excluding uppper and lower bound a, b
  
  x <- init_abscissa(h, a, b)
  
  hx <- sapply(x, h)
  dhx <- grad(h, x)
  
  k <- length(x)
  
  # calculate zs' using formula (1)
  # 2:k -- j+1,  1:k-1 -- j
  z <- (hx[2:k] - hx[1:k-1] - x[2:k] * dhx[2:k] + x[1:k-1] * dhx[1:k-1]) / (dhx[1:k-1] - dhx[2:k])
  
  # z's index start from 0, ends with k, i.e., z[1] represents z_0, z[k+1] represents z_k
  z <- c(a, z, b)
  data <- list(x = x, hx = hx, dhx = dhx, z = z, k = length(x), a = a, b = b)
  return(data) 
}


sampling_step <- function(sample_vector, data, h) {
  
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

