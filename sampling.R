setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('utils/data_utils.R')

exp_sampling = function(n, x_list) {
  
  # Grab the number of rows
  k = x_list$k
  
  # Save the intersection points of our upper hull
  intersections = x_list$z
  
  # First sample from a multinomial with probabilities corresponding to 
  # the area under each piece of the adjusted piecewise upper hull
  
  # Get the area under the curve for each section
  # Maybe we should switch to the integrate function here
  upper = x_list$hx + ((intersections[2:(k+1)] - x_list$x)*x_list$dhx)
  probabilities = abs(diff(exp(upper)/x_list$dhx))
  sum_prob = sum(probabilities)
  
  prob_vector = probabilities / sum_prob
  
  
  # Sample and return a vector of indexes
  sample_multinom = rmultinom(n, 1, prob_vector)
  i_vector = colSums(sample_multinom*(1:(k-1)))
  
  # Next let's draw n uniform variables
  unif_draws = runif(n)
  
  slope_j = x_list$dhx[i_vector]
  z_j = intersections[i_vector]
  constant_j = slope_j*z_j + upper[i_vector]
  prob_indexed = prob_vector[i_vector]
  
  
  sample = log(exp(slope_j*z_j) + 
                 sum_prob*slope_j*prob_indexed*exp(-constant_j)*unif_draws) / slope_j
  
  return(sample)
}

rejection_step = function(sample_vector, x_list, f_x) {
  
  # Find upper and lower shell functions
  
  lower_shell = sapply(sample_vector, lk, data=x_list)
  upper_shell = sapply(sample_vector, uk, data=x_list)
  
  # Uniform Sample
  unif_w = runif(n = length(sample_vector))
  
  # Squeezing step
  
  log_vec = unif_w <= exp(lower_shell - upper_shell)
  
  # Rejection step
  rejection_step_log_vec = unif_w[!log_vec] <= 
    exp(log(f_x(sample_vector[!log_vec])) - 
          upper_shell[!log_vec])
  
  log_vec[which(!log_vec)] = rejection_step_log_vec
  
  # Return accepted values
  return(sample_vector[log_vec])
  
}

set.seed(1)
x = c(-0.2, 0.1, 0.2)
f_x = log(dnorm(x))

func_temp = function(x){
  return(log(dnorm(x)))
}

library(numDeriv)

# Problem - we're seeing gradient always being - f_x, is that correct?
f_p_x = grad(func_temp, x)

normal_test = matrix(c(x, f_x, f_p_x), nrow = length(x))
exp_sampling(10, normal_test)

