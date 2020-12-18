# setting working directory to the file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source('data_utils.R')
# source('steps.R')

ars <- function (g, a, b, N, n_per_step = 100){
  # vector for total results
  total_sample_results <- c()
  # define the h function
  h <- function (x) {
    return(log(g(x)))
  }

  # function checks on g
  
  
  # initialization step
  data <- initialization_step(h, a, b)
  
  while (length(total_sample_results) < N) {
    
    # sampling from the piecewise exponential distribution
    exp_sample_results <- exp_sampling(n_per_step, data, h)
    
    # rejection, and get current step results
    current_step_sample_results <- sampling_step(exp_sample_results, data, h)
    
    # concatenate into sampled results
    total_sample_results <- c(total_sample_results, current_step_sample_results$sample)
    
    # update the data structure
    data <- update_step(data, current_step_sample_results)
  }
  
  return(total_sample_results[1:N])
}


