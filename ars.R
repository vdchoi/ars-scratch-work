# setting working directory to the file location
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source('data_utils.R')
# source('steps.R')

#' Adaptive Rejection Sampler
#'
#' Sample from an arbitrary log-concave function. 
#' 
#'
#' @param g input density function
#' @param a lower bound of function domain D (could be -Inf)
#' @param b upper bound of function domain D (could be Inf)
#' @param N number of observations required 
#' @param n_per_step Default = 100. Number of observations between
#' each update of the hull function
#' @return N independent observations from the the given input density
#' @examples 
#' g <- function (x) {
#' return(dnorm(x))
#' }

#' ars_results <- ars(g, a = -Inf, b= Inf, N = 100000)
#' 
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


