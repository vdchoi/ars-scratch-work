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

# g <- function (x) {
#   return(dnorm(x))
# }
# 
# ars_results <- adaptive_rejection_sampling(g, a = -Inf, 
#                                            b= Inf, N = 100000, n_per_step = 200)
# hist(ars_results, breaks=30)
# print(c(length(ars_results), mean(ars_results), sd(ars_results)))

  ###########
# 
# g <- function (x) {
#   return(x^2)
# }
# 
# ars_results <- adaptive_rejection_sampling(g, a = 1,
#                                            b= 2, N = 10, n_per_step = 200)
# hist(ars_results, breaks=30)
# print(c(mean(ars_results), sd(ars_results)))


# h2 <- function(x) {
#   return(log(g2(x)))
# }
# 
# init_abscissa(h2, -Inf, Inf)
# data2 <- init_data(g2, -Inf, Inf)
# exp_sampling(100, data2, h2)
# hist(rejection_step(exp_sampling(10000, data2, h2), data2, h2)$sample, breaks = 30)



