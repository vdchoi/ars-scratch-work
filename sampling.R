setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('utils/data_utils.R')

exp_sampling <- function(n, data, function_temp) {
  
  # Grab the number of rows
  k <- data$k
  
  # Save the intersection points of our upper hull
  intersections <- data$z
  
  # First sample from a multinomial with probabilities corresponding to 
  # the area under each piece of the adjusted piecewise upper hull
  
  # Get the area under the curve for each section
  # Maybe we should switch to the integrate function here
  upper <- sapply(data$x, uk, data=data)
  
  # To check - do we add the function into our data structure
  
  probabilities <- abs(diff(exp(function_temp(intersections)))/(data$dhx))
  sum_prob <- sum(probabilities)
  prob_vector <- probabilities / sum_prob
  
  # Sample and return a vector of indexes
  sample_multinom <- rmultinom(n, 1, prob_vector)
  
  i_vector <- colSums(sample_multinom*(1:k))
  
  # Next let's draw n uniform variables
  unif_draws <- runif(n)
  slope_j <- data$dhx[i_vector]
  z_j <- intersections[i_vector]
  z_j_1 <- intersections[i_vector+1]
  constant_j <- slope_j*z_j + upper[i_vector]
  
  ac_pi <- (exp(z_j_1*slope_j) - exp(z_j*slope_j)) * exp(constant_j)
  
  
  sample <- log(exp(slope_j*z_j) + 
                 ac_pi*exp(-constant_j)*unif_draws) / slope_j
  
  return(sample)
}

rejection_step <- function(sample_vector, data, h) {
  
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
  return(list(sample <- sample_vector[log_vec], 
              x = eval_x,
              hx = eval_hx,
              dhx = eval_dhx))
}

adaptive_rejection_sampling <- function (g, a, b, N, n_per_step = 100){
  # vector for total results
  total_sample_results <- c()
  # define the h function
  h <- function (x) {
    return(log(g(x)))
  }

  # check log-concavity of g
  
  # initialization step
  data <- init_data(h, a, b)
  
  while (length(total_sample_results) < N) {
  
    # sampling from the piecewise exponential distribution
    exp_sample_results <- exp_sampling(n_per_step, data, h)
    
    # rejection, and get current step results
    current_step_sample_results <- rejection_step(exp_sample_results, data, h)
    
    # concatenate into sampled results
    total_sample_results <- c(total_sample_results, current_step_sample_results$sample)
    
    # update the data structure
    data <- update_data(data, current_step_sample_results)
  }
  
  return(total_sample_results)
}

g <- function (x) {
  return(dnorm(x))
}

ars_results <- adaptive_rejection_sampling(g, a = -Inf, 
                                           b= Inf, N = 10000, n_per_step = 200)
hist(ars_results, breaks=30)
print(c(mean(ars_results), sd(ars_results)))

###########
# 
g <- function (x) {
  return(x^2)
}

ars_results <- adaptive_rejection_sampling(g, a = 1,
                                           b= 2, N = 10000, n_per_step = 200)
hist(ars_results, breaks=30)
print(c(mean(ars_results), sd(ars_results)))




# h2 <- function(x) {
#   return(log(g2(x)))
# }
# 
# init_abscissa(h2, -Inf, Inf)
# data2 <- init_data(g2, -Inf, Inf)
# exp_sampling(100, data2, h2)
# hist(rejection_step(exp_sampling(10000, data2, h2), data2, h2)$sample, breaks = 30)



