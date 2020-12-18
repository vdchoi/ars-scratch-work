setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(assertthat)
library(testthat)
source('../utils.R')
source('../steps.R')
source('../asr.R')


g <- function (x) {
  return(dnorm(x))
}

ars_results <- adaptive_rejection_sampling(g, a = -Inf, 
                                           b= Inf, N = 100000, n_per_step = 200)
hist(ars_results, breaks=30)
print(c(length(ars_results), mean(ars_results), sd(ars_results)))

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