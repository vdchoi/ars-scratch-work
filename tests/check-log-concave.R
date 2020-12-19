library(numDeriv)
g <- function(x){
  return(x*sin(x))
}


mean <- integrate(g, lower = 0, upper = pi)[[1]]


g2 <- function(x){
  return((x-mean)^2*sin(x))
}

variance <- integrate(g2, lower = 0, upper = pi)[[1]]
