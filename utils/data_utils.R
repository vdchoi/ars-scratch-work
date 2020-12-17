library(numDeriv)
# f is the function, a and b are its range
# assuming a, b are real numbers for now
init_data <- function (h, a, b) {
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

init_abscissa <- function(h, a, b, init_k = 20) {
  
  # both finite
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
  # x_value in [z_{j-1}, z_{j}]
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

lk <- function (x_value, data = data) {
  # x_value in [x_{j}, x{j+1}]
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

update_data <- function (data, new_data) {
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

