#### Question 1d ####
# Series with noise #

SimY <- function(n, sigma, omega){
  U = rnorm(n,0,Omega)
  dt = 1/n 
  dW = rnorm(n)*sqrt(dt)
  X = cumsum(sigma*dW)
  Y = X+U
  
  return(Y)
}

