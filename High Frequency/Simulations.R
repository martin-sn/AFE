
SimY <- function(n, sigma, omega){
  Mu = 3
  dt = 1/n 
  dW = rnorm(n)*sqrt(dt)
  #X = cumsum(sigma*dW + Mu/n)
  
  X = array(NA, dim = n)
  
  X[1] = sigma*dW[1] + Mu/n
  for (i in 2:n){
      X[i] = ifelse(i == 1523, X[i-1] + sigma*dW[i] + Mu/n + 1, X[i] = X[i-1] sigma*dW[i] + Mu/n)
  }
  return(X)
}

Y_omega01 = SimY(10000, 1,0)
plot(Y_omega01)

?array



SimY <- function(n, sigma){
  Mu = 3
  dt = 1/n 
  dW = rnorm(n)*sqrt(dt)

  
  X = array(NA, dim = n)
  
  X[1] = sigma*dW[1] + Mu/n

  
  for (i in 2:n){
    X[i] = ifelse(i == 1523, X[i-1] + sigma*dW[i] + Mu/n + 2, X[i-1] + sigma*dW[i] + Mu/n)
  }
  
  return(X)
}

SimY <- function(n, sigma){
  Mu = 3
  dt = 1/n 
  dW = rnorm(n)*sqrt(dt)
  
  
  X = array(NA, dim = n)
  
  X[1] = sigma*dW[1] + Mu/n
  
  
  for (i in 2:n){
    X[i] = ifelse(i == 1523, X[i-1] + sigma*dW[i] + 1, X[i-1] + sigma*dW[i])
  }
  
  return(X)
}



SimY <- function(n, sigma){
  Mu = 3
  dt = 1/n 
  dW = rnorm(n)*sqrt(dt)
  
  
  X = array(NA, dim = n)
  
  X[1] = sigma*dW[1] + Mu/n
  
  
  for (i in 2:n){
    X[i] = ifelse(i == 1523, X[i-1] + sigma*dW[i], X[i-1] + sigma*dW[i])
  }
  
  return(X)
}


Y = SimY(2500000,2)

plot(Y)

RV <- function(Y){
  RV = 0
  for (i in 2:length(Y)){
    RV = RV + (Y[i] - Y[i-1])**2
  }
  return(RV)
}

RV(Y)


2**2 + 3/10000*(3**2) + 1/10000


