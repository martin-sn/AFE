##### EM Algorithm for Zero Inflated Poisson #####


SimZIP <- function(n, p, lambda){
  Y = array(NA, n)
  for (i in 1:n){
    Y[i] = (1-rbinom(1,1, prob = p))*rpois(1,lambda)
  }
  return(Y)
}



Y = SimZIP(1000, 0.15, 15)
plot(Y, type = "l")


### Estimate ZIP by EM Algorithm ###

# E-Step
PMF_ZIP <- function(y,p,lambda){
  if (y == 0){
    prob = p + (1-p)*exp(-lambda)
  }
  else{
    prob = (1-p)*(exp(-lambda)*lambda**(y))/factorial(y)
  }
  return(prob)
}


E_EMZIP <- function(Y,p,lambda){
  
  N = length(Y)
  P = array(NA, c(N,2))
  
  for (i in 1:N){
    P[i,1] <- ifelse(Y[i] != 0, 0, p)
    P[i,2] <- ifelse(Y[i] != 0, (1-p)*(exp(-lambda)*lambda**(Y[i])/factorial(Y[i])), (1-p)*exp(-lambda))
  }
  
  SP <- array(NA, c(N,2))
  
  for (i in 1:N){
    for (j in 1:2){
      SP[i,j] = P[i,j]/sum(P[i,])
    }
  }
  return(SP)

}
P <- E_EMZIP(Y, 4, 1)

#### Probability of a ZIP ####
# Y can be 0, if either poison process outputs zero, or the bernoulli is 1
# So it is the probability that B = 1 + prob that B != 0 and pois = 0
# If y is not zero, then the probability that B is 1 is 0. 


# M-Step 
M_EMZIP <- function(Y,J,P){
  
  N = length(Y)
  
  p = sum(P[,1]) /sum(P)
  
  lambda = sum(P[,2]*Y) / sum(P[,2])
    
  out = matrix(c(p,lambda))
  return(out)
}


RUN_EMZIP <- function(Y,J,Iterations){

  # Init
  p <- runif(1,0,1)
  lambda = mean(Y)
  
  
  for (i in 1:Iterations){
    P <- E_EMZIP(Y,p,lambda)
    M <- M_EMZIP(Y,J,P)
    p = M[1]
    lambda = M[2]
    
  }
  
  print(p)
  print(lambda)
  
}

RUN_EMZIP(Y,2,100)

library(tidyverse)

JMP <- read_csv("Time-varying/Data/JPM.csv", col_names = FALSE)


RUN_EMZIP(abs(JMP$X2), 2, 1000)



###### Mixture ZIP ####### 

# Mixture format

#Matrix

# Col 1 = Switch prob
# Col 2 = lambda

SimMixtureZIP <- function(n,mixture, inflationprob){
  J = nrow(mixture)
  Y = array(NA, n)
  
  for (i in 1:n){
    j = sample(1:J, size = 1, prob=mixture[,1])
    Y[i] = (1-rbinom(1,1, prob = inflationprob))*rpois(1,mixture[j,2])
  }
  return(Y)
}


SimMix <- matrix(c(0.75,0.25, 30, 5), nrow = 2, byrow = FALSE)

Y <- SimMixtureZIP(1000, SimMix, 0.65)

plot(Y, type = "l")



