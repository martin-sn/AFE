##### EM Algorithm for Zero Inflated Poisson #####


SimZIP <- function(n, p, lambda){
  Y = array(NA, n)
  for (i in 1:n){
    Y[i] = (1-rbinom(1,1, prob = p))*rpois(1,1)
  }
  return(Y)
}

plot(SimZIP(100,0.5,10), type = "l")
rpois(1,10)


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


E_EMZIP <- function(Y,J, Mixture){
  
  N = length(Y)
  P = array(NA, c(N,J))
  
  for (i in 1:N){
    for (j in 1:J){
      P[i,j] = dnorm(Y[i], Mixture[j,1], Mixture[j,2])
    }
  }
  return(P)
}

#### Probability of a ZIP ####
# Y can be 0, if either poison process outputs zero, or the bernoulli is 1
# So it is the probability that B = 1 + prob that B != 0 and pois = 0
# If y is not zero, then the probability that B is 1 is 0. 


# M-Step 
M_EMZIP <- function(Y,J,P){
  
  N = length(Y)
  
  Mixture = array(NA, dim = c(J,3))
  
  for (j in 1:J){
    Mixture[j,1] = 1/sum(P[,j]) * sum(P[,j]*Y)
    Mixture[j,2] = sqrt(1/sum(P[,j]) * sum(P[,j]*(Y-Mixture[j,1])**2))
    Mixture[j,3] = sum(P[,j]) / N 
  }
  return(Mixture)
}



Run_EMGaussMix <- function(Y,J, IT){
  
  # Initialize mixture
  
  N = length(Y)
  
  Mixture = array(NA, dim = c(J,3))
  
  seq_init = seq(from = 0.7, to = 1.5, length.out = J)
  
  Mixture[,1] = mean(Y) * seq_init
  Mixture[,2] = sd(Y) * seq_init
  Mixture[,3] = 1/J
  
  for (i in 1:IT){
    P = E_EMGaussMix(Y,J, Mixture)
    Mixture = M_EMGaussMix(Y,J,P)
  }
  
  LLK = 0
  
  for (i in 1:N){
    K_LLK = 0
    for (j in 1:J){
      K_LLK = K_LLK + Mixture[j,3]*dnorm(Y[i], mean = Mixture[j,1], sd = Mixture[j,2])
    }
    
    LLK = LLK + log(K_LLK)
    
    
    
  }
  
  print(paste("LLK:",sum(LLK)))
  
  BIC = (nrow(Mixture)*3-1)*log(N) - 2*sum(LLK)
  
  print(paste("BIC:", BIC))
  return(Mixture)
  
  
}