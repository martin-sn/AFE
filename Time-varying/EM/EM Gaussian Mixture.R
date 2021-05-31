### Generate Gaussian Mixture ###

## Mixture in matrix format. First column is mean, second column is standard deviation, third is mix proportion. ##
Mixture = matrix(c(5,-4,1.5,0.7,0.95,0.05), nrow = 2, byrow = FALSE)


## Simulate Gaussian Mixture ##


SimGaussMix <- function(n, Mixture){
  
  Y = array(NA, n)
  J = nrow(Mixture)
  
  for (i in 1:n){
    Gaus = sample(c(1:J),size = 1, prob = Mixture[,3])
    Y[i] = rnorm(1, mean = Mixture[Gaus,1], sd = Mixture[Gaus,2])

  }
  return(Y)
}

Y = SimGaussMix(1000, Mixture)

plot(Y, type = "l")

### Estimate Gaussian Mixture by EM Algorithm ###

# E-Step
E_EMGaussMix <- function(Y,J, Mixture){
  
  N = length(Y)
  P = array(NA, c(N,J))
  
  for (i in 1:N){
    for (j in 1:J){
      P[i,j] = Mixture[j,3]*dnorm(Y[i], Mixture[j,1], Mixture[j,2]) / sum(Mixture[,3]*dnorm(Y[i], Mixture[,1], Mixture[,2]))
    }
  }
  return(P)
}


# M-Step 
M_EMGaussMix <- function(Y,J,P){
  
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


Run_EMGaussMix(Y,2,1000)



plot(density(Y))


##### Rigorous Edition #####


EStep <- function(Y,J,Mixture){
  
  N = nrow(Y)
  
  # Forward Probabilities #
  
  a = array(NA, dim = c(length(Y), J))
  
  for (i in 1:N){
    for (j in 1:J){
      a[i,j] = dnorm(Y[i], Mixture[j,1], Mixture[j,2], log = FALSE) * (Mixture[j,3])
    }
  }
  
  # Backwards Probabilities # 
  
  # Forward and backward probabilities are equivalent in a gaussian mixture! 
  
  
  # Smoothed Probabilities # 
  
  SP <- array(NA, dim = dim(a))
  
  for (i in 1:N){
    for (j in 1:J){
      SP[i,j] = a[i,j] /sum(a[i,])
    }
  }
  return(SP)
}

test <- EStep(Y, 2, Mixture)

MStep <- function(Y, SP, J){
  
  Mixture = matrix(NA, nrow = J, ncol=3)
  
  N = length(Y)
  
  for(j in 1:J){
    
    # Probability estimate
    Mixture[j,3] = sum(SP[,j]) / N
    
    # Mean estimate
    Mixture[j,1] = sum(SP[,j]*Y)/sum(SP[,j])
    
    # Variance estimate
    Mixture[j,2] = sqrt(sum(SP[,j]*(Y-Mixture[j,1])**2) / sum(SP[,j]))

  }
  return(Mixture)
}


Mix <- MStep(Y, test, 2)

Mix

RunEM <- function(Y, Iterations, J){
  
  Mixture <- matrix(c(3,-2,0.3,1.5,0.3,0.7), nrow = 2, byrow = FALSE)

    
  for (i in 1:Iterations){
    SP = EStep(Y,J,Mixture)
    
  
    Mixture = MStep(Y,SP,J)
    
    }
 
  return(Mixture) 
}

res <- RunEM(Y,1000,2)

res


