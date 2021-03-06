###### Expectation Maximization for Hidden Markov Models #####

# Get functions form Utils.R 

source("Time-varying/HMM/Utils.R")

### Simulate T = 5000 observations from a HMM with J = 2 regimes and transition probability matrix

dat = c(0.95, 0.05, 0.03, 0.97)
p = matrix(dat,nrow = 2, byrow = TRUE)

sim = SimulateHMM(p, 5000, 2)

plot(sim, type = "l")

GaussianMix = matrix(c(5,-4,1.5,0.7), nrow = 2)

y <- GenerateY_GausHMM(sim,GaussianMix)

plot(y, type = "l")


######## EM Algorithm for HMM ##########

# EM is designed to compute the maximum likelihood estimator in the case of missing values
# We consider the unobserved states of the Markov chain as missing data
# 1. Choose a starting value for theta
# 2. E-Step: Compute the conditional expectations of those functions 
# of the missing data appear in the complete-data log-likelihood
# 3. M-step: Maximization of the log-likelihood with respect to the est of parameters
# to be estimated (the missing data are substituted by their conditional expectation)
# 4. Iterate 2. and 3. until convergence. 


EM <- function(Y, J, p, delta, mix){
  
  n = length(Y)
  
  # Init forward probability array
  alpha <- array(NA, c(J, n))

  # First forward probability
  for (j in 1:J){
    alpha[j,1] = delta[1,j]*dnorm(y[1], mix[j,1], mix[j,2], log = FALSE)
  }
  # Remaining forward probabilities
  
  # P(y_1) is a diagonal matrix with j,j-th element p_j(y_1)
  # p_j(y_1) is the density for gaussian j
  

  P = array(0, dim = c(J,J,n))
  
  ### Compute conditional probabilities ###
  
  for (i in 1:n){
    for (j in 1:J){
      P[j,j,i] = dnorm(y[i], mix[j,1],mix[j,2])
    }
  }

  for (i in 2:length(Y)){
    alpha[,i] = alpha[,i-1]%*%p%*%P[,,i]
  }
  
  # Backwards probabilities
  beta = array(NA, dim = c(J, length(Y)))
  
  beta[,length(Y)] = matrix(c(1/n,1/n)) # Initialization. Not sure about this one. 
  
  for (i in (length(Y)-1):1){
    beta[,i] = p%*%P[,,i]%*%beta[,i+1]
  }
  
  out = list()
  
  out[["alpha"]] = alpha
  out[["beta"]] = beta
  
  
  # Multiplying the forward and backward probabilities gives us the likelihood
  
  LK = matrix(NA, nrow = 2)
  
  #LK[1] = sum(alpha[1,]*beta[1,])
  #LK[2] = sum(alpha[2,]*beta[2,])
  
  LK[1] = sum(log(alpha[1,]) + log(beta[1,]))
  LK[2] = sum(log(alpha[2,]) + log(beta[2,]))
  
  #out[["LLK"]] = log(LK)
  
  out[["LLK"]] = log(LK)
  
  
  #### Smoothed probabilities ####
  
  smooth = array(NA, dim = c(J, length(Y)))
  
  for (i in 1:n){
    for (j in 1:J){
      smooth[j,i] = alpha[j,i]*beta[j,i]/(t(alpha[,i])%*%beta[,i])
    }
    
  }
  
  out[["Smooth"]] = smooth
  
  # Smoothed probabilities of switching
  
  # We have (T-1)*J^2  of the switching probabilities
  
  SmoothSwitch = array(NA, dim = c(J,J,n)) # Not sure about this 
  
  SmoothSwitch[,,1] = p
  
  for (i in 2:n){
    for (j in 1:J){
      for (k in 1:J){
        SmoothSwitch[k,j,i] = alpha[j,i-1] * p[j,k] * P[k,k,i] * beta[k,i] / (t(alpha[,i])%*%beta[,i])
      }
    }
  }
  
  
  out[["SmoothSwitch"]] = SmoothSwitch
    
  
  #### Q Func ####
  
  # It seems that we actually do not need to compute the Q func. 
  
  # P(u) is equivalent to smooth prob
  # P(v) is equivalent to smooth switch
  
  
  # Q_func = 0 
  # 
  # E_J = 0
  # E_JJ = 0
  # E_JJJ = 0
  # 
  # for (jj in 1:J){
  #   E_J = E_J + sum(smooth[1,jj]*log(delta[jj]))
  #   for (k in 1:J){
  #     E_JJ = E_JJ + sum(SmoothSwitch[jj,k,]) * log(p[jj,k])
  #   }
  # }
  #   for (iii in 2:n){
  #     for (jjj in 1:j){
  #       E_JJJ = E_JJJ + smooth[jjj,iii]*log(P[jjj,jjj,i])
  #   }
  #   }
  # Q_func = E_J + E_JJ + E_JJJ
  # out[["Q_func"]] = Q_func

  
  return(out)
  
}

####### The M Step ######

# u hat is the smoothed probabilities
# v hat is the smoothed switching


M_Step <- function(EM_out, J,Y){
  
  Smooth=EM_out[["Smooth"]]
  SmoothSwitch=EM_out[["SmoothSwitch"]]
  
  delta = matrix(NA, nrow = J, ncol = J)
  
  for (j in 1:J){
    delta[j,] = Smooth[j,1]
  }
  
  gamma = matrix(NA, nrow = J, ncol = J)
  
  for (j in 1:J){
    for (k in 1:J){
      gamma[j,k] = sum(SmoothSwitch[j,k,]) / sum(SmoothSwitch[j,,])
    }
  }

  mu = matrix(NA, nrow = J, ncol = 1)
  sigma = matrix(NA, nrow = J, ncol = 1)
  
  for (j in 1:J){
    mu[j] = sum(Smooth[j,]*Y) / sum(Smooth[j,])
    sigma[j] = sqrt(sum(Smooth[j,]*(Y-mu[j])**2)/sum(Smooth[j,]))
    
  }
  
  mix = matrix(NA, nrow = J, ncol = J)
  
  mix[,1 ] = mu
  mix[,2 ] = sigma
  
  
  
  out = list()
  
  out[["delta"]] = delta
  out[["gamma"]] = gamma
  out[["mix"]] = mix
  
  
  return(out)
  
  
}


### Run ### 

run <- function(iterations, J, Y){

  # Initial mix
  mu_1 = mean(Y) * 0.7
  mu_2 = mean(Y) * 1.5
  std_1 = var(Y) * 0.7
  std_2 = var(Y) * 1.5
  
  mix = matrix(c(mu_1, mu_2, std_1, std_2), nrow = 2)
  
  # Initial p 
  p = matrix(c(0.95,0.05,0.05,0.95), byrow = TRUE, nrow = 2)
  
  delta = StaDist(p)
  
  for (i in 1:iterations){
    E <- EM(Y, J, p, delta, mix)
    print(E$LLK)
    print(p)
    M <- M_Step(E, J, Y)
    mix = M$mix
    p = M$gamma
    print(paste("Iteration:",i))
  }

  out = list()
  out[["Mix"]] = mix
  out[["p"]] = p
    
  return(out)
}


res = run(10, 2, y)

res$Mix

res$p

# The algorithm works!!!! 
# However, there is an issue. -Inf likelihood with large samples. 
# This is due to numerical underflow. We need a more robust way of computing low probabilities. 




