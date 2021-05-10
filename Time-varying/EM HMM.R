###### Expectation Maximization for Hidden Markov Models #####

### Simulate T = 5000 observations from a HMM with J = 2 regimes and trasition probability matrix

dat = c(0.9, 0.1, 0.1, 0.9)
p = matrix(dat,nrow = 2, byrow = TRUE)


sim_hmm = function(prob_matrix, n, j){
  state = array(NA, n)
  
  state[1] = 1 
  
  for (i in 2:n){
    state[i] = sample(c(1:j),1, prob=prob_matrix[state[i-1],])
    
  }
  
  return(state)
  
}



#testing new sim
dat = c(0.65, 0.35, 0.29, 0.71)
p = matrix(dat,nrow = 2, byrow = TRUE)

sim = sim_hmm(p, 1000, 2)

plot(sim, type = "l")


mu_1 = 1.2
mu_2 = 1.0
sigma_1 = 0.9
sigma_2 = 1.2
s_1 = 0.9
s_2 = 1.2

sim[1] == 1

array(NA, length(sim))

rnorm(1,-4,1.5)

sim

?rnorm

pnorm(1,-4,1.5)

dnorm(1,-4,1.5)




### Generate Y ###
GenY <- function(HMM_Sim, mu_1, mu_2, s_1,s_2){
  
  y = array(NA, length(HMM_Sim))
  
  for (i in 1:length(HMM_Sim)){
    if (HMM_Sim[i] == 1){
      y[i]= rnorm(1, mu_1, s_1)
    }
    else {
      y[i] = rnorm(1, mu_2, s_2)
      
    }
    
  }
  
  return(y)
}



y <- GenY(sim, mu_1, mu_2, sigma_1, sigma_2)

plot(y, type = "l")

###### Numerical Estimation (Not the best way) ######

# Matrix LLK from slide 50 HMM
LLK_Matrix = function(p, y, mu_1, mu_2, s_1, s_2){
  Sta_dist <- p%^% 10000
  
  P_Y = matrix(c(dnorm(y[1], mu_1, s_1, log = TRUE), 0, 0, dnorm(y[1], mu_2, s_2, log = TRUE)), byrow = TRUE, nrow = 2)
  LLK = log(Sta_dist) + P_Y
  
  log_gamma = log(p)
  
  for (i in 2:length(y)){
    
    P_Y = matrix(c(dnorm(y[i], mu_1, s_1, log = TRUE), 0, 0, dnorm(y[i], mu_2, s_2, log = TRUE)), byrow = TRUE, nrow = 2)
    
    LLK = LLK + log_gamma + P_Y
    
  }
  
  return(LLK %*% matrix(1, nrow = 2)) # Not sure about this 1 
}


# Different approach
LLK_Matrix = function(p, y, mu_1, mu_2, s_1, s_2){
  Sta_dist <- p%^% 10000
  
  P_Y = matrix(c(dnorm(y[1], mu_1, s_1, log = FALSE), 0, 0, dnorm(y[1], mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)
  
  LLK = Sta_dist %*% P_Y
    
  #LOG_LLK = log(LLK)
  
  #log_gamma = log(p)
  
  for (i in 2:length(y)){
    
    P_Y = matrix(c(dnorm(y[i], mu_1, s_1, log = FALSE), 0, 0, dnorm(y[i], mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)
    
    LLK = LLK %*% p %*% P_Y
    
    #LOG_LLK = LOG_LLK + log(LLK)
    
  }
  
  return(LLK %*% array(1, 2)) # Not sure about this 1 
}



LLK_Matrix(p, y, mu_1, mu_2, s_1, s_2)

LLK_Matrix(p, y, -3, 2, 1, 1) # <- Using other means than what the data was generated with gives better results...





## Testing below.... ##

P_Y = matrix(c(dnorm(1, mu_1, s_1, log = FALSE), 0, 0, dnorm(1, mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)

max(P_Y)

P_Y

Sta_dist * P_Y



P_Y + log(p%^% 10000)


P_Y = matrix(c(dnorm(1, mu_1, s_1, log = FALSE), 0, 0, dnorm(1, mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)


log(p%^% 10000 * P_Y)




######## EM with forward and backwards ##### 

EM <- function(Y, J, p,delta, mix){
  
  n = length(Y)
  
  # Init forward probability array
  alpha <- array(NA, c(J, n))
  
  
  #### Gaussian Mixture ###
#  mu_1 = -4
#  mu_2 = 5
#  std_1 = 1.5
#  std_2 = 0.7
  
 # mix = matrix(c(mu_1, mu_2, std_1, std_2), nrow = 2)
  
  mu_1 = mix[1,1]
  mu_2 = mix[2,1]
  std_1 = mix[1,2]
  std_2 = mix[2,2]
  
  # We find the limiting (stationairy) distrubtion of the transition matrix
  require(expm)
  sta_p = delta

  # First forward probability
  for (j in 1:J){
    alpha[j,1] = sta_p[1,j]*dnorm(y[1], mix[j,1], mix[j,2])
  }
  
  
  # Remaining forward probabilities
  
  # P(y_1) is a diagonal matrix with j,j-th element p_j(y_1)
  # p_j(y_1) is the density for gaussian j
  
  for (i in 2:length(Y)){
    P_1 = dnorm(y[i], mix[1,1],mix[1,2])
    P_2 = dnorm(y[i], mix[2,1], mix[2,2])
    
    P = matrix(c(P_1,0,0,P_2), nrow = 2)
    
    alpha[,i] = alpha[,i-1]%*%p%*%P
  }
  
  # Backwards probabilities
  
  beta = array(NA, dim = c(J, length(Y)))
  
  beta[,length(Y)] = matrix(c(1/n,1/n)) # Initialization. Not sure about this one. 
  
  for (i in (length(Y)-1):1){
    
    P_1 = dnorm(y[i], mix[1,1],mix[1,2])
    P_2 = dnorm(y[i], mix[2,1], mix[2,2])
    P = matrix(c(P_1,0,0,P_2), nrow = 2)
    
    beta[,i] = p%*%P%*%beta[,i+1]
  }
  
  out = list()
  
  out[["alpha"]] = alpha
  out[["beta"]] = beta
  
  
  # Multiplying the forward and backward probabilities gives us the likelihood
  
  LK = matrix(NA, nrow = 2)
  
  LK[1] = sum(alpha[1,]*beta[1,], axis = 0)
  LK[2] = sum(alpha[2,]*beta[2,], axis = 0)
  
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
  
  P_matrix = array(NA, c(J,J,n))

  
  for (i in 2:n){
    P_1 = dnorm(y[i], mix[1,1],mix[1,2])
    P_2 = dnorm(y[i], mix[2,1], mix[2,2])
    P_matrix[,,i] = matrix(c(P_1,0,0,P_2), nrow = 2)
  
  
  SmoothSwitch = array(NA, dim = c(J,J,n)) # Not sure about this 
  
  SmoothSwitch[,,1] = p
  
  for (i in 2:n){
    P_1 = dnorm(y[i], mix[1,1],mix[1,2])
    P_2 = dnorm(y[i], mix[2,1], mix[2,2])
    P = matrix(c(P_1,P_2), nrow = 2)
    for (j in 1:J){
      for (k in 1:J){
        
        
        SmoothSwitch[k,j,i] = alpha[j,i-1] * p[j,k] * P[k] * beta[k,i] / (t(alpha[,i])%*%beta[,i])
        
      
        
      }
    }
  }
  
  
  out[["SmoothSwitch"]] = SmoothSwitch
    
  
  #### Q Func ####
  
  # It seems that we actually do not need to compute the Q func. 
  
  # P(u) is equivalent to smooth prob
  # P(v) is equivalent to smooth switch
  
  
  Q_func = 0 

  E_J = 0
  E_JJ = 0
  E_JJJ = 0
  
  for (jj in 1:J){
    E_J = E_J + sum(smooth[1,jj]*log(sta_p[jj]))
    for (k in 1:J){
      E_JJ = E_JJ + sum(SmoothSwitch[jj,k,]) * log(p[jj,k])
    }
  }
    
    for (iii in 2:n){
      for (jjj in 1:j){
        E_JJJ = E_JJJ + smooth[jjj,iii]*log(P_matrix[jjj,jjj,i])
        
    }
    
    
    }
  
  
    
    
  }
  
  Q_func = E_J + E_JJ + E_JJJ

  
  out[["Q_func"]] = Q_func

  
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

M = M_Step(res, 2, y)

M$mix

### Run ### 


seq(0.7,1.5, length.out = 2)

run <- function(iterations, J, Y){
  
  # Inital mix
  mu_1 = mean(Y) * 0.7
  mu_2 = mean(Y) * 1.5
  std_1 = var(Y) * 0.7
  std_2 = var(Y) * 1.5
  
  mix = matrix(c(mu_1, mu_2, std_1, std_2), nrow = 2)
  
  # Inital p 
  
  p = matrix(c(0.95,0.05,0.05,0.95), byrow = TRUE, nrow = 2)
  
  delta = p %^% 10000
  
  for (i in 1:iterations){
    E <- EM(Y, J, p, delta, mix)
    print(E$LLK)
    print(p)
    M <- M_Step(E, J, Y)
    mix = M$mix
    p = M$gamma
  }

  out = list()
  out[["Mix"]] = mix
  out[["p"]] = p
    
  return(out)
}


res = run(50, 2, y)

res$Mix

res$p

# The algorithm works!!!! 
# However, there is an issue. -Inf likelihood with large samples. 
# Maybe we need a robust way of computing the log likelihood. 


######## EM Algorithm for HMM ######## 

# EM is designed to compute the maximum likelihood estimator in the case of missing values
# We consider the unobserved states of the Markov chain as missing data
# 1. Choose a starting value for theta
# 2. E-Step: Compute the conditional expectations of those functions 
# of the missing data appear in the complete-data log-likelihood
# 3. M-step: Maximization of the log-likelihood with respect to the est of parameters
# to be estimated (the missing data are substituted by their conditional expectation)
# 4. Iterate 2. and 3. until convergence. 

