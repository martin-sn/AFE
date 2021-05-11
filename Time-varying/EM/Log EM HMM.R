###### Expectation Maximization for Hidden Markov Models #####

### Implementation with log likelihoods to avoid numerical underflow. !NOT DONE!

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


foo = rep(2,2)
dens = matrix(2, nrow = 1, ncol=2)

gam = matrix(1, nrow =2, ncol=2)

gam %*% t(dens*foo)


######## EM Algorithm for HMM ##########

# EM is designed to compute the maximum likelihood estimator in the case of missing values
# We consider the unobserved states of the Markov chain as missing data
# 1. Choose a starting value for theta
# 2. E-Step: Compute the conditional expectations of those functions 
# of the missing data appear in the complete-data log-likelihood
# 3. M-step: Maximization of the log-likelihood with respect to the est of parameters
# to be estimated (the missing data are substituted by their conditional expectation)
# 4. Iterate 2. and 3. until convergence. 


###### LOG Forward probabilities (zucchini et al 2019) ######

dnorm(1,c(1,2),c(2,3))

LogForwardProb <- function(Y, J, mix, Gamma, Delta, Densities){
  
  n = length(Y)
  alpha = matrix(0, nrow = J, ncol = n)
  foo = Delta*dnorm(Y[1], mean = mix[,1], mix[,2])
  sumfoo = sum(foo)
  lscale = log(sumfoo)
  foo <- foo/sumfoo
  alpha[,1] = lscale + log(foo)
  
  for (i in 2:n){
    
    dense = matrix(NA, nrow = 1, ncol = 2)
    dense[1] = Densities[1,1,i]
    dense[2] = Densities[2,2,i]
    foo = foo %*% Gamma * dense
    sumfoo = sum(foo)
    lscale = lscale + log(sumfoo)
    foo = foo/sumfoo
    alpha[,i] = lscale + log(foo)
  }
  
  return(alpha)
  
}

test = matrix(NA, nrow = 1, ncol = 2)


test[1] = 1
test[2] = 2

test
###### LOG Backwards probabilities (zucchini et al 2019) ######

LogBackwardProb <- function(Y, J, mix, Gamma, Delta, Densities){
  n = length(Y)
  beta = matrix(0, nrow = J, ncol = n)
  foo = rep(1/J,J)
  scale = log(J)
  
  for (i in (n-1):1){
    
    # Debugging
  #  print(i)
  #  print("densities")
  #  print(Densities[,,i+1])
  #  print("gamma")
  #  print(Gamma)
  #  print("foo")
  #  print(foo)
    dense = matrix(NA, nrow = 1, ncol = 2)
    dense[1] = Densities[1,1,i+1]
    dense[2] = Densities[2,2,i+1]
  #  print("dense")
  #  print(dense)
    foo = Gamma %*% (matrix((dense*t(foo)), nrow = 2))
    beta[,i] = log(foo) + scale
    sumfoo = sum(foo)
    foo = foo/sumfoo
    scale = scale + log(sumfoo)
    
  }
  
  return(beta)
  
  
}

rep(2,2) * rep(2,2)

rep(2,2)[1]

EM <- function(Y, J, p, delta, mix){
  
  n = length(Y)
  

  # Init forward probability array

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
  
  alpha <- LogForwardProb(Y,J,mix,p,delta,P)
  
  # Backwards probabilities
  
  beta = LogBackwardProb(Y,J,mix,p,delta,P)
  
  out = list()
  
  out[["alpha"]] = alpha
  out[["beta"]] = beta
  
  
  # Multiplying the forward and backward probabilities gives us the likelihood
  
  LK = matrix(NA, nrow = 2)
  
  #LK[1] = sum(alpha[1,]*beta[1,])
  #LK[2] = sum(alpha[2,]*beta[2,])
  
  LK[1] = sum(exp(alpha[1,]) * exp(beta[1,]))
  LK[2] = sum(exp(alpha[2,]) *exp(beta[2,]))
  
  
  out[["LLK"]] = log(LK)
  
  
  #### Smoothed probabilities ####
  
  smooth = array(NA, dim = c(J, length(Y)))
  
  for (i in 1:n){
    for (j in 1:J){
      smooth[j,i] = exp(alpha[j,i] + beta[j,i] - LogSumExp(alpha[,i] + beta[,i]))
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
        SmoothSwitch[k,j,i] = exp(alpha[j,i-1]) + log(p[j,k]) + log(P[k,k,i]) + beta[k,i] - LogSumExp(alpha[,i] + beta[,i])
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


matrix(c(1,2,3,4), ncol = 2) %*% t(matrix(c(0.02638976, 0.01271048), ncol=2) * rep(1/2,2))

res = run(10, 2, y)

res$Mix

res$p

# The algorithm works!!!! 
# However, there is an issue. -Inf likelihood with large samples. 
# Maybe we need a robust way of computing the log likelihood. 




