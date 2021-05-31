### Simulate Markov Switching GARCH Model ###

# Simulate GARCH #

SimGARCH <- function(n, Omega, Alpha, Beta){
  sd = vector(length = n)
  sd[1] = 1
  epsilon = rnorm(n)
  y = vector(length = n)
  
  y[1] = epsilon[1]
  for (i in 2:n){
    sd[i] = Omega + Alpha*y[i-1]**2 + Beta*sd[i-1]
    y[i] = epsilon[i]*sqrt(sd[i])
  }
  return(y)
}


sim <- SimGARCH(5000, 0.1, 0.3, 0.5)

plot(sim, type = "l")

# Source GARCH functions to estimate
source("Time-varying/MSGARCH/Utils/GARCH.r")

# Source Markov Chain Utils # 

source("Time-varying/HMM/Utils.r")

EstimateGARCH(sim)


### Simulate MSGARCH ###

# mGARCH matrix of GARCH Models #
SimMSGARCH <- function(n,mGARCH, P){
  J = ncol(mGARCH)
  
  sim = array(NA, c(n,J))

  for (j in 1:J){
    sim[,j] = SimGARCH(n,mGARCH[1,j], mGARCH[2,j], mGARCH[3,j])
  }
  HMM_SIM = SimulateHMM(P,n,J)
}


SimMSGARCH <- function(n,mGARCH, P){
  J = ncol(mGARCH)
  Y = vector(length = n)
  Y[1] = rnorm(1)
  
  HMM_SIM = SimulateHMM(P,n,J)
  
  SD <- array(NA, c(n,2))
  SD[1,] = 1
  
  
  
  for (i in 2:n){
    G = HMM_SIM[i]
    SD[i,1] = mGARCH[1,1] + mGARCH[2,1]*Y[i-1]**2 + mGARCH[3,1]*SD[i-1,1]
    SD[i,2] = mGARCH[1,2] + mGARCH[2,2]*Y[i-1]**2 + mGARCH[3,2]*SD[i-1,2]

    Y[i] = rnorm(1)*sqrt(SD[i,G])
  }
  return(Y)
}

GARCH_MIX = matrix(c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45), byrow=TRUE, nrow=3)
P = matrix(c(0.90, 0.05, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 0.90), byrow = TRUE, nrow = 3)

GARCH_MIX = matrix(c(0.05,0.1,0.2,0.25,0.35,0.4), byrow=TRUE, nrow=3)
P = matrix(c(0.95, 0.05, 0.15, 0.85), byrow = TRUE, nrow = 2)

P

GARCH_MIX



GARCH_MIX
sim = SimMSGARCH(1000, GARCH_MIX, P)

plot(sim, type = "l")


### Estimate MSGARCH by package ### 
library(MSGARCH)
spec <- CreateSpec()
fit <- FitMCMC(spec = spec, data = sim, ctr = list(nburn = 500L, nmcmc = 500L, nthin = 1L))
fit <- FitML(spec = spec, data = sim, ctr = list(nburn = 500L, nmcmc = 500L, nthin = 1L))
summary(fit)


### Estimate Manually ###

# Likelihood of a gausian MS Garch(1,1) #

MSGARCH_GAUSS_LLK <- function(Y,J,P,Omega1, Omega2, Alpha1, Alpha2, Beta1, Beta2){
  
  require(expm)
  
  n = length(Y)
  LLK = 0
  
  N = n
  
  sigma1 <- array(NA, n)
  sigma2 <- array(NA, n)
  
  # initialize # 
  
  StaDist = P %^% 1000
  
  
  sigma1[1] = 1.5
  sigma2[1] = 0.7
  
  #Omega1 = 0.05
  #Omega2 = 0.1
  
  #Alpha1 = 0.2
  #Alpha2 = 0.25
  
  #Beta1 = 0.35
  #Beta2 = 0.4
  
  
  for (i in 2:n){
    sigma1[i] <- Omega1 + Alpha1*Y[i-1]**2 + Beta1*sigma1[i-1]
    sigma2[i] <- Omega2 + Alpha2*Y[i-1]**2 + Beta2*sigma2[i-1]
  }
  
  Prob1 <- dnorm(Y, 0, sigma1) + 0.0001
  Prob2 <- dnorm(Y, 0, sigma2) + 0.0001
  
  Prob <- array(NA, dim = c(n,2))
  Prob[,1] = Prob1
  Prob[,2] = Prob2
  
  
  FilterProb1 <- Prob1 / (Prob1 + Prob2) + 0.0001
  FilterProb2 <- Prob2 / (Prob1 + Prob2) + 0.0001
  
  FilterProb <- array(NA, dim = c(n,2))
  
  FilterProb[1,1] = StaDist[1,1]*Prob[1,1] / (StaDist[1,1]*Prob[1,1] + StaDist[2,2]*Prob[1,2])
  FilterProb[1,2] = StaDist[2,2]*Prob[1,2] / (StaDist[1,1]*Prob[1,1] + StaDist[2,2]*Prob[1,2])

  for (i in 2:N){
    FilterProb[i,1] = FilterProb[i-1,1] * P[1,1] * Prob[i,1] + FilterProb[i-1,2] * P[2,1] * Prob[i,1] 
    FilterProb[i,2] = FilterProb[i-1,2] * P[2,2] * Prob[i,2] + FilterProb[i-1,1] * P[1,2] * Prob[i,1]
    FilterProb[i,1] = FilterProb[i,1] / sum(FilterProb[i,]) + 0.0001
    FilterProb[i,1] = FilterProb[i,2] / sum(FilterProb[i,]) + 0.0001
  }
  
  
  
  LLK = 0
  
  for (i in 2:N){
    for (j in 1:2){
      for (k in 1:2)
        LLK = LLK + log(P[j,k]*FilterProb[i,j]*Prob[i,k])
    }
  }
  return(LLK)
}

GARCH_MIX = matrix(c(0.05,0.1,0.2,0.25,0.35,0.4), byrow=TRUE, nrow=3)




Obj_func_MSGARCH <- function(vPar, Y){
  Omega1 = vPar[1]
  Omega2 = vPar[2]
  Alpha1 = vPar[3]
  Alpha2 = vPar[4]
  Beta1 = vPar[5]
  Beta2 = vPar[6]
  P11 = vPar[7]
  P12 = vPar[8]
  P21 = vPar[9]
  P22 = vPar[10]
  
  P = matrix(c(P11,P12, P21, P22), nrow = 2, byrow = TRUE)
  
  llk = MSGARCH_GAUSS_LLK(Y,2,P,Omega1, Omega2, Alpha1, Alpha2, Beta1, Beta2)
  return(-llk)
  
}





MSGARCH_est <- function(Y){
  #alpha = 0.8 
  #beta = 0.1
  #omega = var(y) * (1 - alpha - beta)
  
  require(Rsolnp)
  
  Omega1 = 0.05
  Omega2 = 0.01
  Alpha1 = 0.1
  Alpha2 = 0.8
  Beta1 = 0.2
  Beta2 = 0.1

  P11 = 0.95
  P12 = 0.05
  P21 = 0.20
  P22 = 0.80
  
  vPar = c(Omega1, Omega2, Alpha1, Alpha2, Beta1, Beta2, P11, P12, P21, P22)
  
  
  #New optimizer with fixed constraints 
  
  optimizer = solnp(vPar, fun = Obj_func_MSGARCH, ineqfun = function(vPar,...) {
    sum1 <- sum(vPar[3],vPar[5])
    sum2 <- sum(vPar[4], vPar[6])
    return(c(sum1,sum2))}
 , eqfun = function(vPar, ...){sum1 <- sum(vPar[7] + vPar[8])
 sum2 <- sum(vPar[9], vPar[10])
 return(c(sum1, sum2))}, eqB = c(1,1),
  
  
  ineqLB = c(1e-12,1e-12), ineqUB = c(0.99999999,0.99999999), 
  LB = c(1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12,1e-12), UB = c(0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999),Y=Y)
  
  
  
  par = optimizer$par
  llk = -optimizer$value
  
  output = list()
  output[["Par"]] = par
  output[["llk2"]] = llk
  
  
  return(output)
}


MSGARCH_est(sim)

GARCH_MIX


func_1 <- function(x,y){
  return(x+y)
}

func_2 <- function(x,y){
  return(x*y)
  
}
c(func_1, func_2)


### Estimate with MS GARCH Package ###

library(tidyverse)

SP500 <- read_csv("Time-varying/Data/sp500ret.csv", col_names = FALSE)

# Estimate
library(MSGARCH)
spec_gaus <- CreateSpec(variance.spec =list(model = c("sGARCH")), distribution.spec = list(distribution = c("norm")), switch.spec = list(K=2))
fit_gaus <- FitML(spec = spec_gaus, data = sim)
summary(fit_gaus)



Pred
# Plot conditional volatility 

plot(Volatility(fit_gaus))

plot(State(fit_gaus)$SmoothProb[,1,2,drop = TRUE])




spec_t <- CreateSpec(variance.spec =list(model = c("sGARCH")), distribution.spec = list(distribution = c("std")), switch.spec = list(K=2))
fit_t <- FitML(spec = spec_t, data = SP500$X2)
summary(fit_t)

