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
  SD = vector(length = n)
  
  SD[1] = 1
  Y[1] = rnorm(1)
  
  HMM_SIM = SimulateHMM(P,n,J)
  
  for (i in 2:n){
    G = HMM_SIM[i]
    SD[i] = mGARCH[1,G] + mGARCH[2,G]*Y[i-1]**2 + mGARCH[3,G]*SD[i-1]
    Y[i] = rnorm(1)*sqrt(SD[i])
  }
  return(Y)
}

GARCH_MIX = matrix(c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45), byrow=TRUE, nrow=3)
P = matrix(c(0.90, 0.05, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 0.90), byrow = TRUE, nrow = 3)

GARCH_MIX = matrix(c(0.05,0.1,0.2,0.25,0.35,0.4), byrow=TRUE, nrow=3)
P = matrix(c(0.90, 0.1, 0.1, 0.90), byrow = TRUE, nrow = 2)


GARCH_MIX
sim = SimMSGARCH(10000, GARCH_MIX, P)

plot(sim, type = "l")


### Estimate MSGARCH by package ### 
library(MSGARCH)
spec <- CreateSpec()
fit <- FitMCMC(spec = spec, data = sim, ctr = list(nburn = 500L, nmcmc = 500L, nthin = 1L))
fit <- FitML(spec = spec, data = sim, ctr = list(nburn = 500L, nmcmc = 500L, nthin = 1L))
summary(fit)


### Estimate Manually ###

# Likelihood of a gausian MS Garch(1,1) #

MSGARCH_GAUSS_LLK <- function(Y,J,P){
  
  n = length(Y)
  LLK = 0
  
  for (i in 1:N){
    for (j in 1:J){
      for (k in 1:J){
        LLK = LLK + log(P[j,k]) + log()
      }
    }
  }
}





### Estimate with MS GARCH Package ###

library(tidyverse)

SP500 <- read_csv("Time-varying/Data/sp500ret.csv", col_names = FALSE)

# Estimate
library(MSGARCH)
spec_gaus <- CreateSpec(variance.spec =list(model = c("sGARCH")), distribution.spec = list(distribution = c("norm")), switch.spec = list(K=2))
fit_gaus <- FitML(spec = spec_gaus, data = SP500$X2)
summary(fit_gaus)



Pred
# Plot conditional volatility 

plot(Volatility(fit_gaus))

plot(State(fit_gaus)$SmoothProb[,1,2,drop = TRUE])




spec_t <- CreateSpec(variance.spec =list(model = c("sGARCH")), distribution.spec = list(distribution = c("std")), switch.spec = list(K=2))
fit_t <- FitML(spec = spec_t, data = SP500$X2)
summary(fit_t)

