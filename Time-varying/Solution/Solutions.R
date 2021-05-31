vY = read.csv("Time-varying/Data/sp500ret.csv", sep = ",", row.names = 1, header = FALSE, stringsAsFactors = FALSE)[[1]]

### Functions for EM of Gaussian mixture model
# this function computes log(sum(exp(vX))) in a robust way
LogSumExp <- function(vX) {
  dC = max(vX)
  dLogSumExp = dC + log(sum(exp(vX - dC)))
  return(dLogSumExp)
}

## vY is the data
## iJ is the number of mixture components
## itermax is the maximum number of iterations
## tol is the tollerance of the algorithm, when the log
## likelihood increment is lower than tol, then the algorithm stops.
EM_Gauss <- function(vY, iJ, itermax = 1e3, tol = 1e-6) {
  
  ## Initialization
  vMu     = mean(vY) * seq(0.7, 1.5, length.out = iJ)
  vSigma2 = var(vY) * seq(0.7, 1.5, length.out = iJ)
  vOmega  = rep(1/iJ, iJ)
  
  ## updated variables
  vMu_next = vMu
  vSigma2_next = vSigma2
  vOmega_next = vOmega
  
  # vector to store the likelihood values
  vLLK = numeric(itermax)
  
  iT = length(vY)
  # log(u_ij) variables
  mLogU = matrix(NA, iT, iJ)
  # counter of the algorithm
  iC = 1
  
  #likelihood increment
  eps = 1
  
  if(iJ == 1) {
    vMu = mean(vY)
    vSigma2 = var(vY)
    vOmega = 1
    dLLK = sum(dnorm(vY, vMu, sqrt(vSigma2), log = TRUE))
    vLLK = NULL
  } else {
    
    while(eps > tol & iC < itermax) {
      
      # E step
      ## compute the unnormalized log(u_jt) variables (the numerator of the expression for u_{jt})
      for (j in 1:iJ) {
        mLogU[, j] = dnorm(vY, vMu[j], sqrt(vSigma2[j]), log = TRUE) + log(vOmega[j])
      }
      ## compute the log likelihood
      dLLK = sum(apply(mLogU, 1, function(x) LogSumExp(x)))
      
      # normalize the u variables (divide by the sum)
      mLogU = t(apply(mLogU, 1, function(x) x - LogSumExp(x)))
      
      # the U variables
      mU = exp(mLogU)
      
      # the denominator in the M step
      vNormalizing = colSums(mU)
      # M step
      vMu_next = colSums(mU*vY)/vNormalizing
      for (j in 1:iJ) {
        vSigma2_next[j] = sum(mU[, j] * (vY-vMu_next[j])^2)
      }
      vSigma2_next = vSigma2_next/vNormalizing
      vOmega_next  =  vNormalizing/iT
      
      vLLK[iC] = dLLK
      if(iC > 1) {
        eps = vLLK[iC]-vLLK[iC-1] 
      }
      
      iC = iC + 1 
      
      ## updated variables
      vMu = vMu_next
      vSigma2 = vSigma2_next
      vOmega = vOmega_next
    }
    
    vLLK = vLLK[1:(iC-1)]
    dLLK = tail(vLLK, 1)
  }
  ## number of (free) parameters
  iN = 3*iJ  - 1
  BIC = iN * log(iT) - 2*dLLK
  
  return(list(vMu = vMu, vSigma2 = vSigma2, vOmega = vOmega, BIC = BIC, dLLK = dLLK, vLLK = vLLK))
  
}

##

# Estimate a Gaussian mixture model assuming that the returns are iid. Choose the
# number of components using BIC.

lFit = list()

for (iJ in 1:5) {
  lFit[[iJ]] = EM_Gauss(vY, iJ)
}

vBIC = sapply(lFit, function(x) x$BIC)

vBIC

vBIC = sapply(lFit, function(x) x$vMu)


vBIC


which.min(vBIC)
# we select J = 3
iJ = 3

hist(vY, nclass = 50, freq = FALSE)
lines(sort(vY), dnorm(sort(vY), mean(vY), sd(vY)), col = "red", lwd = 2)

vMu = lFit[[iJ]]$vMu
vSigma2 = lFit[[iJ]]$vSigma2
vOmega = lFit[[iJ]]$vOmega

vDensityMixture = sapply(sort(vY), function(y) sum(vOmega * dnorm(y, vMu, sqrt(vSigma2))))
lines(sort(vY), vDensityMixture, col = "blue", lwd = 2)
legend("topright", legend = c("Gaussian", "Mixture"), col = c("red", "blue"), lwd = 2)

### HMM

## compute the stationary distribution of the MC
getDelta <- function(mGamma) {
  
  iJ = ncol(mGamma)
  mI = diag(iJ)
  mU = matrix(1, iJ, iJ)
  vU = rep(1, iJ)
  
  foo = t(mI - mGamma + mU)
  
  vDelta = solve(foo) %*% vU;
  
  return(vDelta)
}

# Simulate from a HMM with J regimes conditional Gaussian distribution  
SimulateHMM <- function(iT, vMu, vSigma2, mGamma) {
  
  # compute the stationary distribution of the MC
  vDelta = getDelta(mGamma)
  
  # vector of states
  vS = numeric(iT)
  # vector of observations
  vY = numeric(iT)
  
  vS[1] = sample(1:iJ, 1, prob = vDelta)
  
  for (t in 2:iT) {
    vS[t] = sample(1:iJ, 1, prob = mGamma[vS[t-1], ]) 
    vY[t] = rnorm(1, mean = vMu[vS[t]], sd = sqrt(vSigma2[vS[t]]))
  }
  
  return(list(vS = vS, vY = vY))
  
}

###############################################################################
## Log forward probabilitites log(alpha_jt), from zucchini et al 2019
norm.HMM.lforward <- function(vY, iJ, vMu, vSigma2, mGamma, vDelta, mDensities){
  
  vSigma = sqrt(vSigma2)
  iT = length(vY)
  lalpha <- matrix(NA, iJ, iT) 
  foo    <- vDelta * dnorm(vY[1], mean = vMu, sd = vSigma)
  sumfoo <-sum(foo)
  lscale <- log(sumfoo)
  foo    <- foo/sumfoo
  lalpha[,1] <-lscale + log(foo)
  
  for (i in 2:iT){
    foo<- foo %*% mGamma * mDensities[i, ]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo<- foo/sumfoo
    lalpha[,i]<- lscale + log(foo)
  }
  return(lalpha)
  
}

###############################################################################
## ## Log backward probabilitites log(beta_jt), from zucchini et al 2019
norm.HMM.lbackward <- function(vY, iJ, vMu, vSigma2, mGamma, vDelta, mDensities){
  vSigma = sqrt(vSigma2)
  iT = length(vY)
  lbeta <- matrix(NA, iJ, iT)
  lbeta[,iT] <- rep(0, iJ)
  foo <- rep(1/iJ, iJ)
  lscale <- log(iJ)
  for (i in (iT-1):1){
    foo<- mGamma %*%( mDensities[i+1, ] * foo)
    lbeta[,i]<-log(foo) +  lscale
    sumfoo <- sum(foo)
    foo<- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
  
}

# EM algorithm
EM_HMM <- function(vY, iJ, itermax = 1e3, tol = 1e-8) {
  
  ## initialize
  mGamma = matrix(NA, iJ, iJ)
  diag(mGamma) = 0.95
  mGamma[!diag(iJ)] = 0.05/(iJ-1)
  
  vMu     = mean(vY) * seq(0.7, 1.5, length.out = iJ)
  vSigma2 = var(vY) * seq(0.7, 1.5, length.out = iJ)
  vDelta  = rep(1/iJ, iJ)
  
  ## update
  vMu_next     = vMu
  vSigma2_next = vSigma2
  vDelta_next  = vDelta
  mGamma_next  = mGamma
  
  # vector to store the likelihood values
  vLLK = numeric(itermax)
  
  iT = length(vY)
  # matrix with log(p(y_t | S_t = j))
  mLogDensities = matrix(NA, iT, iJ)
  
  # log(u_ij) variables
  mLogU = matrix(NA, iT, iJ)
  # counter of the algorithm
  iC = 1
  
  #likelihood increment
  eps = 1
  
  if(iJ == 1) {
    vMu = mean(vY)
    vSigma2 = var(vY)
    vDelta = 1
    mGamma = 1
    dLLK = sum(dnorm(vY, vMu, sqrt(vSigma2), log = TRUE))
    vLLK = NULL 
  } else {
    
    while(eps > tol & iC < itermax) {
      
      ## compute densities 
      for (j in 1:iJ) {
        mLogDensities[, j] = dnorm(vY, vMu[j], sqrt(vSigma2[j]), log = TRUE)
      }
      
      ################################### E step ###################################
      lalpha = norm.HMM.lforward(vY, iJ, vMu, vSigma2, mGamma, vDelta, exp(mLogDensities))
      lbeta  = norm.HMM.lbackward(vY, iJ, vMu, vSigma2, mGamma, vDelta, exp(mLogDensities))
      
      # Compute the \hat u_jt variables 
      
      mU = matrix(0, iT, iJ) 
      
      #compute the likelihood
      dLLK = LogSumExp(lalpha[, 1] + lbeta[, 1])
      
      for(j in 1:iJ) {
        mU[, j] = exp(lalpha[j, ] + lbeta[j, ] - dLLK)
      }
      
      # compute the v_ijt variables
      aV = array(0, dim = c(iJ, iJ, iT))
      
      for (t in 2:iT) {
        for (j in 1:iJ) {
          for (k in 1:iJ) {
            aV[j, k, t] = exp(lalpha[j, t-1] + lbeta[k, t] + log(mGamma[j, k]) + mLogDensities[t, k] - dLLK)
          }
        }
      }
      
      ################################### M step ###################################
      for (i in 1:iJ) {
        for (j in 1:iJ) {
          mGamma_next[i,j] = sum(aV[i,j, ])
        }
      }
      
      for (i in 1:iJ) {
        mGamma_next[i, ] = mGamma_next[i, ] / sum(mGamma_next[i, ])
      }
      
      vDelta_next = mU[1, ]
      
      # the denominator in the M step
      vNormalizing = colSums(mU)
      # M step
      vMu_next = colSums(mU*vY)/vNormalizing
      
      for (j in 1:iJ) {
        vSigma2_next[j] = sum(mU[, j] * (vY-vMu_next[j])^2)
      }
      
      vSigma2_next = vSigma2_next/vNormalizing
      
      vLLK[iC] = dLLK
      if(iC > 1) {
        eps = vLLK[iC]-vLLK[iC-1] 
      }
      
      iC = iC + 1 
      
      ## updated variables
      vMu     = vMu_next
      vSigma2 = vSigma2_next
      vDelta  = vDelta_next
      mGamma  = mGamma_next
    }
    
    vLLK = vLLK[1:(iC-1)]
    dLLK = tail(vLLK, 1)
    
  }
  
  ## number of (free) parameters
  #    mu and sigma       delta       gamma
  iN = 2*iJ          +   (iJ-1) +   iJ * (iJ - 1)
  BIC = iN * log(iT) - 2*dLLK
  
  return(list(vMu = vMu, vSigma2 = vSigma2, vDelta = vDelta, mGamma = mGamma, dLLK = dLLK, vLLK = vLLK, iter = iC, eps = eps, mU = mU, BIC = BIC))
  
}

# Data Generation 
iJ = 2

mGamma = matrix(c(0.97, 0.03, 0.01, 0.99), iJ, iJ, byrow = TRUE)

vMu = c(-4, 5)
vSigma2 = c(1.5, 0.7)^2

iT = 5000

## Simulate the data

Sim = SimulateHMM(iT, vMu, vSigma2, mGamma)

# Estimate a two regimes HMM model on the simulated data.
Fit = EM_HMM(Sim$vY, iJ)

## Compare estimated parameters with true ones
Fit$vMu
vMu

Fit$vSigma2
vSigma2

Fit$mGamma
mGamma

## Estimate the model and select J

lFit = list()
for (iJ in 1:5) {
  lFit[[iJ]] = EM_HMM(vY, iJ)
}

vBIC = sapply(lFit, function(x) x$BIC)

which.min(vBIC)
# we select J = 2
iJ = 2

Fit = lFit[[iJ]]

# Estimated parameters
Fit$vSigma2
Fit$vMu
Fit$mGamma

# plot smoothed probabilitites 
plot.ts(Fit$mU, main = "Smoothed probabilities")

## Price changes

## I noted that for some reasons price changes are not rounded, here I use round()
## absolute price changes
vY = abs(round(read.csv("Time-varying/Data/JPM.csv", sep = ",", row.names = 1, header = FALSE, stringsAsFactors = FALSE)[[1]]))

## EM for the ZIP model

# Density zip model
pZIP <- function(dY, dLambda, dPi, log = TRUE) {
  
  if (dY == 0) {
    dLPDF = log(dPi + (1-dPi) * exp(-dLambda))
  } else {
    dLPDF = log(1-dPi) - dLambda + dY * log(dLambda) - lfactorial(dY)
  }
  
  if (!log) {
    dLPDF = exp(dLPDF)
  }
  return(dLPDF)
  
}

EM_ZIP <- function(vY, itermax = 1e3, tol = 1e-8) {
  
  iT = length(vY)
  
  #initializations
  dLambda = mean(vY)
  dPi = sum(vY == 0)/iT
  
  ## update
  dLambda_next = dLambda
  dPi_next = dPi
  
  # vector to store the likelihood values
  vLLK = numeric(itermax)
  
  # counter of the algorithm
  iC = 1
  
  #likelihood increment
  eps = 1
  
  ## observations that are zero
  vZeros = which(vY == 0)
  
  while(eps > tol & iC < itermax) {
    
    ## log densities
    vLogDensities = sapply(vY, pZIP, dLambda = dLambda, dPi = dPi, log = TRUE)
    
    ## log likelihood
    dLLK = sum(vLogDensities)
    
    ## E step
    vB = numeric(iT)
    
    vB[vZeros] = exp(log(dPi) - vLogDensities[vZeros])
    
    ## MStep
    
    dPi_next = sum(vB)/iT
    dLambda_next = sum((1 - vB)*vY)/sum(1 - vB)
    
    vLLK[iC] = dLLK
    if(iC > 1) {
      eps = vLLK[iC]-vLLK[iC-1] 
    }
    
    iC = iC + 1 
    
    ## updated variables
    dLambda = dLambda_next
    dPi = dPi_next
    
  }
  
  vLLK = vLLK[1:(iC-1)]
  dLLK = tail(vLLK, 1)
  
  ## number of parameters
  iN = 2
  BIC = iN * log(iT) - 2*dLLK
  
  return(list(dPi = dPi, dLambda = dLambda, dLLK = dLLK, vLLK = vLLK, iter = iC, eps = eps, BIC = BIC))
  
}

Fit = EM_ZIP(vY)
Fit$dPi
Fit$dLambda

# MSGARCH

library(MSGARCH)

vY = read.csv("sp500ret.csv", sep = ",", row.names = 1, header = FALSE, stringsAsFactors = FALSE)[[1]]

# Specify the MSGARCH model
spec_MSGARCH = CreateSpec(variance.spec = list(model = "sGARCH"), 
                  distribution.spec = list(distribution = "norm"), 
                  switch.spec = list(K = 2))

# estimate by ML
Fit_MSGARCH = FitML(spec_MSGARCH, vY)

# extract the filtered volatility
vSigma_MSGARCH = Volatility(Fit_MSGARCH)

plot.ts(vSigma_MSGARCH, las = 1)

# Compare the estimated conditional volatility with that of a GARCH(1,1) model.

# Specify the GARCH model
spec_GARCH = CreateSpec(variance.spec = list(model = "sGARCH"), 
                          distribution.spec = list(distribution = "norm"), 
                          switch.spec = list(K = 1))

# estimate by ML
Fit_GARCH = FitML(spec_GARCH, vY)

# extract the filtered volatility
vSigma_GARCH = Volatility(Fit_GARCH)

lines(vSigma_GARCH, col = "red")
legend("topright", legend = c("MSGARCH", "GARCH"), col = c("black", "red"), lty = 1)

# Specify the MSGARCH model with t shocks
spec_MSGARCH_t = CreateSpec(variance.spec = list(model = "sGARCH"), 
                          distribution.spec = list(distribution = "std"), 
                          switch.spec = list(K = 2))

# estimate by ML
Fit_MSGARCH_t = FitML(spec_MSGARCH_t, vY)

which.min(c("t" = BIC(Fit_MSGARCH_t), "norm" = BIC(Fit_MSGARCH)))

# we select the norm specifications



