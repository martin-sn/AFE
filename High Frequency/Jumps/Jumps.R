library(tidyverse)

# 5 sec prices
SPY = read_csv("High Frequency/Data/SPY.csv", col_names = FALSE)
# we need to resample to 5 min

sample <- seq(from =1, to= 4681, by = 60)


SPY = SPY[sample,]

### Measuring Jumps ###

# Realized Volatility

RV = vector(length= ncol(SPY))

for (j in 1:ncol(SPY)){
  RV[j] = 0
 for (i in 2:nrow(SPY)){
   RV[j] = RV[j] + (SPY[[i,j]] - SPY[[i-1,j]])**2 
 } 
}



# Bipower Variation
BV = vector(length = ncol(SPY))

for (j in 1:ncol(SPY)){
  BV[j] = 0
  for (i in 3:nrow(SPY)){
    BV[j] = BV[j] + abs(SPY[[i-1,j]] - SPY[[i-2,j]])*abs(SPY[[i,j]] - SPY[[i-1,j]])
  } 
  BV[j] = pi/2*BV[j]
  #pi/2 is normalizing factor.
}



# Proportion of Jumps # 

PJ = vector(length = ncol(SPY))

for (j in 1:ncol(SPY)){
  PJ[j] = 1-BV[j]/RV[j]
}

PJ # We see negative values
# We can truncate to avoid the negative values. 

### Testing for Jumps ###

# BNS (2006) Jump test # (Slide 17 in Jumps Lecture)

# We need to select an destimator for integrated quarticity. 
# We can for example use tripower quarticity, quadpower quarticity
# or truncated realized volatility

# Quadpower quarticity

QPQ <- vector(length = ncol(SPY))

# Constants #
n = nrow(SPY)
m = sqrt(2/pi) # Expected value of abs(Norm(0,1))
theta = pi**2/4+pi-5

for (j in 1:ncol(SPY)){
  QPQ[j] = 0
  for (i in 5:nrow(SPY)){
    QPQ[j] = QPQ[j] + abs(SPY[[i,j]]-SPY[[i-1,j]])*abs(SPY[[i-1,j]]-SPY[[i-2,j]])*abs(SPY[[i-2,j]]-SPY[[i-3,j]])*abs(SPY[[i-3,j]]-SPY[[i-4,j]])
  } 
  QPQ[j] = theta*n/(m**4)*QPQ[j]
}

TStat <- array(NA, dim = ncol(SPY))


for (j in 1:ncol(SPY)){
  TStat[j] = sqrt(n)*(RV[j]-BV[j])/sqrt(QPQ[j]) 
}

TStat > qnorm(0.95)
TStat > qnorm(0.99)
# 5 Days at 5% significance, 3 days at 1%
# One sided test!!! So we cant just do dnorm.



