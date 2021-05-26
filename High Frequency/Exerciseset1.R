library(tidyverse)

# 5 Sec prices
SPY <- read_csv("High Frequency/Data/SPY.csv", col_names = FALSE)

SPY <- read_csv("High Frequency/Data/SPY2013.csv", col_names = FALSE)


# We need to resample to 5 min
sample <- seq(from =1, to= 4681, by = 60)
SPY_5min = SPY[sample,]

SPY_5min
# Realized Quarticity # 


RealizedQuarticity <- function(X){
  X = log(X)
  RQ = 0
  n = length(X)
  
  for (i in 2:n){
    RQ = RQ + abs(X[i]-X[i-1])**4
  }
  RQ = RQ*(n-1)
  
  return(RQ)
}


RQ_1 <- RealizedQuarticity(SPY_5min$X1)

PreAvg <- function(Y){
  Y = log(Y)
  n = length(Y)
  omega_sq = 0
  for (i in 2:length(Y)){
    omega_sq = omega_sq + (Y[i]-Y[i-1])**2
  }
  omega_sq = 1/(2*(n-1))*omega_sq
  
  
  return(omega_sq)
  
}


Omega <- PreAvg(SPY$X1)


RQ_1/3/(4*Omega**2)


(RQ_1/(4*Omega**2))^(1/3)


Omega

RQ = 1.5846e-11
Omega = 1.6575e-10

(RQ/3)/(4*(Omega)**2)


(RQ/(4*Omega**2))^(1/3)

Omega_4 <- Omega**2

(RQ/3)/(4*Omega_4)

(RQ/(4*Omega_4))**(1/3)


(RQ/3)/(4*(Omega**2))



length(SPY_5min[,1, drop = TRUE])

### Find n for each of the 21 days ###

n = array(NA,dim = 252)

for(i in 1:252){
  RQ <- RealizedQuarticity(SPY_5min[,i, drop = TRUE])
  Omega = PreAvg(SPY[,i, drop = TRUE])
  n[i] = (RQ/(4*Omega**2))**(1/3)

}

n

mean(n)

(RQ_1/(4*Omega*Omega))**(1/3)


SPY_5min[,"X1"]$X1

col = paste("X",1, sep = "")

SPY_5min$(paste("X",1, sep = ""))


SPY_5min[,1, drop = TRUE]



