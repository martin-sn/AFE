library(tidyverse)

### Read Data ###

#The file moint.txt contains monthly data on zero-coupon interest rates (U.S. Interbank
# Rates, London) corresponding to maturities 1 through 12 months.

dat <- read_delim("Term Structure/Data/moint.txt", delim = " ", )[-1,]


moint <- read_table2("Term Structure/Data/moint.txt", col_types = cols(Code = col_skip(), BBUSD1M = col_double(), BBUSD2M = col_double(), BBUSD3M = col_double(), BBUSD4M = col_double(), BBUSD5M = col_double(), 
BBUSD6M = col_double(), BBUSD7M = col_double(), 
BBUSD8M = col_double(), BBUSD9M = col_double(), 
BBUSD10 = col_double(), BBUSD11 = col_double(), BBUSD12 = col_double()))


moint <- moint[-1,]       

plot(matrix(dat[1,-1]), type = "l")

plot(matrix(dat[157,-1]), type = "l")


### Data Exercise ### 

# We wish to test the expectations hypothesis using the Yield Spread test# #

# Work with equations 10.2.14 10.2.15 and 10.2.16 in CLM.
# Page 421

# Compute Yield Spread #
sn <- matrix(NA,nrow = nrow(dat), ncol=ncol(dat)-1)

for (i in 1:nrow(sn)){
  for (j in 1:ncol(sn)){
    sn[i,j] = log(moint[[i,j]]) - log(moint[[i,1]])
  }
}



# McCulloch and Kwon (1993) #
# https://www.asc.ohio-state.edu/mcculloch.2/ts/mcckwon/mccull.htm
# We need to get zero zield data from 1952 to 1991. The data is in a very bad structure, must be cleaned. 
library(readr)
zeroyld1 <- read_table2("Term Structure/Data/zeroyld1.txt", skip = 8, col_names = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"))
zeroyld2 <- read_table2("Term Structure/Data/zeroyld2.txt", skip = 16, col_names = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9"))


# Start at 1952
Zero_trimmed <- zeroyld1[-(1:488),]

# Get indicies of metadata
MetaIndicies <- sequence(nrow(Zero_trimmed)/8, from = 1, by = 8)
MetaIndicies_2 <- sequence(nrow(zeroyld2)/8, from = 1, by = 8)


# Put metadata in frame
Zero_Metadata <- Zero_trimmed[MetaIndicies, ]
Zero_Metadata_2 <- zeroyld2[MetaIndicies_2, ]

# Maturities in years
Maturities = c(0.000, 0.083, 0.167,0.250, 0.333, 0.417, 0.500, 
              0.583, 0.667,0.750, 0.833, 0.917, 1.000, 1.083,
              1.167, 1.250, 1.333,1.417, 1.500, 1.750, 2.000,
              2.500, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000,
              9.000, 10.000, 11.000, 12.000, 13.000, 14.000, 15.000,
              16.000, 17.000, 18.000, 19.000, 20.000, 21.000, 22.000,
              23.000, 24.000, 25.000, 26.000, 27.000, 28.000, 29.000,
              30.000, 31.000, 32.000, 33.000, 34.000, 35.000, 40.000)



# Remove metadata from yield data
Zero_nometa <- Zero_trimmed[-MetaIndicies,]
Zero_nometa_2 <- zeroyld2[-MetaIndicies_2,]

# Create 3d array for yield data
Zero_arr <- array(NA, dim = c(nrow(Zero_Metadata), 7,9))
Zero_arr_2 <- array(NA, dim = c(nrow(Zero_Metadata_2), 7,9))

# Get indicies for 3d arrays
Zero_seq <- sequence(nrow(Zero_nometa)/7, from =1, by =7)
Zero_seq_2 <- sequence(nrow(Zero_nometa_2)/7, from =1, by =7)


# Convert 2d array to 3d array
for (i in 1:dim(Zero_arr)[1]){
  for (j in 1:dim(Zero_arr)[2]){
    for (k in 1:dim(Zero_arr)[3]){
      Zero_arr[i,j,k] = Zero_nometa[[Zero_seq[i]+j-1,k]]
    }
  }
}

for (i in 1:dim(Zero_arr_2)[1]){
  for (j in 1:dim(Zero_arr_2)[2]){
    for (k in 1:dim(Zero_arr_2)[3]){
      Zero_arr_2[i,j,k] = Zero_nometa_2[[Zero_seq_2[i]+j-1,k]]
    }
  }
}


# Convert 3d array back to 2d array
New_arr = array(Zero_arr, dim = c(dim(Zero_arr)[1], dim(Zero_arr)[2]*dim(Zero_arr)[3]))
New_arr_2 = array(Zero_arr_2, dim = c(dim(Zero_arr_2)[1], dim(Zero_arr_2)[2]*dim(Zero_arr_2)[3]))


for (i in 1:nrow(New_arr)){
  New_arr[i,] = array(t(Zero_arr[i,,]))
}

for (i in 1:nrow(New_arr_2)){
  New_arr_2[i,] = array(t(Zero_arr_2[i,,]))
}

# FINAL ARRAY THAT WORKS! # 

Final_arr <- rbind(New_arr, New_arr_2)


# Compute Yield Spread #
sn <- matrix(NA,nrow = nrow(New_arr), ncol=1)

for (i in 1:nrow(sn)){
  for (j in 1:ncol(sn)){
    sn[i,j] = log(New_arr[[i,3]]) - log(New_arr[[i,2]])
  }
}


### Beta test ### 

# Equation (10.2.16)
LongYield <- function(y1, yn, yn_minus, n){
  Y = array(NA, dim = (length(y1)-1))
  
  for (i in 1:(length(Y))){
    Y[i] = yn_minus[i+1] - yn[i]
  }
  
  sn = (yn-y1)/(n-1) # Compute spread and divide by n-1
  reg = lm(reg <- lm(Y ~ sn[1:(length(sn)-1)])) #  Equation 10.2.16
  
  return(reg)
}

LongYield(Final_arr[,2], Final_arr[,13], Final_arr[,12], 12)


### Gamma Test ###  From slides

ShortRate <- function(y1,yn,yn_minus, n){
  sn = array(NA, dim = (length(y1))) 
  
  for (i in 1:(length(y1))){
    sn[i] = yn[i] - y1[i] # Yield spread
  }
  sn_star <- array(NA, dim = (length(y1)-n))
  
  for (i in 2:(length(sn_star)+1)){
    sn_star_i = array(NA, dim = (n))
    for (j in 0:(n-1)){
      sn_star_i[j+1] = (1-j/n)*(y1[i+j] - y1[i+j-1]) 
    }
    
    sn_star[i-1] = sum(sn_star_i)
    
  }
  reg <- lm(sn_star~sn[2:(length(sn_star)+1)]) 
  return(reg)
}

ShortRate(Final_arr[,2], Final_arr[,3], Final_arr[,2], 2)



### Gamma Test ### From book

ShortRate <- function(y1,yn,yn_minus, n){
  sn = array(NA, dim = (length(y1)))
  
  for (i in 1:(length(sn))){
    sn[i] = yn[i] - y1[i] 
  }
  sn_star <- array(NA, dim = (length(y1)-n))
  
  #sn_star = equation 10.2.17
  for (i in 1:(length(sn_star))){
    sn_star_i = array(NA, dim = (n-1))
    for (j in 1:(n-1)){
      sn_star_i[j] = (1-j/n)*(y1[i+j] - y1[i+j-1]) 
    }
    sn_star[i] = sum(sn_star_i)
  }
  
  reg <- lm(sn_star~sn[1:length(sn_star)]) # Equation 10.2.18 
  
  return(reg)
}


# The equations from the book work

reg <- ShortRate(Final_arr[,2], Final_arr[,3], Final_arr[,2], 2)

summary(reg)

plot(reg$residuals, type = "l")


##### Vasicek vs CIR ##### 


# Fix the moint format because it is weird #

Moint_1 <- as.vector(moint[,1])
moint_arr <- array(NA, dim = nrow(Moint_1))

for (i in 1:nrow(Moint_1)){
  moint_arr[i] = Moint_1[[i,1]]
}


moint_matrix <- array(NA, dim = dim(moint))

for (i in 1:nrow(moint)){
  for (j in 1:ncol(moint)){
    moint_matrix[i,j] = moint[[i,j]]
    
    
  }
  
}


# We regress Y1,t+1 - y1,t on y_1t, square the residuals

HM <- function(y1){
  sn <- diff((y1))
  y1 = y1[1:(length(y1)-1)]
  
  reg <- lm(sn ~ y1)
  
  reg_resid <- lm(reg$residuals**2 ~ (y1**2))
  
  
  return(reg_resid)
  
}

TestMoint <- HM(moint_arr)

summary(TestMoint)

# Coefficient for y1 is insigificant (barely) The Vasicek model is favored


TestLargeData <- HM(Final_arr[,2])

summary(TestLargeData)

# coefficient is way more significant in this one! 


### Model forward rates? using the Vasicek (Homoscedastic model) ###

HomoModel <- function(y_matrix, maturities){
  
  # y1,t+1 - y1,t and run equation 8 (lecture 2 slides) regression
  sn <- diff(y_matrix[,1])
  y1 = y_matrix[1:(nrow(y_matrix)-1),1]
  reg <- lm(sn ~ y1)
  
  # Compute observed forward rates (slide 7 lecture 1)
  
  ForwardRate <- function(yn1, yn, n){
    fw <- yn + (n+1)*(yn1-yn)
    return(fw)
  }
  
  FW_Matrix <- array(NA, dim = dim(y_matrix))
  
  for (i in 1:(ncol(y_matrix)-1)){
    FW_Matrix[,i] = ForwardRate(y_matrix[,i+1], y_matrix[,i],maturities[i]*12)
  }
  
  n = ncol(y_matrix)
  b = reg$coefficients[[2]]
  phi = b+1
  sigma2 <- var(reg$residuals)
  avg_diff <- mean(y_matrix[,n] - y_matrix[,1])
  
  beta_func <- function(beta){
    res = (1/(1-phi))**2*sigma2/2+(beta/(1-phi))*sigma2 + avg_diff
    return(res)
  }
  beta <- uniroot(beta_func, lower = -1000, upper = 1000)$root
  
  
  mu_func <- function(mu){
    res = (1-phi)*(mu-0.5*beta**2*sigma2) - reg$coefficients[[1]]
  }
  
  mu <- uniroot(mu_func, lower = -100, upper = 100)$root
  
  fnt <- function(n, mu, beta, phi, sigma2, xt){
    fnt <- array(NA, dim = length(xt))
    fnt <- mu - (beta + (1-phi**n)/(1-phi))**2*sigma2/2+(phi**n)*(x_t-mu)
    return(fnt)
  }
    
  
  x_t <- y_matrix[,1] + 0.5*beta**2*sigma2
  
  fnt_arr <- array(NA, dim = c(length(x_t), length(maturities)))
  

  for (i in 1:length(maturities)){
    fnt_arr[,i] <- fnt(maturities[i], mu,beta,phi,sigma2,x_t)
  }

  out = list()
  
  out[["fnt_arr"]] = fnt_arr
  out[["FW_Matrix"]] = FW_Matrix
  
  return(out)
}


Hom <- HomoModel(Final_arr[,2:13], Maturities[2:13])

plot(y=Hom$fnt_arr[150,1:12],x = Maturities[1:12], type = "l")
lines(Hom$FW_Matrix[150,1:12], x = Maturities[1:12], type = "l", col = "red")

plot(Final_arr[1,2:30])
lines(Hom$FW_Matrix[1,1:12], col = "red")

MeanForward = array(NA, dim = 12)
MeanForwardPred = array(NA, dim = 12)


for (i in 1:12){
  MeanForwardPred[i] = mean(Hom$fnt_arr[,i])
}

for (i in 1:12){
  MeanForward[i] = mean(Hom$FW_Matrix[,i])
}


plot(MeanForward, type = "l", ylim = c(5,15))
lines(MeanForwardPred, type = "l", col = "red")

Moint_Maturities <- c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12, 7/12, 8/12, 9/12, 10/12, 11/12, 12/12)

Hom <- HomoModel(moint_matrix, Moint_Maturities)


plot(y=Hom$fnt_arr[150,1:12],x = Maturities[1:12], type = "l", xlim=c(0,1), ylim =c(4,10))
lines(Hom$FW_Matrix[150,1:12], x = Maturities[1:12], type = "l", col = "red")

MeanForward = array(NA, dim = 12)
MeanForwardPred = array(NA, dim = 12)


for (i in 1:12){
  MeanForwardPred[i] = mean(Hom$fnt_arr[,i])
}

for (i in 1:12){
  MeanForward[i] = mean(Hom$FW_Matrix[,i])
}


plot(MeanForward, type = "l", ylim = c(5,10))
lines(MeanForwardPred, type = "l", col = "red")

MeanForwardPred
MeanForward



plot(moint_matrix[150,1:12])
lines(Hom$FW_Matrix[150,1:12], col = "red")




#### Nelson-Siegel ####

NelsonSiegel <- function(y_matrix, maturities){
  
  
  # Here t actually needs to be estimated. But for simplicity we can leave this out (Diebold & Li 2006)
  # https://www.sas.upenn.edu/~fdiebold/papers/paper49/Diebold-Li.pdf
  t = 0.3
  maturities = maturities*12
  b_1 <- array(NA, dim = ncol(y_matrix))
  b_2 <- array(NA, dim = ncol(y_matrix))
  
  for (i in 1:ncol(y_matrix)){
    b_1[i] = (1-exp(-maturities[i]/t))/(maturities[i]/t)
    b_2[i] = (1-exp(-maturities[i]/t))/(maturities[i]/t) - exp(-maturities[i]/t)
  }
  
  fw_matrix <- array(NA, dim = dim(y_matrix))
  
  for (i in 1:nrow(fw_matrix)){
    reg <- lm(y_matrix[i,]~b_1 + b_2)
    
    fw_matrix[i,] <- reg$fitted.values
  }

  
  
  
  
  return(fw_matrix)
  
}



NelsonSiegelPackage <- function(y_matrix, maturities){
  require(YieldCurve)
  
  fit <- Nelson.Siegel(y_matrix,maturities)
  
  Fit_arr <- array(NA, dim = dim(y_matrix))
  
  beta_0 <- fit[,1]
  beta_1 <- fit[,2]
  beta_2 <- fit[,3]
  t <- fit[,4]
  
  for (j in 1:ncol(y_matrix)){
    m = maturities[j]
    
    Fit_arr[,j] = beta_0 + beta_1*(1-exp(-m/t))/(m/t) + beta_2*((1-exp(-m/t))/(m/t) - exp(-m/t))
    
    
  }
  
  return(Fit_arr)
}




#### Nelson Siegel With Moint ####

reg <- NelsonSiegel(moint_matrix, c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12, 7/12, 8/12, 9/12, 10/12, 11/12, 12/12))

test <- c(1/12, 2/12, 3/12, 4/12, 5/12, 6/12, 7/12, 8/12, 9/12, 10/12, 11/12, 12/12)

test*12
reg

MeanForward = array(NA, dim = 12)
MeanForwardPred = array(NA, dim = 12)
MeanForwardPredNelson = array(NA, dim = 12)

for (i in 1:12){
  MeanForwardPredNelson[i] = mean(reg[,i])
}

for (i in 1:12){
  MeanForward[i] = mean(reg[,i])
}

MeanForward
MeanForwardPred
plot(MeanForward, type = "l", ylim = c(5,10))
lines(MeanForwardPred, type = "l", col = "red")
lines(MeanForwardPredNelson, type = "l", col = "blue")


plot(reg$fitted.values, ylim = c(8,10), type = "l")

lines(Hom$FW_Matrix[1,], col = "red")

lines(Hom$fnt_arr[1,], col = "blue")

# With large data # 
 
reg <- NelsonSiegel(Final_arr[,2:30], Maturities[2:30])


for (i in 1:12){
  MeanForwardPredNelson[i] = mean(reg[,i])
}

for (i in 1:12){
  MeanForwardPred[i] = mean(Hom$fnt_arr[,i])
}

for (i in 1:12){
  MeanForward[i] = mean(Hom$FW_Matrix[,i])
}

Hom <- HomoModel(Final_arr[,2:30], Maturities[2:30])

Maturities[25]*12

plot(Hom$FW_Matrix[123,], type = "l")
lines(reg[123,])



plot(Final_arr[123,2:30])

lines(reg[123,])



reg

reg


plot(MeanForward, type = "l", ylim = c(5,10))
lines(MeanForwardPred, type = "l", col = "red")
lines(MeanForwardPredNelson, type = "l", col = "blue")




Fit <- NelsonSiegelPackage(moint_matrix, Moint_Maturities*12)

plot(Fit[1,], type = "l", ylim = c(0,10))
lines(reg[1,], type = "l")
plot(Hom$FW_Matrix[1,], type = "l")
plot(Hom$fnt_arr[1,])




lines(moint_matrix[1,])

moint_matrix[1,]



library(YieldCurve)

nel_curve <- Nelson.Siegel(Final_arr[,2:30], Maturities[2:30]*12)

coef <- c(manel_curve[1,][[1]], nel_curve[1,][[2]], nel_curve[1,][[3]], nel_curve[1,][[4]])

?NSrates

arr <- array(NA, dim = dim(nel_curve))
NSrates(nel_curve, 10)

for (i in 1:nrow(nel_curve)){
  for (j in 1:ncol(nel_curve)){
    arr[i,j] <- nel_curve[[i,j]]
  }
  
}
nel_curve[[1,2]]
time(nel_curve)

NSrates(arr,Maturities[2:30]*12)
