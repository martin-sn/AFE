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
  
  for (i in 1:(length(Y))){
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
  reg <- lm(sn_star~Y[2:(length(sn_star)+1)]) 
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


