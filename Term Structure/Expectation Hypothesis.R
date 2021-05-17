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
dat[-1,]$


datplot(dat$BBUSD1M)

plot(matrix(dat[1,-1]), type = "l")

plot(matrix(dat[157,-1]), type = "l")

dat[157,]

nrow(dat)

dat[1,-1]


array(dat[,-1])


100/((1+0.08)^6)

dat[,-1]




### Data Exercise ### 

# Compute Yield Spread # 

yield <- matrix(data = dat[,-1], )

as.matrix(as.double(dat[,-1]))


dat

yield

# Work with equations 10.2.14 10.2.15 and 10.2.16 in CLM.
# Page 421


sn <- matrix(NA,nrow = nrow(dat), ncol=ncol(dat)-1)

for (i in 1:nrow(sn)){
  for (j in 1:ncol(sn)){
    sn[i,j] = moint[[i,j]] - moint[[i,1]] 
  }
}




