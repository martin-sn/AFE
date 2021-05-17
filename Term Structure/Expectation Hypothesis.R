library(tidyverse)

### Read Data ###

#The file moint.txt contains monthly data on zero-coupon interest rates (U.S. Interbank
# Rates, London) corresponding to maturities 1 through 12 months.

dat <- read_delim("Term Structure/Data/moint.txt", delim = " ")[-1,]

plot(dat$BBUSD1M)

plot(matrix(dat[1,-1]), type = "l")

plot(matrix(dat[157,-1]), type = "l")

dat[157,]

nrow(dat)

dat[1,-1]


array(dat[1,-1])


100/((1+0.08)^6)
