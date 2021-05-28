### Simple Term Structure Model (From video) ###

mean_moint = mean(moint$BBUSD1M)

y1 <- moint[2:(nrow(moint)),1]$BBUSD1M - mean_moint

lagy1 <- moint[1:(nrow(moint)-1),1]$BBUSD1M - mean_moint

reg <- lm(y1 ~ lagy1)

summary(reg)

phi <- reg$coefficients[[2]]

FW_rates <- array(NA, dim = 12)

for (i in 2:12){
  FW_rates[i] = mean_moint + phi**i*(moint[[80,1]]-mean_moint)

}

FW_rates[1] = moint[[80,1]]

plot(FW_rates, type = "l")


dim(moint)

plot(FW_rates, type = "l")

Yield_rates <- array(NA, dim = 12)

for(i in 1:12){
  Yield_rates[i] = 1/i*sum(FW_rates[1:i])

}

plot(Yield_rates, type = "l")
lines(FW_rates)
lines(moint_80)


lines(moint[80,], col = "red")

Yield_rates
moint_80 <- array(NA, dim = 12)

for (i in 1:12){
  moint_80[i] = moint[[80,i]]
}



##### Discrete time vasicek ##### 

# y1t+1 - y1t = xt+1 - xt

deltaY <- diff(moint$BBUSD1M)


reg <- lm(deltaY ~ moint[1:(nrow(moint)-1),1]$BBUSD1M)

summary(reg)

b <- reg$coefficients[[2]]
a <- reg$coefficients[[1]]

phi = b + 1

sigma2 <- var(reg$residuals)

forward_rates <- 1/12*moint$BBUSD12 - 1/11*(moint$BBUSD11)

forward_rates <- -11*moint$BBUSD11 +12*(moint$BBUSD12)

plot(forward_rates)

y_f_diff <- mean(forward_rates - moint$BBUSD1M)

beta_func <- function(beta){
  res <- -(1/(1-phi)**2*sigma2/2-(beta/(1-phi))*sigma2) + y_f_diff
  return(res)
  
}


beta <- uniroot(beta_func, lower = -100, upper = 100)$root

x_t <- moint$BBUSD1M + 0.5*beta**2*sigma2

mu = mean(moint$BBUSD1M) + 0.5*beta**2*sigma2

forward_model <- array(NA, dim = dim(moint))


for (i in 1:ncol(forward_model)){
  t = i/12
  forward_model[,i] = mu - (beta + ((1-phi**t)/(1-phi)))**2*sigma2/2+phi**t*(x_t-mu)
  
}

for (i in 1:ncol(forward_model)){
  t = i
  forward_model[,i] = mu - (beta+1/(1-phi))**2*sigma2/2
  
}

y_forecast <- array(NA, dim = (nrow(moint)-1))

for (i in 1:length(y_forecast)){
  y_forecast[i] = (1-phi)*((mu - 0.5*beta**2*sigma2)-moint[[i,1]]) + moint[[i,1]]
}

plot(y_forecast[1:25], type = "l")
lines(moint[2:26,1]$BBUSD1M, col = "red")


forward_model[,12]

forward_rates

forward_rates


plot(y1)

forward_model[1,]

x_t 

mu


library(maRketSim)

library(remotes)
install_github("cran/maRketSim")

sim <- vasicek.discrete(8, 0.1, 1, 0.97, 7,25)

sim$rates



?vasicek.discrete



### Simple Term Structure Model FINAL ARR DATA ###

mean_moint = mean(Final_arr[,2])

y1 <- Final_arr[2:(nrow(Final_arr)),2] - mean_moint

lagy1 <- Final_arr[1:(nrow(Final_arr)-1),2] - mean_moint

reg <- lm(y1 ~ lagy1)

summary(reg)

phi <- reg$coefficients[[2]]


FW_rates <- array(NA, dim = 12)


mu - (beta + 1/(1-phi))**2*sigma2/2


for (i in 2:12){
  FW_rates[i] = mean_moint + phi**i*(moint[[80,1]]-mean_moint)
  
}

FW_rates[1] = moint[[80,1]]


dim(moint)

plot(FW_rates, type = "l")

Yield_rates <- array(NA, dim = 12)

for(i in 1:12){
  Yield_rates[i] = 1/i*sum(FW_rates[1:i])
  
}

plot(Yield_rates, type = "l")
lines(FW_rates)
lines(moint_80)


lines(moint[80,], col = "red")

Yield_rates
moint_80 <- array(NA, dim = 12)

for (i in 1:12){
  moint_80[i] = moint[[80,i]]
}



##### Discrete time vasicek ##### 

# y1t+1 - y1t = xt+1 - xt

deltaY <- diff(Final_arr[,2])


reg <- lm(deltaY ~ Final_arr[1:(nrow(Final_arr)-1),2])

summary(reg)

b <- reg$coefficients[[2]]
a <- reg$coefficients[[1]]

phi = b + 1
phi
sigma2 <- var(reg$residuals)

sigma2

forward_rates <- 1/12*moint$BBUSD12 - 1/11*(moint$BBUSD11)

forward_rates <- -60*Final_arr[,25] +120*(Final_arr[,31])

forward_rates

Maturities
plot(forward_rates)

y_f_diff <- mean(forward_rates - Final_arr[,2])

beta_func <- function(beta){
  res <- -(1/(1-phi)**2*sigma2/2-(beta/(1-phi))*sigma2) + y_f_diff
  return(res)
  
}


beta <- uniroot(beta_func, lower = -1000, upper = 1000)$root

beta
x_t <- Final_arr[,2] + 0.5*beta**2*sigma2

mu = mean(Final_arr[,2]) + 0.5*beta**2*sigma2

forward_model <- array(NA, dim = dim(Final_arr[,2:31]))

x_t
dim(Final_arr)



length(x_t)
ncol(forward_model)
for (i in 1:ncol(forward_model)){
  t = Maturities[i+1]*12
  forward_model[,i] = mu - (beta + ((1-phi**t)/(1-phi)))**2*sigma2/2+phi**t*(x_t-mu)
  
}

plot(forward_model[60,1:10], x = Maturities[1:10], type = "l")

Final_arr[60,2:11]


for (i in 1:ncol(forward_model)){
  t = i
  forward_model[,i] = mu - (beta+1/(1-phi))**2*sigma2/2
  
}

y_forecast <- array(NA, dim = (nrow(Final_arr)-1))

for (i in 1:length(y_forecast)){
  y_forecast[i] = (1-phi)*((mu - 0.5*beta**2*sigma2)-Final_arr[[i,1]]) + Final_arr[[i,1]]
}
plot(y_forecast[1:25], ylim = c(0.8,3))
lines(Final_arr[2:26,2])





