
library(geometry)


f_SVHeston = function(T, n, theta, kappa, xi, rho){
  
  N = T*n + 1 
  dt = 1/n 
  
  # Initialization
  v0 = theta 
  
  
  dW = rnorm(N)*sqrt(dt) # increments of Brownian motion W
  v = rep(1,N)*v0
  dB = rho*dW + sqrt((1-rho^2))*rnorm(N)*sqrt(dt)
  
  for (j in 2:N){
    v[j] = v[j-1] + kappa*(theta-v[j-1])*dt + xi*sqrt(v[j-1])*dB[j-1]
  }
  
  
  for (i in 1:length(v)){
    v[i] = max(v[i],0)
  }
  sigma = sqrt(v)
  X = cumsum(sigma*dW)

  
  
  out = list()

  out[["X"]] = X
  out[["sigma"]] = sigma
  out[["v"]] = v

  
  return(out)  
}




n = 2340
dt = 1/n
T = 10
theta = 0.04/250
kappa = 5/250
xi = 0.5/250
rho = -0.5

hest = f_SVHeston(T, n, theta, kappa, xi, rho)

X = hest[[1]]
sigma = hest[[2]]
plot(hest[[2]], type = "l")




X



a = 250*100


diff_X = diff(X)

X_matrix = matrix(diff_X, nrow = n, ncol = T)

RV = colSums(X_matrix^2)

diff_X

X

RV

Sigma_matrix = matrix(sigma[1:length(sigma)-1], nrow=n,ncol=T)

IV = colMeans(Sigma_matrix^2)


plot(a*hest[[2]]^2, type = "l")
points(x=n*(1:T),y=a*RV, pch = 4)
points(x=n*(1:T),y=a*IV, pch = 0)

# Realized Quarticity

RQ = n*colSums(X_matrix^4)
CI = 1.96*sqrt(2*RQ/3)/sqrt(n)

plot(x = (1:T), y = a*RV, pch = 0, ylim = c(3,6))
points(x=1:T,y=a*IV, pch = 1, col = "red")
points(x=1:T,y=a*(RV-CI), pch = 2)
points(x=1:T,y=a*(RV+CI), pch = 2)
legend("toplef", legend = c("RV", "IV", "CI"), bty = "n",
       lwd = 2, cex = 1.2, col = c("black", "red", "black"), lty = c(NA, NA, NA), pch = c(0, 1, 2))


