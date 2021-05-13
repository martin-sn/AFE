nsim = 10000
n = 23400 # 1 sec sampling in 6.5 trading hours
theta = 1 # tuning parameter of pre-averaging

sigma0 = 0.04/250
kappa = 5/250
xi = 0.5/250
rho = -0.5
gamma = 0.5 # Noise to signal ratio

P = array(NA, dim = nsim)
IV = array(NA, dim = nsim)

for (i in 1:nsim){
  print(i)
  
  heston = f_SVHeston(1,n,sigma0, kappa,xi,rho)
  X = heston[["X"]]
  sigma = heston[["sigma"]]
  IV[i] = mean(sigma[1:(length(sigma)-1)]**2)
  
  
  ## Adding noise to prices
  omega = gamma*sqrt(IV[i]/n) #omega**2 is variance of noise
  noise = omega*rnorm((n+1)) # add iid gaussian npise
  Y = X + noise
  
  # Pre-averaging
  
  K = round(theta*sqrt(n))
  P[i] = f_preav(Y,K)
  
  out = list()
  
  out[["P"]] = P 
  out[["IV"]] = IV
}

K

bias = mean(P-IV)
bias

Rbias = mean(P/IV-1)
Rbias

MSE = mean((P-IV)**2)
MSE

1.866552977303856e-10

IV
P
