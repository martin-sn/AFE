#### Preaveraging estimator
# Y: Observed prices
# K: Pre-averaged window, i.e. k_n
# n: Sample size 




f_preav = function(Y,K){
  kernel = function(x){
    mini = min(x,1-x)
    return(mini)
    }
  
  #r = diff(Y)
  r = diff(Y)
  
  n = length(Y)
  
  # We first find pre-averaged returns, (e.e., /bar(y)_i/n)
  r_pa = 0
  for (i in 1:(K-1)){
    r_pa = r_pa + kernel(i/K)*r[i:(length(r)+i-(K-1))]
  #  print(r_pa[i])
  }
  idx = (0:K)/K
  
  g = array(data = NA, dim = length(idx))
  for (i in 1:length(idx)){
    g[i] = kernel(idx[i])

  }
  dg = diff(g)
  psi1 = sum(dg**2)*K
  psi2 = sum(g**2)/K
  
  P = sum((r_pa)**2)
  P = P*n/((n-K+2)*K*psi2)
  P = P-sum(r**2)*psi1/(2*psi2*K**2)
  
  return(P)
  
}

f_preav(X,10)

length(X)


R_P