pi

arr = array(data = NA, dim = c(2,3))

arr

Gaussian_Dist <- function(x, mean,var){
  log_prob = -1/2 * log(2*pi*var) - 0.5*((x-mean)^2) /var
  return(exp(log_prob))
  
}

Gaussian_Dist(0,0,1)

library(quantmod)

SP = getSymbols("^GSPC", from = "2002-01-02", to = "2008-08-29", auto.assign = FALSE)

GSPC_rt <- as.numeric(diff(log(SP$GSPC.Adjusted))*100)[-1]

vY = GSPC_rt

iJ = 5

## Initialization
vMu     = mean(vY) * seq(0.7, 1.5, length.out = iJ)
vSigma2 = var(vY) * seq(0.7, 1.5, length.out = iJ)
vOmega  = rep(1/iJ, iJ)

init_gaus = list(vMu, vSigma2, vOmega)

init_gaus[[1]]

E_step = function(X, GaussianMixture, J){
  

  
  mean = GaussianMixture[[1]]
  var = GaussianMixture[[2]]
  prob = GaussianMixture[[3]]
  K = J 
  
  
  arr = array(NA, dim = c(n,K))
  
  LLK = 0
  
  for (i in 1:n){
    for (j in 1:K){
      arr[i,j] = prob[j]*dnorm(X[i],mean[j],var[j])
      
      total = sum(arr[i])
      
    }
    
    arr[i] = arr[i] / total
    
    LLK = LLK + log(total)
    
    
    
    
  }
  
  out = list()
  out[["post"]] = arr
  out[["llk"]] = LLK
  
  return(out)
}


E_step(vY, init_gaus, 5)



M_step = function(X, post, J){
  n = len(X)
  K = J
  
  mu_hat = array(NA, dim=J)
  var_hat = array(NA, dim = J)
  prob_hat = array(NA, dim = J)
  
  
  for (j in 1:K){
    mu_hat[j] = sum(X*post[,j]) / sum(post[,j])
    prob_hat[j] = 1/n*sum(post[,j])
    var_hat[j] = (sum((X-mu_hat[j])^2) * post[,j])/sum(post[j])
    

    
  }
  
  GaussianMixture = list(mu_hat, prob_hat,var_hat)
  
  return(GaussianMixture)
}




Run = function(){
  
  
  
  
  
}