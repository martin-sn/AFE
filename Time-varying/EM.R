##### EXPECTATION MAXIMIZATION ALGORITHM (GAUSSIAN DIST) ####

### GET DATA ###

library(quantmod)

SP = getSymbols("^GSPC", from = "2002-01-02", to = "2008-08-29", auto.assign = FALSE)

GSPC_rt <- as.numeric(diff(log(SP$GSPC.Adjusted))*100)[-1]


arr = array(data = NA, dim = c(2,3))

arr

# Dist 

Gaussian_Dist <- function(x, mean,var){
  log_prob = -1/2 * log(2*pi*var) - 0.5*((x-mean)^2) /var
  return(exp(log_prob))
  
}

Gaussian_Dist(0,0,1)

dnorm(0,0,1)





E_step = function(X, GaussianMixture, J){
  
  mean = GaussianMixture[[1]]
  var = GaussianMixture[[2]]
  prob = GaussianMixture[[3]]
  K = J 
  
  n = length(X)
  
  arr = array(NA, dim = c(n,K))
  
  LLK = 0
  
  for (i in 1:n){
    for (j in 1:K){
      arr[i,j] = prob[j]*dnorm(X[i],mean[j],sqrt(var[j]))
      
      total = sum(arr[i,])
      
    }
    
    arr[i,] = arr[i,] / total
    
    LLK = LLK + log(total)
    
    
    
    
  }
  
  out = list()
  out[["post"]] = arr
  out[["llk"]] = LLK
  
  return(out)
}


e = E_step(vY, init_gaus, 5)

post = e[["post"]]
?dnorm

M_step = function(X, post, J){
  n = length(X)
  K = J
  
  mu_hat = array(NA, dim=J)
  var_hat = array(NA, dim = J)
  prob_hat = array(NA, dim = J)
  
  
  for (j in 1:K){
    mu_hat[j] = sum(X*post[,j]) / sum(post[,j])
    prob_hat[j] = 1/n*sum(post[,j])
    var_hat[j] = sum((X-mu_hat[j])^2 * post[,j])/sum(post[,j])
    

    
  }
  
  GaussianMixture = list(mu_hat, prob_hat,var_hat)
  
  return(GaussianMixture)
}



M_step(vY,post,5)



Run = function(X,J){
  
  
  vY = X
  iJ = J
  
  ## Initialization
  vMu     = mean(vY) * seq(0.7, 1.5, length.out = iJ)
  vSigma2 = var(vY) * seq(0.7, 1.5, length.out = iJ)
  vOmega  = rep(1/iJ, iJ)
  
  loops = 0

  init_gaus = list(vMu, vSigma2, vOmega)
  
  
  llk_prev = -10000
  
  e = E_step(X, init_gaus, J)
  
  post = e[["post"]]
  llk_new = e[["llk"]]
  
  not_conv = TRUE
  
 # while (not_conv == TRUE){
  for (i in 1:3){
    mixture = M_step(X, post, J)
    llk_prev = llk_new
    e = E_step(X, mixture, J)
    post = e[["post"]]
    llk_new = e[["llk"]]
    loops = loops + 1
    
    
   # if (llk_new - llk_prev <= 0.000001*abs(llk_new)){
  #    not_conv = FALSE

   # }
  }

  
  out = list()
  out[["Mixture"]] = mixture
  out[["Post"]] = post
  out[["LLK"]] = llk_new
  out[["loops"]] = loops
  out[["conv"]] = not_conv
  out[["llk_prev"]] = llk_prev
  out[["llk_new"]] = llk_new
  return(out)
  
}

Run(GSPC_rt,5)

GSPC_rt

GSPC_rt




