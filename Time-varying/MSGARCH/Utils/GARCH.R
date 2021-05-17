###### ESTIMATE GARCH ###### 

likelihood_est_garch <- function(y,omega,alpha,beta){
  t = length(y)
  sigma2 = numeric(t)
  
  sigma2[1] = omega/(1.0-alpha-beta)
  
  for (i in 2:t){
    sig_t = sigma2[(i-1)]
    sigma2[i] = omega + alpha*y[(i-1)]^2 + beta*sig_t
  }
  llk = sum(dnorm(y,0,sqrt(sigma2), log = TRUE))
  
  return(llk)
}



####Objective function 

Obj_func_garch <- function(vPar, y){
  omega = vPar[1]
  alpha = vPar[2]
  beta = vPar[3]
  llk = likelihood_est_garch(y,omega,alpha,beta)
  return(-llk)
  
}



EstimateGARCH <- function(y){
  require(Rsolnp)
  alpha = 0.8 
  beta = 0.1
  omega = var(y) * (1 - alpha - beta)
  vPar = c(omega, alpha, beta)
  AB = c(alpha,beta)
  
  #New optimizer with fixed constraints 
  
  optimizer = solnp(vPar, fun = Obj_func_garch, ineqfun = function(vPar,...) {
    sum(vPar[2],vPar[3])
  }, ineqLB = 1e-12, ineqUB = 0.99999999, 
  LB = c(1e-12, 1e-12, 1e-12), UB = c(0.99999999, 0.99999999, 0.99999999),y=y)
  
  
  
  par = optimizer$par
  llk = -optimizer$value
  
  output = list()
  output[["Par"]] = par
  output[["llk2"]] = llk
  
  
  return(output)
}
