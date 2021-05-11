##### Simulate HMM ####

SimulateHMM = function(prob_matrix, n, j){
  state = array(NA, n)
  
  # We initialize the first state to be 1. 
  state[1] = 1 
  
  for (i in 2:n){
    state[i] = sample(c(1:j),1, prob=prob_matrix[state[i-1],])
    
  }
  
  return(state)
  
}


##### Generate Y for Gaussian HMM ###
# Accepts a HMM_Simulation and a Gaussian mixture as input. 
# The gaussian mixture must be a matrix where the first column are the means and the second column are the standard deviations. 

GenerateY_GausHMM <- function(HMM_Sim, mix){
  y = array(NA, length(HMM_Sim))
  J = nrow(mix)
  for (i in 1:length(HMM_Sim)){
    for (j in 1:J)
      if (HMM_Sim[i] == j){
        y[i]= rnorm(1, mix[j,1], mix[j,2])
      }
  }
  return(y)
}




### Computes the stationairy distrubtion of a Markov Chain ###
StaDist <- function(p){
  j = nrow(p)
  U = matrix(1, nrow=j, ncol=j)
  ONE = matrix(1, nrow=j, ncol=1)
  identity = diag(j)
  
  # We find the inverse of (I_j - Gamma + U)
  invs = solve(identity - p + U)
  
  # We multiply both sides of the equation with the inverse we just found
  # and then we get the stationary distribution (Delta')
  statio = t(ONE) %*% invs
  
  return(statio)
  
}

