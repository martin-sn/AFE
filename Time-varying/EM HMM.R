###### Expectation Maximization for Hidden Markov Models #####

### Simulate T = 5000 observations from a HMM with J = 2 regimes and trasition probability matrix

dat = c(0.97, 0.03, 0.01, 0.99)
p = matrix(dat,nrow = 2, byrow = TRUE)


sim_hmm = function(prob_matrix, n, j){
  state = array(NA, n)
  
  state[1] = 1 
  
  for (i in 2:n){
    state[i] = sample(c(1:j),1, prob=prob_matrix[state[i-1],])
    
    
    
    
  }
  
  return(state)
  
}

sim = sim_hmm(p, 5000, 2)

plot(sim, type = "l")


mu_1 = -4
mu_2 = 5
sigma_1 = 1.5
sigma_2 = 0.7


######## EM Algorithm for HMM ######## 

### E-Step
E_step <- function(X, Mixture){
  
  
  
  
  
  
  
  
  
}


