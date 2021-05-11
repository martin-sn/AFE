###### Numerical Estimation (Not the best way) ######

# Matrix LLK from slide 50 HMM
LLK_Matrix = function(p, y, mu_1, mu_2, s_1, s_2){
  Sta_dist <- p%^% 10000
  
  P_Y = matrix(c(dnorm(y[1], mu_1, s_1, log = TRUE), 0, 0, dnorm(y[1], mu_2, s_2, log = TRUE)), byrow = TRUE, nrow = 2)
  LLK = log(Sta_dist) + P_Y
  
  log_gamma = log(p)
  
  for (i in 2:length(y)){
    
    P_Y = matrix(c(dnorm(y[i], mu_1, s_1, log = TRUE), 0, 0, dnorm(y[i], mu_2, s_2, log = TRUE)), byrow = TRUE, nrow = 2)
    
    LLK = LLK + log_gamma + P_Y
    
  }
  
  return(LLK %*% matrix(1, nrow = 2)) # Not sure about this 1 
}


# Different approach
LLK_Matrix = function(p, y, mu_1, mu_2, s_1, s_2){
  Sta_dist <- p%^% 10000
  
  P_Y = matrix(c(dnorm(y[1], mu_1, s_1, log = FALSE), 0, 0, dnorm(y[1], mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)
  
  LLK = Sta_dist %*% P_Y
  
  #LOG_LLK = log(LLK)
  
  #log_gamma = log(p)
  
  for (i in 2:length(y)){
    
    P_Y = matrix(c(dnorm(y[i], mu_1, s_1, log = FALSE), 0, 0, dnorm(y[i], mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)
    
    LLK = LLK %*% p %*% P_Y
    
    #LOG_LLK = LOG_LLK + log(LLK)
    
  }
  
  return(LLK %*% array(1, 2)) # Not sure about this 1 
}



LLK_Matrix(p, y, mu_1, mu_2, s_1, s_2)

LLK_Matrix(p, y, -3, 2, 1, 1) # <- Using other means than what the data was generated with gives better results...





## Testing below.... ##

P_Y = matrix(c(dnorm(1, mu_1, s_1, log = FALSE), 0, 0, dnorm(1, mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)

max(P_Y)

P_Y

Sta_dist * P_Y



P_Y + log(p%^% 10000)


P_Y = matrix(c(dnorm(1, mu_1, s_1, log = FALSE), 0, 0, dnorm(1, mu_2, s_2, log = FALSE)), byrow = TRUE, nrow = 2)


log(p%^% 10000 * P_Y)

