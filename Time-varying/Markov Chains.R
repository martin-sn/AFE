#### Stationary distribution of a Markov chain ####


# The stationary distribution of a Markov chain is equal to the
# limiting distribution of the Markov chain

# There are two ways to find it. We can raise the matrix to a high power As seen below. 

dat = c(0.97, 0.03, 0.01, 0.99)
p = matrix(dat,nrow = 2, byrow = TRUE)

library(expm)

p%^% 10000

# Or we can solve it numerically. 

# Delta'(I_j - Gamma + U) = 1'
# Where I_j is the JxJ identity matrix
# U is a JxJ matrix of ones
# and 1 is a vector of ones. 
# For proof see Grimmet and Stikzaker (2001)


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

StaDist(p)



