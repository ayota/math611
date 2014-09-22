###part a###
##initialize all the things!##
p1 <- .7
p2 <- .3
mu1 <- 0
mu2 <- 3
var1 <- 1
var2 <- 4

GMM <- function(N = 100, p1 = .7, p2 = .3, mu1 = 0, mu2 = 3, var1 = 1, var2 = 4) {
  X <- matrix(0, N, 1)
  n1 = n2 = 0
  
  for(i in 1:N) {
    prob <- runif(1)
    if(prob < p1) {
      X[i, 1] <- rnorm(1, mean = mu1, sd = sqrt(var1) )
    } else {
      X[i, 1] <- rnorm(1, mean = mu2, sd = sqrt(var2) )
    }
  }
  return(X)
}

samples <- GMM(100, p1, p2, mu1, mu2, var1, var2)

histsmaples <- hist(samples, breaks = 10)

###part b, i###



###part b, ii: steepest ascent algorithm###

##test function for f(x1,x2) = -(10-x1)^2 - (5-x2)^2

##initialize the function
#x is a vector of param x1 and x2
testfunction <- function(x) {
  y <- -(10-x[1])^2 - (5 - x[2])^2
  return(y)
}

##function to calculate the gradient
#x is a vector of params x1 and x2
gradient <- function(x) {
  df.dx1 <- (20 - 2*x[1])
  df.dx2 <- (10 - 2*x[2])
  grad <- c(df.dx1,df.dx2)
  return(grad)
}

##function to calculate the norm of the gradient
norm <- function(grad) {
  gradnorm = sqrt(sum(grad^2))
  return(gradnorm)
}

##steepest ascent algorithm for test function
##inputs:
#starts <- start values for all params, as a vector
#target <- some sufficiently small value
#step <- how far each step moves

##output: a path for x1/x2 to maximizers

ascent <- function(starts = c(0,0), target = .01, step = .1) {
  #initializing
  x.new <- starts
  grad.new <- gradient(x.new)
  path <- matrix(0,1,2)
  gradients <- matrix(0,1,2)
  
  #loop will continue testing norm of gradient until it is smaller than the target
  while (norm(grad.new) > target) {
    
    #test to make sure gradient is continuing to get smaller, if not, we decrease step
    if (grad.new[1] %in% gradients[,1] && grad.new[2] %in% gradients[,2]) {
      step <- step/2
    }
    
    path <- rbind(path, x.new)  #add new x to previous path
    
    #make room for new x/gradient of new x, and store old gradients in matrix of past values
    x.old <- x.new 
    grad.old <- grad.new 
    gradients <- rbind(gradients,grad.old)
    
    #calculate the direction of the new x
    dir <- grad.old/norm(grad.old)
    
    #get new x and gradient
    x.new <- x.old + dir*step
    grad.new <- gradient(x.new)
  }

  return(rbind(path,x.new))
}

test <- ascent(testfunction)
test
