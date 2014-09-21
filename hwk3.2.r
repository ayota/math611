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
#this doesn't work
#fxn <- function to be tested (as a function)
#starts <- start values for all params
#target <- some sufficiently small value
#step <- how far each step moves

ascent <- function(fxn) {
  #initializing
  x.new <- c(0,0)
  target <- .01
  step <- 1
  grad.new <- gradient(x.new)
  path <- NULL
  
  while (norm(grad.new) > target) {
    path <- c(path, x.new)
    x.old <- x.new
    grad.old <- grad.new
    dir <- grad.old/norm(grad.old)
    move <- dir*step
    x.new <- x.old + move
    grad.new <- gradient(x.new)
    print(x.new)
  }

  return(c(path,x.new))
}
