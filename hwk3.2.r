###part a###
##this function generates a gaussian mixture model.
#Input: N, the number of RVs desired; a set of parameters for each normal, in the form c(p1,p2,mu1,mu2,var1,var2)
#Output: A vector of Xi's

params <- c(0, 3, 1, 4,.7)

GMM <- function(N = 100, params = c(0, 3, 1, 4,.7)) {
  #generate 100 uniform RVs to represent test for normal 1 or normal 2
  set.seed(100)
  X <- as.vector(runif(N))
  
  #now apply the condition if uniform < p1, we assign normal 1, if greater assign normal 2
  set.seed(100)
  sapply(X, FUN = function(x) {if (x < params[5]) rnorm(1, mean = params[1], sd = sqrt(params[3])) else rnorm(1, mean = params[2], sd = sqrt(params[4]))})
  
  return(X)
}

samples <- GMM(100, params)

histsamples <- hist(samples, breaks = 10)

###part b, ii: steepest ascent algorithm###
#this finds the mle estimator for a given function

##test function for f(x1,x2) = -(10-x1)^2 - (5-x2)^2
##initialize the function
#Input: x is a vector of param x1 and x2
#Output: single value y
testfunction <- function(x) {
  y <- -(10-x[1])^2 - (5 - x[2])^2
  return(y)
}

##function to calculate the gradient
#Input: x is a vector of params x1 and x2
#Output: a vector of two values representing the partial for x1 and x2

gradient <- function(x) {
  df.dx1 <- (20 - 2*x[1])
  df.dx2 <- (10 - 2*x[2])
  grad <- c(df.dx1,df.dx2)
  return(grad)
}

##function to calculate the norm of the gradient
#Input: a gradient vector
#Output: single numeric value
norm <- function(grad) {
  gradnorm = sqrt(sum(grad^2))
  return(gradnorm)
}

##steepest ascent algorithm for test function
##inputs:
#starts <- start values for all params, as a vector
#target <- some sufficiently small value
#step <- how far each step moves
##output: a path for x1/x2 to maximizers, last value is the max

testascent <- function(starts = c(0,0), target = .01, step = .1) {
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

test <- testascent()

##steepest ascent for two component gmm
##gmm of two normal RVs: N(mu1,var1) w/ prob of p1 & N(mu2,var2) w/ prob p2

##log likelihood function
#params to optimize: mu1, mu2, var1, var2, p1, p2
#because p1 + p2 = 1, we define p1 = p, p2 = 1 - p
#Input: Vector of params in form (mu1,mu2,var1,var2,p); vector of Xi's
#loglikelihood estimate L(theta)

loglike <- function(params, samples = samples) {
  #(this is not necessary, it is just so I don't get confused)
  mu1 <- params[1]
  mu2 <- params[2]
  var1 <- params[3]
  var2 <- params[4]
  P1 <- params[5]
  P2 <- 1-P1
  
  #define the log-likelihood function without the summation
  loglike.fxn <- function(x) {
    Y <- log( P1*( (sqrt(2*pi*var1))^(-1) * exp( -(x-mu1)^(2) * (2*var1)^(-1)) ) + P2*( (sqrt(2*pi*var2))^(-1) * exp( -(x-mu2)^(2) * (2*var2)^(-1)) ) )
    return(Y)
  }
  
  #apply the new function for all the Xi's
  L.theta <- sum(sapply(samples, FUN = loglike.fxn))
  
  #return the sum (aka L(theta))
  return(L.theta)
}

##function to calculate the gradient
#samples is a vector of 100 x's calculated using GMM fxn
#params is a vector of 5 initial values
#df.dx represents derivative w/ respect to x
#the partial derivative fxns are compiled to make things go faster

#Partial derivative finding function
#inputs: single x value
#output: vector of partials df.mu1, df.mu2, df.var1, df.var2, df.P

#loading the compiler package to hopefully make this run faster
library(compiler)

#numerator-finding fxns
mu.numerator <- function(P,mu,var,x) { 
    numerator <- P*( ((2*pi*var)^-.5) * exp(-(x-mu)^2/(2*var)) * 2*(x-mu))
    return (numerator)
  }




var.numerator <- function(P,mu,var, x) { 
    numerator <- P/(sqrt(2*pi)*var^4)*exp(-(x-mu)^2/(2*var))*(var^2.5*(-x+mu)-.5*var^2.5 +.5*var^1.5*(x-mu)^2)
    return(numerator)
  }


p.numerator <- function(mu,var,x) {
    numerator <- (2*pi*var)^(-.5)*exp(-(x-mu)^2 * (2*var)^(-1))
    return(numerator)
  }

partializer <- function(x=1,params = c(0,3,1,4,.7)) {
  #initialize the params
  mu1 <- params[1]
  mu2 <- params[2]
  var1 <- params[3]
  var2 <- params[4]
  P1 <- params[5]
  P2 <- 1-P1
  #denominator for all partials
  R <- (P1*(1/sqrt(2*pi*var1))*exp((-(x-mu1)^2)/(2*var1)) + (1-P1)*(1/sqrt(2*pi*var2))*exp((-(x-mu2)^2)/(2*var2)))
  
  #calculate all the partials
  df.dmu1 <- mu.numerator(P1,mu1,var1,x)/R
  df.dmu2 <- mu.numerator(P2,mu2,var2,x)/R
  df.dvar1 <- var.numerator(P1,mu1,var1,x)/R
  df.dvar2 <- var.numerator(P2,mu2,var2,x)/R
  df.dp <- (p.numerator(mu1,var1,x)/R - p.numerator(mu2,var2,x)/R)
  grad <- c(df.dmu1,df.dmu2, df.dvar1,df.dvar2,df.dp)
  return(grad)  
}

gradient <- function(x.new = c(0,3,1,4,.3), X = samples) {
  #put the params in the right spots
  
  #find the partials for all Xi's
  grad <-vapply(X, FUN = partializer, params = x.new, c(mu1=0, mu2=0, var1=0, var2=0, P=0))
  
  #sum the partials for each parameter
  grad.sums <- apply(grad, 1, sum)
  
  #return vector of the partials for each parameter
  return(grad.sums)
}


#test to make sure it works
gradparams1 <- c(0, 3, 1, 4,.3)
testgrad3 <- gradient(gradparams1,samples)


c(0, 3, 1, 4,.3)
gradparams2 <- c(0,1,3,4,.7)
testgrad4 <- gradient(gradparams2,samples)

##function to calculate the norm of the gradient
norm <- function(grad) {
  gradnorm = sqrt(sum(grad^2))
  return(gradnorm)
}

#making sure it works
testnorm <- norm(testgrad)

##steepest ascent algorithm for test function
##inputs:
#starts <- start values for all params, as a vector
#target <- some sufficiently small value
#step <- how far each step moves

##output: a path for x1/x2 to maximizers
samples <- GMM(100, params)

ascent <- function(starts = c(-.5,2,.01,3,.4), target = .01, step = .1) {
  #initializing
  params.new <- starts
  grad.new <- gradient(params.new, samples)
  path <- matrix(0,1,5)
  norms <- matrix(0,1,1)
  count <- 0
  grad.old <- 0
  
  #loop will continue testing norm of gradient until it is smaller than the target
  while (abs(norm(grad.new) - norm(grad.old)) > target) {
    
    #test to make sure gradient is continuing to get smaller, if not, we decrease step
    #if (norm(grad.new) %in% norms) {
     # step <- step/10
    #}
    
    path <- rbind(path, params.new)  #add new x to previous path
    
    #make room for new x/gradient of new x, and store old gradients in matrix of past values
    params.old <- params.new 
    grad.old <- grad.new 
    norms <- rbind(norms,norm(grad.old))
    
    #calculate the direction of the new x
    dir <- grad.old/norm(grad.old)
    
    #get new x and gradient
    params.new <- params.old + dir*step
    #x.new[5] <- ifelse(abs(x.new[5]) >= 1, abs(x.new[5]/10), abs(x.new[5]))
    
    #count <- count + 1
    #if(count > 10000) {break}
    grad.new <- gradient(params.new, samples)
  }
  return(path)
}

#compile all the functions
cmpfun(loglike)
cmpfun(mu.numerator)
cmpfun(var.numerator)
cmpfun(p.numerator)
cmpfun(partializer)
cmpfun(gradient)
cmpfun(norm)
cmpfun(ascent)

samples <- GMM(1000, params)

testascent <- ascent(starts = c(-.1,2.8,.8,3.8,.5), target = .001, step = .005)


#let's see how we did
logs <- NULL

for(i in 2:length(testascent[,1])) {
  log <- loglike(testascent[i,],samples) 
  logs <- c(logs,log) 
}

plot(logs)
tail(testascent)

mean(samples[1:70])
mean(samples[30:100])
sd(samples[1:70])
sd(samples[30:100])
