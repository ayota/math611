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

gradient <- function(params = c(0,3,1,4,.3), samples) {
  mu1 <- params[1]
  mu2 <- params[2]
  var1 <- params[3]
  var2 <- params[4]
  P <- params[5]
  
  grad <-matrix(0, 1, 5)
  
  for (i in 1:length(samples)) {
  #define denominator for all partials
  R <- (P*(1/sqrt(2*pi*var1))*exp((-(samples[i]-mu1)^2)/(2*var1)) + (1-P)*(1/sqrt(2*pi*var2))*exp((-(samples[i]-mu2)^2)/(2*var2)))

  #numerators
  #mu1.numerator <- P*(1/sqrt(2*pi*var1))*((samples[i]-mu1)/var1)*exp((-(samples[i]-mu1)^2)/(2*var1))
  #mu2.numerator <- (1-P)*(1/sqrt(2*pi*var2))*((samples[i]-mu2)/var2)*exp((-(samples[i]-mu2)^2)/(2*var2))

  mu1.numerator <- P*((samples[i]-mu1)*exp((-(mu1-samples[i])^2)/(2*var1))*(1/(sqrt(2*pi)*var1^(3/2))))
  mu2.numerator <- (1-P)*((samples[i]-mu2)*exp((-(mu2-samples[i])^2)/(2*var2))*(1/(sqrt(2*pi)*var2^(3/2))))
  
  #numerators for the var's (need to simplify)
  var1.numerator <- P*((exp(-(samples[i]-mu1)^2/(2*var1)))*(((-1/2)*(2*pi)^(-1/2)*(var1)^(-3/2)) + (1/sqrt(2*pi*var1))*((samples[i]-mu1)^2/(2*var1^2))))
  var2.numerator <- (1-P)*((exp(-(samples[i]-mu2)^2/(2*var2)))*(((-1/2)*(2*pi)^(-1/2)*(var2)^(-3/2)) + (1/sqrt(2*pi*var2))*((samples[i]-mu2)^2/(2*var2^2))))
  
  #var1.numerator <- P*(exp((-(samples[i]-mu1)^2)/(2*var1))*(-var1 + (-samples[i] + mu1)^2)*(1/(2*sqrt(2*pi)*samples[i]^(5/2))))
  #var2.numerator <- (1-P)*(exp((-(samples[i]-mu2)^2)/(2*var2))*(-var2 + (-samples[i] + mu2)^2)*(1/(2*sqrt(2*pi)*samples[i]^(5/2))))
    
  #numerator for p
  p.numerator <- (1/sqrt(2*pi*var1))*exp(-(samples[i]-mu1)^2/(2*var1)) - (1/sqrt(2*pi*var2))*exp(-(samples[i]-mu2)^2/(2*var2))
  
  df.dmu1 <- mu1.numerator/R
  df.dmu2 <- mu2.numerator/R
  df.dvar1 <- var1.numerator/R
  df.dvar2 <- var2.numerator/R
  df.dp <- p.numerator/R
  grad <- rbind(grad,c(df.dmu1,df.dmu2, df.dvar1,df.dvar2,df.dp))
}
  final <- c(sum(grad[,1]),sum(grad[,2]),sum(grad[,3]),sum(grad[,4]),sum(grad[,5]))

  return(final)
}

debug(gradient)
samples <- GMM(100, p1, p2, mu1, mu2, var1, var2)
testgrad <- gradient(samples = samples[1:100])
testgrad2 <- gradient(samples = samples[1:100])
testgrad
testgrad2

##function to calculate the norm of the gradient
norm <- function(grad) {
  gradnorm = sqrt(sum(grad^2))
  return(gradnorm)
}

debug(norm)
testnorm <- norm(testgrad)

##steepest ascent algorithm for test function
##inputs:
#starts <- start values for all params, as a vector
#target <- some sufficiently small value
#step <- how far each step moves

##output: a path for x1/x2 to maximizers
samples <- GMM(100, p1, p2, mu1, mu2, var1, var2)

ascent <- function(starts = c(-.5,2,0,3,.4), target = .01, step = .1) {
  #initializing
  x.new <- starts
  grad.new <- gradient(x.new, samples[1:100])
  path <- matrix(0,1,5)
  norms <- matrix(0,1,1)
  
  #loop will continue testing norm of gradient until it is smaller than the target
  while (norm(grad.new) > target) {
    
    #test to make sure gradient is continuing to get smaller, if not, we decrease step
    if (norm(grad.new) %in% norms) {
      step <- step/10
    }
    
    path <- rbind(path, x.new)  #add new x to previous path
    
    #make room for new x/gradient of new x, and store old gradients in matrix of past values
    x.old <- x.new 
    grad.old <- grad.new 
    norms <- rbind(norms,norm(grad.old))
    
    #calculate the direction of the new x
    dir <- grad.old/norm(grad.old)
    
    #get new x and gradient
    x.new <- x.old + dir*step
    #x.new[5] <- ifelse(abs(x.new[5]) >= 1, abs(x.new[5]/10), abs(x.new[5]))
    grad.new <- gradient(x.new, samples[1:100])
    print(x.new)
    print(norm(grad.new))
  }
  
  return(rbind(path,x.new))
}

testascent <- ascent(starts = c(-.1,2.8,.8,3.8,.5), target = .05, step = .0001)

fhtestascent
debug(ascent)

