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

histsamples <- hist(samples, breaks = 10)

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

test <- ascent()
test

##steepest ascent for two component gmm
##gmm of two normal RVs: N(mu1,var1) w/ prob of p1 & N(mu2,var2) w/ prob p2


##initialize the function
#params to optimize: mu1, mu2, var1, var2, p1, p2
#because p1 + p2 = 1, we define p1 = p, p2 = 1 - p

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
      step <- step/2
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

testascent <- ascent(starts = c(0,0,1,1,0), target = .05, step = .0001)

testascent
debug(ascent)

