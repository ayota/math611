library(compiler)
library(microbenchmark)

##611 HWK 2##
qnorm(.95)
qnorm(.975)
qnorm(.99)

###HWK 2, 3A###
#p is the desired probability we are looking for x for
inverse <- function(p, step, start) {
  #initialize all the things!
  norm<-function(x){1/sqrt(2*pi)*exp(-x^2/2)} 
  i <- start + step
  area <- 0
  
  #integrate function
  while (area < p) {
    integrand <- integrate(norm,lower=-Inf,upper=i)
    area <- integrand$value
    xvalue <- i
    i <- i + step
  }
  
  return(xvalue)
}

##timing##
compare <- microbenchmark(inverse(.975), qnorm(.975), times = 10)

###Question 3b###

sampler <- function(N) {
  #Set seed and generate vector of N random RVs
  set.seed(100)
  results <- matrix(0,N,3)
  results[,1] <- as.vector(runif(N,min=0,max=1))
  results[,2] <- sapply(results[,1], FUN = inverse, step=.01,start=-3)
  results[,3] <- sapply(results[,1], FUN =qnorm)
  return(results)
}

qqnorm(results[,2])
qqline(results[,2],col="red")

###Question 3c###
#Overly favors postive Q over negative Q
#C is the constant g(x) is multiplied by
#N is the desired number of samples
genNorm <- function(C, N){
  nvars <- c()
  while(length(nvars) < N){
    u = C*runif(1); 
    y = rexp(1,rate=1)
    q = -rexp(1, rate=1)
    if(u < exp(y-(y^2/2))/sqrt(2*pi)) {
      nvars <- c(nvars,y)
    }
    
    if(u < exp(q-(q^2/2))/sqrt(2*pi)) {
      nvars <- c(nvars,q)
    }
  }
  return(nvars)
}

vars <- genNorm(1,10000)
nvars <- rnorm(10000)

hist(vars)
hist(nvars)
qqnorm(vars)
qqline(vars, col="red")
