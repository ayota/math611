library(compiler)
library(microbenchmark)
library(sm)

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
  results[,2] <- sapply(results[,1], FUN = inverse, step=.005,start=-5)
  results[,3] <- sapply(results[,1], FUN =qnorm)
  return(results)
}

sample <- sampler(10000)

hist(sample[,2], col="red")
hist(sample[,3], col="blue", add=T)
plot(density(sample[,2]),col="blue")
plot(density(sample[,3]), col="red")

class(sample[,2])

qqnorm(sample[,2])
qqline(sample[,2],col="red")

###Question 3c###
#Overly favors postive Q over negative Q
#C is the constant g(x) is multiplied by
#N is the desired number of samples
genNorm <- function(C, N){
  nvars <- c()
  
  while(length(nvars) < N){
    u <- C*runif(1); 
    expvars <- runif(1)
    
    if (expvars > .5) {
      y = rexp(1,rate=1)
    } else {
      y = -rexp(1, rate=1)
    }
    
    ##redo the math here
    if(u < exp((y^2/2)-y)/sqrt(2*pi)) {
      nvars <- c(nvars,y)
    }

  }
  return(nvars)
}

vars <- genNorm(2,10000)
rnormals <- rnorm(10000)

hist(vars, col = "red")
hist(rnormals, add=T)

plot(density(vars), col="red")
plot(density(rnormals))

qqnorm(vars)
qqline(vars, col="red")

