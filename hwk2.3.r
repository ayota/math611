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
