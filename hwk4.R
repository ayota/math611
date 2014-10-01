###2A : Write a function that samples from X ~ QY

#generates N RVs X where X ~ QY

xsampler <- function(N) {
  #enter the cov matrix sigma and find its eigenvalues/eigenvectors
  sigma <- matrix(c(2, 1, 1, 10), nrow=2,ncol=2,byrow =TRUE)
  sigma.eig <- eigen(sigma)
  
  #Q is the matrix of eigenvectors of sigma
  Q <- sigma.eig$vectors
  
  #generate Y's with distritbutions N(0, sigma.sq = eigenvalues)
  set.seed(100)
  y.pairs <- cbind(rnorm(N,mean = 0,sd=sqrt(sigma.eig$values[1])), rnorm(N,mean = 0,sd=sqrt(sigma.eig$values[2])))
  
 #
  X <- apply(y.pairs, 1, function(x) Q %*% x)
  
  
  return(X)
}

x <- xsampler(5)

pvalue <- function (N) {
  x <- xsampler(N)
  p.value <- sum(ifelse(((10*x[1,]^2)/19 - (2*x[1,]*x[2,])/19 + (2*x[2,]^2)/19) > 65.789, 1, 0))/N
  return(p.value)
}

pvalue(100)
#this is 0

##3a: K-means problem##

kmeanest <- function(data, means) {
  
  
}


