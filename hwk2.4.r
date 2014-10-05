##Q 4##
#Inputs: P is the cdf vector, A is numbers from 1:10^6, N is the desired number of A's

A <- seq(1:10^6)
CDF <- pgeom(A,10^-5)
N <- 10

discrete <- function(A,CDF,N) {
  M <- length(A)
  U <- runif(N)
  X <- numeric(N)
  P <- numeric(M)
  for(i in 1:M){
    if (i == 1) {
      P[i] = CDF[1]
      W <- which(U<= P[i])
    } else {
      P[i]<-CDF[i]
      W <- which(P[i-1]<U & U<= P[i])
    }
    X[W]=A[i]
    }
  return(X)
  }

discrete(A,CDF,N)

#a = a list of numbers
#p = the probabilities of drawing each number; these sum to 1
a <- 1:10^6
p <- 1e-5 *(.99999)^(a-1)

rvsampler <- function(a,p,n) {
  x <- runif(1)
  y <- numeric(n)
  P <- numeric(length(a))
  for(i in 1:length(a)) {
    if(i==1) {
      P[i] = p[1]
      w <- which(x <= P[i])
    } else {
      p[i] <- sum(p[1:i])
      w <- which(P[i-1]<x & x<=P[i])
    }
    y[w] <- a[i]
    }
  return(y)
  }


rvsampler(a,p,1)
