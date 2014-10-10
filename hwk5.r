###PROBLEM 2###


###PROBLEM 3###
###PART A###

params <- data.frame(rho = .1, time = 100)

markov <- function(params = data.frame(rho = .1, time = 5)) {
  P <- matrix(c(1-3*params$rho, params$rho, 2*params$rho, params$rho, 1 - 5*params$rho, 4*params$rho, 2*params$rho, 4*params$rho, 1-6*params$rho),3,3)
  x <- 1
  i <- 0
  path <- data.frame(time = i, X = x)
  
  while (i <= params$time) {
    i <- i + 1
    q <- runif(1)
    if (q < P[x,1]) { 
      x.new = 1 } 
    else if (q >= P[x,1] && q < (P[x,1] + P[x,2])) { 
      x.new = 2} 
    else { x.new = 3}
    x <- as.numeric(x.new)
    path <- rbind(path, data.frame(time = i, X = x.new))
  }
  return(path)
}

library(ggplot2)
path <- markov(params)
qplot(data=path, x=time, y=X, geom="step")

###PART B###
##part i: pi's##
rhos <- data.frame(rho1 <- .1, rho2 <- .01, rho3 <- .0001)

eig.method <- function (rho) {
  P <- matrix(c(1-3*params$rho, params$rho, 2*params$rho, params$rho, 1 - 5*params$rho, 4*params$rho, 2*params$rho, 4*params$rho, 1-6*params$rho),3,3)  
  Q <- eigen(P)
  pi <- sapply(Q$vectors[,which.max(Q$values)], function (x) x/sum(Q$vectors[,which.max(Q$values)]))
  return(pi)
}

for (i in 1:3) {
  pi <- eig.method(rho[1,i])
  print(pi)
}

##part ii: monte carlo##
params1 <- data.frame(rho = .1, time = 1000)
params2 <- data.frame(rho = .01, time = 1000)
params3 <- data.frame(rho = .0001, time = 1000)

monte.carlo <- function (N, params) {
  results <- data.frame(X1 = 0, X2 = 0, X3 =0)
  for (i in 1:N) {
  sample <- markov(params)
  x.1 <- sample[which(sample$X == 1),]
  x.2 <- sample[which(sample$X == 2),]
  x.3 <- sample[which(sample$X == 3),]
  results <- rbind(results, data.frame(X1 = nrow(x.1)/params$time, X2=nrow(x.2)/params$time, X3=nrow(x.3)/params$time))
  }
  final <- data.frame(X1 = mean(results$X1), X2 = mean(results$X2), X3 = mean(results$X3))
  return(final)
}

#run these later
monte.carlo(10000, params1)
monte.carlo(10000, params2)
monte.carlo(10000, params3)

##part iii: using P^t##
library(expm)

params1 <- data.frame(rho = .1, time = 1000)
params2 <- data.frame(rho = .01, time = 1000)
params3 <- data.frame(rho = .0001, time = 100000)

P.t <- function(params = data.frame(rho = .1, time = 5)) {
    P <- matrix(c(1-3*params$rho, params$rho, 2*params$rho, params$rho, 1 - 5*params$rho, 4*params$rho, 2*params$rho, 4*params$rho, 1-6*params$rho),3,3)  
    P.t <- P%^%params$time
    return(P.t)
}

P.t(params1)
P.t(params2)
P.t(params3)
