###PROBLEM 2b###

#set starting parameters
params <- data.frame(alpha = .2, beta =.5, mu = 5, lambda = 10)

loglike <- function(n = hiv$frequency, i = hiv$encounters, params = data.frame(alpha = .01, beta =.02, mu = .7, lambda = .3)) {
  #using dpois() and 0! is 1, the factorials are integrated into the loglikelihood function for better readibility and to reduce human error
  pi <- ifelse(i == 0, (params$alpha + params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda)),(params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda)))
  fxn <- sum(n * log(pi))
  return(fxn)
}

EM.hiv <- function(n = hiv$frequency, i=hiv$encounters, params = data.frame(alpha = .01, beta =.02, mu = .7, lambda = .3)) {
  iter <- 2
  L.theta <- loglike(params)
  path <- rbind(list(iter = 0, alpha = 0, beta = 0, mu = 0, lambda = 0, L.theta = 0),list(iter=0, alpha = params$alpha, beta = params$beta, mu = params$mu, lambda = params$lambda, L.theta = L.theta))
  
  #the function runs until L.theta is no longer increasing
  while (abs(path[iter,]$L.theta - path[iter-1,]$L.theta) > 0) {
    
    #expectation step
    z <- params$alpha/(params$alpha + params$beta*dpois(0, params$mu) + (1-params$alpha-params$beta)*dpois(0,params$lambda))
    t <- ifelse(i == 0, ((params$beta*dpois(i, params$mu)) /(params$alpha + params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))),((params$beta*dpois(i, params$mu))/(params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))))
    p <- ifelse(i == 0, (((1 - params$alpha - params$beta)*dpois(i, params$lambda)) /(params$alpha + params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))),(((1 - params$alpha - params$beta)*dpois(i, params$lambda))/(params$beta*dpois(i, params$mu) + (1-params$alpha-params$beta)*dpois(i,params$lambda))))
    
    
    #maximization step
    params$alpha <- n[1]*z/sum(n)
    params$beta <- sum(n*t/sum(n))
    params$mu <- sum(i*n*t)/sum(n*t)
    params$lambda <- sum(i*n*p)/sum(n*p)
    
    L.theta <- loglike(params)
    iter <- iter + 1
    path <- rbind(path,list(iter = (iter-2), alpha = params$alpha, beta = params$beta, mu = params$mu, lambda = params$lambda, L.theta = L.theta))
  }
  
  return(as.list(path[2:nrow(path),]))
}

path <- EM.hiv(params)

library(ggplot2)
qplot(x=as.numeric(path[,1]), y=as.numeric(path[,6]), xlab = "Iteration", ylab="L.theta")

###PROBLEM 3###
###PART A###

#set starting parameters
params <- data.frame(rho = .1, time = 100)

markov <- function(params = data.frame(rho = .1, time = 5)) {
  #enter probability matrix P
  P <- matrix(c(1-3*params$rho, params$rho, 2*params$rho, params$rho, 1 - 5*params$rho, 4*params$rho, 2*params$rho, 4*params$rho, 1-6*params$rho),3,3)
  
  #start x in space 1 at time (i) = 0
  x <- 1
  i <- 0
  path <- data.frame(time = i, X = x)
  
  while (i <= params$time) {
    i <- i + 1
    
    q <- runif(1)
    
    #select path of X based on uniform rv q
    if (q < P[x,1]) {x.new = 1} 
    else if (q >= P[x,1] && q < (P[x,1] + P[x,2])) {x.new = 2} 
    else {x.new = 3}
    
    x <- as.numeric(x.new)
    path <- rbind(path, data.frame(time = i, X = x.new))
  }
  return(path)
}

library(ggplot2)
path <- markov(params)
qplot(data=path, x=time, y=X, geom="step")

###PART B###
##PART i: Pi's##

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

##PART ii: MONTE CARLO##
params1 <- data.frame(rho = .1, time = 1000)
params2 <- data.frame(rho = .01, time = 10000)
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

library(compiler)
cmpfun(markov)
cmpfun(monte.carlo)

monte.carlo(1000, params1)
monte.carlo(1000, params2)
monte.carlo(1000, params3)

##PART iii: USING P^t##
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
