G(t) = GDP at time t, where t is yearly quarters
GDP <- GDP[,2]

#apply formula 1 to GDP data
Y.data <- 100*(log(GDP[2:267]) - log(GDP[1:266]))

#functions to find epsilons, transition probabilities, and m-h ratio
epsilon <- function(Z=Z, phi.1 = params$phi.1, phi.2 = params$phi.2, phi.3 = params$phi.3, phi.4 = params$phi.4, sigma = params$sigma) {
  errors <- (Z[5:266] - phi.1*Z[4:265]- phi.2*Z[3:264] - phi.3*Z[2:263] - phi.4*Z[1:262])
  errors <- dnorm(errors, mean = 0, sd = sigma)
  return(errors)
}

trans.prob <- function(s.1, s.2, p = params$p, q = params$q) {
  if(s.1 == s.2) {
    prob <- ifelse(s.2 == 1, p, q)
  } else {
    prob <- ifelse(s.2 == 1, 1-p, 1-q)
  }
  return(prob)
}

nu <- function(eps, probs, sigma = params$sigma) {
  nu <- sum(log(c(eps, probs)))
  return(nu)
}

#initializing params
params <- data.frame("alpha.1" = 1.522, "alpha.0" = -.3577, 
                     "p" = .9049, "q" = .7550, "sigma" = .7690, "phi.1" = .014,
                     "phi.2" = -.058, "phi.3" = -.247, "phi.4" = -.213)

#initializing states
#1 = recession; 0 = OK
initial.mix <- sample(c(0,1), 266, replace = TRUE) #mix
initial.g <- rep(0,266) #all good times
initial.r <- rep(1,266) #all recession


#M-H RATION WRONG: NEEDS TO BE NU.J/NU.I
mh.sampler <- function(N=10000, initial=initial, params=params, Z=Z, Y.data=Y.data) {
  current <- initial
  Z.current <-  Y.data-params$alpha.1*current-params$alpha.0 #find Z's
  p.current <- sapply(2:266, function(x) trans.prob(s.1=current[x-1], s.2=current[x])) #transitions
  ep.current <- epsilon(Z=Z.current, phi.1=params$phi.1, phi.2=params$phi.2, phi.3=params$phi.3, phi.4=params$phi.4, sigma=params$sigma) #find error terms
  nu.i <- nu(eps=ep.current, probs=p.current, sigma=params$sigma)
  
  i <- 1 #counter
  #states <- matrix(0, 266, N) #this will store the states over the course of the simulation
  metric <- states <- matrix(0, 1, N) #the metric for judging detailed balance will be the % of time in recession (sum(current)/266)
  accept <- 0 #this will count the number of time the proposal was accepted over the course of the simulation
  
  while ( i <= N) {
    #make proposal
    flip <- sample(1:266, 1) #flip 1 state
    proposal <- current
    proposal[flip] <- abs(1-current[flip])
    #proposal <- sample(c(0,1), 266, replace = TRUE)
    
    #calculate all the parts of nu
    Z.proposal <-  Y.data-params$alpha.1*proposal-params$alpha.0 #find Z's
    p.proposal <- sapply(2:266, function(x) trans.prob(s.1=proposal[x-1], s.2=proposal[x])) #transitions
    ep.proposal <- epsilon(Z=Z.proposal, phi.1=params$phi.1, phi.2=params$phi.2, phi.3=params$phi.3, 
                           phi.4=params$phi.4, sigma=params$sigma) #find error terms
    nu.j <- nu(eps=ep.proposal, probs=p.proposal, sigma = params$sigma)
    
    
    U <- runif(1)
    
    if (U < min(1,exp(nu.j-nu.i))) {
      current <- proposal
      nu.i <- nu.j
      accept <- accept + 1
    }
    
    #states[,i] <- current
    metric[,i] <- sum(current)/266
    i <- i + 1
    if(i%%1000==0) {print(i)}
  }
  final <- list("% Acceptance" = accept/N, "% in Recession" = metric)
  return(final)
}

####CONFIGURING M-H SAMPLER
##PROPOSALS
#Acceptance Rate: .2484, up from less than .0001, when switched proposal to 1 flip instead of 10
#initial state is mix of recession and growth
test1 <- mh.sampler(N=10000, initial=initial.mix, params=params, Z=Z, Y.data=Y.data)

#Acceptance Rate: .001, when completely reshuffle data each time
#initial state is mix of recession and growth
test2 <- mh.sampler(N=10000, initial=initial.mix, params=params, Z=Z, Y.data=Y.data)

#Decision: Flip 1 at a time

##DETERMINING DETAILED BALANCE
library(ggplot2)
#Mixed initial state
test.mix <- mh.sampler(N=50000, initial=initial.mix, params=params, Z=Z, Y.data=Y.data)
test.mix <- as.numeric(unlist(test.mix[2]))
plot.mix <- qplot(1:50000,y = test.mix, geom = "line")

#All 0's
test.r <- mh.sampler(N=50000, initial=initial.g, params=params, Z=Z, Y.data=Y.data)
test.g <- as.numeric(unlist(test.r[2]))
plot.g <- qplot(1:50000,y = test.g, geom = "line")

#All 1's
test.g <- mh.sampler(N=50000, initial=initial.r, params=params, Z=Z, Y.data=Y.data)
test.r <- as.numeric(unlist(test.g[2]))
plot.r <- qplot(1:50000,y = test.r, geom = "line")

#INSET 3-PLOT FIGURE HERE

###M-H RATION WRONG: NEEDS TO BE NU.J/NU.I
##So all plots show a steady state really quickly, but just to be safe we will set burn time = 10,000 steps and run the simulation 1,000,000 times
mh.configure <- function(N=1000000, burn=10000, initial=initial, params=params, Z=Z, Y.data=Y.data) {
  current <- initial
  Z.current <-  Y.data-params$alpha.1*current-params$alpha.0 #find Z's
  p.current <- sapply(2:266, function(x) trans.prob(s.1=current[x-1], s.2=current[x])) #transitions
  ep.current <- epsilon(Z=Z.current, phi.1=params$phi.1, phi.2=params$phi.2, phi.3=params$phi.3, phi.4=params$phi.4, sigma=params$sigma) #find error terms
  nu.i <- nu(eps=ep.current, probs=p.current, sigma=params$sigma)
  
  i <- 1 #counter
  states <- matrix(0, 266, 1) #this will store the states over the course of the simulation in the form of a sum
  #metric <- states <- matrix(0, 1, N) #the metric for judging detailed balance will be the % of time in recession (sum(current)/266)
  accept <- 0 #this will count the number of time the proposal was accepted over the course of the simulation
  
  while ( i <= N) {
    #make proposal
    flip <- sample(1:266, 1) #flip 1 state
    proposal <- current
    proposal[flip] <- abs(1-current[flip])
    #proposal <- sample(c(0,1), 266, replace = TRUE)
    
    #calculate all the parts of nu
    Z.proposal <-  Y.data - params$alpha.1*proposal-params$alpha.0 #find Z's
    p.proposal <- sapply(2:266, function(x) trans.prob(s.1=proposal[x-1], s.2=proposal[x])) #transitions
    ep.proposal <- epsilon(Z=Z.proposal, phi.1=params$phi.1, phi.2=params$phi.2, phi.3=params$phi.3, 
                           phi.4=params$phi.4, sigma=params$sigma) #find error terms
    nu.j <- nu(eps=ep.proposal, probs=p.proposal, sigma = params$sigma)
    
    
    U <- runif(1)
    if (U < min(1,exp(nu.j-nu.i))) {
      current <- proposal
      nu.i <- nu.j
      accept <- accept + 1
    }
    
    if(i >= burn) { states <- states + current }
    #metric[,i] <- sum(current)/266
    i <- i + 1
    if(i%%1000==0) {print(i)}
  }
  final <- list("% Acceptance" = accept/N, "States" = states/N)
  return(final)
}

#test
test.config <- mh.configure(N=200000, burn=10000,initial=initial.r, params=params, Z=Z, Y.data=Y.data)
acceptance <- test.config[1]
test.config <- as.numeric(unlist(test.config[2]))
plot.config <- qplot(1:152,y = test.config[1:152], geom = "line")


##PART B: ADDING PRIORS
library(msm)
##All normal priors, specifically:
alpha.1.star <- rtnorm(1, mean=1.5, sd=1, lower = 0, upper = 4)
alpha.0.star <-rtnorm(1, mean=-.36, sd=.75, lower = -1, upper = 1)
p.star <- rtnorm(1, mean=.6, sd = .4, lower = 0, upper = 1)
q.star <- rtnorm(1, mean=.6, sd = .4, lower = 0, upper = 1)
sigma.star <- rtnorm(1, .8, sd=.1, lower = 0, upper = 4)
phi.1.star <- rtnorm(1, 0, sd=.5, lower = -2, upper = 2)
phi.2.star <- rtnorm(1, 0, sd=.5, lower = -2, upper = 2)
phi.3.star <- rtnorm(1, -.25, sd=.5, lower = -2, upper = 2)
phi.4.star <- rtnorm(1, -.25, sd=.5, lower = -2, upper = 2)

#functions to find epsilons, transition probabilities, and m-h ratio

priors <- function(alpha.0, alpha.1, p, q, sigma, phi.1, phi.2, phi.3, phi.4) {
  priors <- dnorm(alpha.1, 1.5, sd=.1) * dnorm(alpha.0, -.36, sd=.1)*dnorm(p, .6, sd=.1)*dnorm(q, .6, sd=.1)*dnorm(sigma, .8, sd=.1)*dnorm(phi.1, 0, sd=.1)*dnorm(phi.2, 0, sd=.1)*dnorm(phi.3, -.25, sd=.1)*dnorm(phi.4, -.25, sd=.1)
  return(priors)
}

priors <- priors(alpha.0.star, alpha.1.star, p.star, q.star, sigma.star, phi.1.star, phi.2.star, phi.3.star, phi.4.star)

nu.norm <- function(eps, probs, sigma = params$sigma, priors=priors) {
  nu <- sum(log(c(eps,probs,priors)))
  return(nu)
}

#ALWAYS REJECTING -- HOW TO GET ACCEPTANCE RATE UP?
#initializing params
params <- data.frame("alpha.1" = 1, "alpha.0" = -.2, 
                     "p" = .8, "q" = .6, "sigma" = .5, "phi.1" = 0,
                     "phi.2" = 0, "phi.3" = 0, "phi.4" = 0)
initial <- sample(c(0,1), 266, replace = TRUE) #mix

mh.norm <- function(N=1000000, burn=10000, initial=initial, params=params, Y.data=Y.data) {
  current <- initial
  Z.current <-  Y.data-params$alpha.1*current-params$alpha.0 #find Z's
  p.current <- sapply(2:266, function(x) trans.prob(s.1=current[x-1], s.2=current[x])) #transitions
  ep.current <- epsilon(Z=Z.current, phi.1=params$phi.1, phi.2=params$phi.2, phi.3=params$phi.3, phi.4=params$phi.4, sigma=params$sigma) #find error terms
  priors.current <- priors(params$alpha.0, params$alpha.1, params$p, params$q, params$sigma, params$phi.1, params$phi.2, params$phi.3, params$phi.4)
  nu.i <- nu.norm(eps=ep.current, probs=p.current, sigma=params$sigma, priors=priors.current)
  
  alpha.1 <- params$alpha.1
  alpha.0 <- params$alpha.0
  p <- params$p
  q <- params$q
  sigma <- params$sigma
  phi.1 <- params$phi.1
  phi.2 <- params$phi.2
  phi.3 <- params$phi.3
  phi.4 <- params$phi.4
  
  
  i <- 1 #counter
  states <- matrix(0, 266, 1) #this will store the states over the course of the simulation in the form of a sum
  params.est <- params #this will store the new parameters
  #metric <- states <- matrix(0, 1, N) #the metric for judging detailed balance will be the % of time in recession (sum(current)/266)
  accept <- 0 #this will count the number of time the proposal was accepted over the course of the simulation
  
  while ( i <= N) {
    #make proposal
    flip <- sample(1:266, 1) #flip 1 state
    proposal <- current
    proposal[flip] <- abs(1-current[flip])
    #proposal <- sample(c(0,1), 266, replace = TRUE)
    
    #propose new parameters based on normal priors
    alpha.1.star <- rtnorm(1, mean=1.5, sd=.005, lower = 1, upper = 2)
    alpha.0.star <-rtnorm(1, mean=-.36, sd=.005, lower = -2, upper = 0)
    p.star <- rtnorm(1, mean=.9, sd = .001, lower = 0, upper = 1)
    q.star <- rtnorm(1, mean=.7, sd = .001, lower = 0, upper = 1)
    sigma.star <- rtnorm(1, .8, sd=.005, lower = .5, upper = 2)
    phi.1.star <- rtnorm(1, .03, sd=.01, lower = 0, upper = 1)
    phi.2.star <- rtnorm(1, -.03, sd=.01, lower = -2, upper = 0)
    phi.3.star <- rtnorm(1, -.25, sd=.005, lower = -2, upper = 0)
    phi.4.star <- rtnorm(1, -.25, sd=.005, lower = -2, upper = 0)
    
    #calculate all the parts of nu
    Z.proposal <-  Y.data - alpha.1.star*proposal-alpha.0.star #find Z's
    p.proposal <- sapply(2:266, function(x) trans.prob(s.1=proposal[x-1], s.2=proposal[x])) #transitions
    ep.proposal <- epsilon(Z=Z.proposal, phi.1=phi.1.star, phi.2=phi.2.star, phi.3=phi.3.star, 
                           phi.4=phi.4.star, sigma=sigma.star) #find error terms
    priors.proposal <- priors(alpha.0.star, alpha.1.star, p.star, q.star, sigma.star, phi.1.star, phi.2.star, phi.3.star, phi.4.star)
    nu.j <- nu.norm(eps=ep.proposal, probs=p.proposal, sigma = sigma.star, priors = priors.proposal)
    
    U <- runif(1)
    
    if (U < min(1,exp(nu.j-nu.i))) {
      current <- proposal
      nu.i <- nu.j
      accept <- accept + 1
      alpha.1 <- alpha.1.star
      alpha.0 <- alpha.0.star
      p <- p.star
      q <- q.star
      sigma <- sigma.star
      phi.1 <- phi.1.star
      phi.2 <- phi.2.star
      phi.3 <- phi.3.star
      phi.4 <- phi.4.star
    }
    
    if(i >= burn) {
    params.est <- data.frame("alpha.1" = alpha.1, "alpha.0" = alpha.0, "p" = p, "q" = q, "sigma" = sigma, "phi.1" = phi.1, "phi.2" = phi.2, "phi.3" = phi.3, "phi.4" = phi.4)
    write.table(params.est, file = "paramest.txt", append = TRUE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
    
  write.table(matrix(current, nrow=1, ncol=266), file = "states.txt", append = TRUE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
    }
      if(i%%1000==0) {
        print(i)
        final <- data.frame(params.est)
      }
    i <- i + 1
    #metric[,i] <- sum(current)/266
  }

  return(accept/N)
}

library(ggplot2)
library(compiler)
cmpfun(priors)
cmpfun(nu.norm)
cmpfun(mh.norm)

test.norm <- mh.norm(N=1000, burn=0, initial=initial, params=params, Y.data=Y.data)
test.norm

final.norm <- mh.norm(N=2000000, burn=50000, initial=initial, params=params, Y.data=Y.data)
final.norm

params.norm <- as.data.frame(paramest[1:50000,])
plot.alpha1 <- qplot(1:50000,y = params.norm$V1, geom="line")
plot.alpha0 <- qplot(1:50000,y = params.norm$V2, geom="line")
plot.p <- qplot(1:50000,y = params.norm$V3, geom="line")
plot.q <- qplot(1:50000,y = params.norm$V4, geom="line")
plot.sigma <- qplot(1:50000,y = params.norm$V5, geom="line")
plot.phi1 <- qplot(1:50000,y = params.norm$V6, geom="line")
plot.phi2 <- qplot(1:50000,y = params.norm$V7, geom="line")
plot.phi3 <- qplot(1:50000,y = params.norm$V8, geom="line")
plot.phi4 <- qplot(1:50000,y = params.norm$V9, geom="line")

plot.alpha1
plot.alpha0
plot.p
plot.q
plot.sigma
plot.phi1
plot.phi2
plot.phi3
plot.phi4
