2. 
(a)

jumps <- function(n, rate) {
  intervals <- rexp(n, rate = rate)
  total <- cumsum(intervals)
  return(total)
}

(b)
lambda <- function(time) {
  max <- 5*10^(-4) + 7.5858*10^(-5)*exp(log(1.09144)*time)
  return(max)
}

death <- function (time) {
  max <- lambda(time)
  samples <- jumps(1000,max)
  proposals <- samples[which(samples <= time)]
  i <- 1
  alive <- TRUE
  time.of.death <- 10000 #this is the null so it won't matter in calculation
  
  while ( alive & i <= length(proposals) ) {
    U <- runif(1)
    if(U <= lambda(proposals[i])/max) {
      time.of.death <- proposals[i]
      alive <- FALSE
    }
    i <- i + 1
  }
  return(time.of.death)
}

death(100)

(d)
##integral
integrand <- function(x) {return((5*10^(-4)+7.5858*10^(-5)*exp(log(1.09144)*x)))}
exponent <- integrate(Vectorize(integrand), 60, 90)
P.90 <- integrate(Vectorize(integrand), 0, 90)
P.60 <- integrate(Vectorize(integrand), 0, 60)

# this is P(T>90 | T >60) = P(T>90) = [1 - P(T<90)] / [1 - P(T<60)]
integrate.death <- exp(-P.90$value)/exp(-P.60$value)

##monte carlo
deaths <- sapply(1:100000, function(x) death(90))
mc.death <- length(deaths[which(deaths > 90)])/length(deaths[which(deaths > 60)])

mc.death
integrate.death


3. Casino cheating

(a) Simulations for X(t) / Y(t)


x.sim <- function(rolls, alpha) {
  #we say F = 1 and C = -1
  x <- 1
  path <- x
  for (i in 1:rolls) {
  U <- runif(1)
  x.new <- ifelse(U <= alpha, x <- -x, x <- x)
  path <- c(path,x.new)
  }
  
  final <- ifelse(path == 1, x <- "F", x <- "C")
  return(final)
}

y.sim <- function(x.state) {
  ifelse(x.state == "F", y <- sample(1:6,1), y <- sample(c(6, rep(1:5,3)),1))
  return(y)
}

chain.sim <- function(rolls, alpha) {
  x<-x.sim(rolls, alpha)
  chain <- sapply(x, function(x) y.sim(x))
  chain.final <- data.frame("X.State" = x, "Y.Roll" = chain)
  return(chain.final)
}

chain.sim(100, .2)

(b)
sim.10000 <- chain.sim(10000, .2)
