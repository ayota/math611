###PART C###
#Input Variables##
##lambda-> rate for interarrival times
##mu-> rate for service times
##N-> Desired number of simulations
##customer-> desired number of customers in the line

queue <- function(lambda, mu, N, customer) {
  #list to hold results
  #first column of each matrix is time since first person arrived, second is the wait time of each     customer in line
  queueTimes <- list()
  
  #loop that conducts desired number of simulations, aka N
  for (i in 1:N) {
    
    #matricies to hold wait times to be generated
    #names <- c("Current Time", "Waiting Time", "Departure Time", "No. in Line", "Interarrival", "Service")
    waits <- matrix(0,customer,6,)
    #dimnames(waits) <- list(names,seq(1,customer,1))
    
    #generating exponential RVs to model interarrival times (arrivals) and service times (services)
    arrivals <-rexp(customer,rate=lambda)
    services <- rexp(customer,rate=mu)
    
    #loop to generate desired number of wait times
    for (j in 2:customer) {
      #calculating and recording cumulative times
      waits[j,1] <- sum(arrivals[1:j])
      
      #calculating and recording wait times
      #this formula is derived by translating the base formula, D - A, into terms of W, S, and T
      #W1 = 0
      #W2 = D1 - A2 = W1 + S1 + T1 - T1 -T2 = 0 + S1 - T2
      #W3 = D2 - A3 = W2 + S2 + T1 + T2 - T1 - T2 - T3 = W2 + S2 - T3
      #W4 = D3 - A4 = W3 + S3 + T1 + T2 + T3 -T1 -T2 -T3 - T4 = W3 + S3 - T4
      #We see all terms cancel except for W[i-1], S[i-1], and T[i]
      waits[j,2] <- max(0,(waits[j-1,] + services[j-1] - arrivals[j]))
      
      #calculating the number of customers at the given time
      #departure times
      departures <- services[j-1] + sum(arrivals[1:j-1])
      waits[j-1,3] <- departures
      
      #number of people in line
      people <- NULL
      
      for (k in 1:j) {
        if (departures > (sum(arrivals[1:j]))) {
          people[k] <- 1 } else if (departures <= (sum(arrivals[1:j]))) {people[k] <- 0}
      }
      waits[j,4] <- sum(people)
      
      waits[j,5] <- arrivals[j]
      waits[j,6] <- services[j]
      
    }
    queueTimes[[i]] <- waits[1:customer,]
    }
  
  return(queueTimes)
}


###PART D###
###Rewrite of Part C function to be more efficient###
queue <- function(lambda, mu, N, customer) {
  #matrix to hold times for desired nth customer
  queueTimes <- matrix(0, N)
  
  #loop that conducts desired number of simulations, aka N
  for (i in 1:N) {
    
    #matricies to hold wait times to be generated
    waits <- matrix(0,customer,1)
    
    #generating exponential RVs to model interarrival times (arrivals) and service times (services)
    arrivals <-rexp(customer,rate=lambda)
    services <- rexp(customer,rate=mu)
    
    #loop to generate desired number of wait times
    for (j in 2:customer) {
      waits[j,] <- max(0,(waits[j-1,] + services[j-1] - arrivals[j]))
    }
    
    queueTimes[i] <- waits[customer,]
  }
  
  return(queueTimes)
}
###MEAN, VARIANCE, AND CI FUNCTION###
meanvar <- function (results, N) {
  
  mean <- sum(results)/N
  variances <- matrix(0,N,1)
  for (i in 1:N) {
    variances[i] <- (results[i] - mean)^2
  }
  var <- (1/N)*sum(variances)
  
  sd <- sqrt(var)
  error <- 1.96*sd/sqrt(N)
  left <- mean - error
  right <- mean + error
  
  totals <- cbind(mean, var, error, left, right)
  
  return(totals)
}

##calculate P(W2 = 0) for N = 100##
sample <- queue(.9,1,100,2)
results <- NULL
for (i in 1:100) {
  if (sample[i,] == 0) {
    results[i] = 1
  } else if (sample[i,] > 0) {
    results[i] = 0
  }
}

summary <- meanvar(results,100)

##calculate P(W2 = 0) for N = 1000##
sample2 <- queue(.9,1,1000,2)
results2 <- NULL
for (i in 1:1000) {
  if (sample2[i,] == 0) {
    results2[i] = 1
  } else if (sample2[i,] > 0) {
    results2[i] = 0
  }
}

summary2 <- meanvar(results2,1000)

##calculate P(W2 = 0) for N = 10000##
sample3 <- queue(.9,1,10000,2)
results3 <- NULL
for (i in 1:10000) {
  if (sample3[i,] == 0) {
    results3[i] = 1
  } else if (sample3[i,] > 0) {
    results3[i] = 0
  }
}

summary3 <- meanvar(results3,10000)

##calculate P(W100 = 0) for N = 100##
sample4 <- queue(.9,1,100,100)
results4 <- NULL
for (i in 1:100) {
  if (sample4[i,] <= 0) {
    results4[i] = 1
  } else if (sample4[i,] > 0) {
    results4[i] = 0
  }
}

summary4 <- meanvar(results4,100)

##calculate P(W100 = 0) for N = 1000##
sample5 <- queue(.9,1,1000,100)
results5 <- NULL
for (i in 1:1000) {
  if (sample5[i,] <= 0) {
    results5[i] = 1
  } else if (sample5[i,] > 0) {
    results5[i] = 0
  }
}

summary5 <- meanvar(results5,1000)

##calculate P(W100 = 0) for N = 10000##
sample6 <- queue(.9,1,10000,100)
results6 <- NULL
for (i in 1:10000) {
  if (sample6[i,] <= 0) {
    results6[i] = 1
  } else if (sample6[i,] > 0) {
    results6[i] = 0
  }
}

summary6 <- meanvar(results6,10000)

###PART E###
w100s <- queue(.9,1,10000,100)

#histogram of 10,000 simulations of W100#
w100hist <-hist(w100s)

#qqplot of 10,000 simulations of W100#
w100qq <- qqnorm(w100s)
w100qq <- qqline(w100s,col = 2,lwd=2,lty=2)

##PART F##
m <- NULL
for (i in 1:10000) {
  w100 <- queue(.9,1,1000,100)
  m[i] <- sum(w100)/1000
}

mhist <- hist(m)
mqq <- qqnorm(m)
mqq <- qqline(m,col = 2,lwd=2,lty=2)
