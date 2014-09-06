#input variables
#lamda = rate of T (Interarrival time)
#mu = rate of S (Service time)
#S = Service time for ith person
#T = Interarrival time for ith person
#W = Wait time for ith person
#customer = number of customer of interest
#N = number of samples
#wait = wait time is > than (i.e. c)


genQueue <- function (mu,lamda, N) {

#create matrix for S/T RVs
arrivals <- matrix(0, 100, 2, byrow =TRUE)
services <- matrix(0, 99, 2, byrow =TRUE)
#mu <- .01
#lamda <- .09

#fill out values for service time
services <- rexp(100,rate=mu)

#fill out values for interarrival time
arrivals <- rexp(100,rate=lamda)

#Loop over service/interrival times to calculate wait times

#generate matrix to store wait times
waiting <- matrix(0, 100, 1, byrow =TRUE)

#skipping W1, we calculate the wait times for subsequent customers

for (i in 2:100){
  waiting[i] <- (waiting[i-1] + services[i-1] - arrivals[i])
}

#Repeat function for desired number of times and store all results in n x 100 matrix
#create the matrix to hold the wait times
waitimes <- NULL

#fill the matrix with values
for (i in 1:N) {
  sample <- generateTimes(mu,lamda)
  waitimes[i] <- sample[100]
}

return(waitimes)
}

#this function takes 3.242 seconds with 10,000 samples, .32 seconds with 1000
system.time(sample <- generatewaits(10000,1,.9))


#Test list of waiting times for desired customer to see whether it is over/under desired waiting period
#takes params wait (length of wait testing) / customer (whose wait time calculating) / N (no. of simulations)
#also takes a 100 by n matrix of values to test

sample <-genQueue(1000,1,.9)

testResults <- function (wait, customer, N, sample) {

for (i in 1:N) {
  if (sample[customer,i] > wait){
    results[i] <- 1
  } else {
    results[i] <- 0
  }
}

return(results)
}


results2 <- NULL
y <- matrix(1, 100, 1)
results2 <- y[sample > 0]




#This function calculates the mean, variance, error, and confidence interval (95%)
#it takes a matrix of result values
meanvar <- function (results, N) {
  
  mean <- sum(results)/N
  variances <- matrix(0,N,1)
  for (i in 1:N) {
    variances[i] <- (results[i] - mean)^2
  }
  var <- (1/N)*sum(variances)
  
  sd <- sqrt(var)
  error <- qnorm(.95)*sd/sqrt(N)
  left <- mean - error
  right <- mean + error
  
  totals <- cbind(mean, var, error, left, right)
  
  return(totals)
}

#histogram of the distribution of 10,000 simulations for W100
sample <-genQueue(10000,1,.9)
w100 <- sample

w100hist <- hist(w100)

#QQplot
w100qq <- qqnorm(w100)
w100qq <- qqline(w100,col = 2,lwd=2,lty=2)

#histogram of m
w100s <- NULL

for(i in 1:100) {
  sample <- genQueue(1000, 1, .9)
  w100s <- sum(sample)/1000
}
