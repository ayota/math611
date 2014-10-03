####2A: rv sampler####

sigma <- matrix(c(2, 1, 1, 10), nrow=2,ncol=2,byrow =TRUE)

xsampler <- function(N) {
  #enter the cov matrix sigma and find its eigenvalues/eigenvectors
  sigma <- matrix(c(2, 1, 1, 10), nrow=2,ncol=2,byrow =TRUE)
  sigma.eig <- eigen(sigma)
  #Q is the matrix of eigenvectors of sigma
  Q <- sigma.eig$vectors
  #generate Y's with distritbutions N(0, sigma.sq = eigenvalues)
  set.seed(100)
  y.pairs <- cbind(rnorm(N,mean = 0,sd=sqrt(sigma.eig$values[1])), rnorm(N,mean = 0,sd=sqrt(sigma.eig$values[2])))
  X <- apply(y.pairs, 1, function(x) Q %*% x)
  return(X)
}

####2B: Caculate p-value####
pvalue <- function (N) {
  x <- xsampler(N)
  p.value <- sum(ifelse(((10*x[1,]^2)/19 - (2*x[1,]*x[2,])/19 + (2*x[2,]^2)/19) > 65.789, 1, 0))/N
  return(p.value)
  }

pvalue(100)

###Problem 3###

##3A: K-means##
heights.all <- read.table("~/Desktop/MTH611/heights.txt", header=TRUE, quote="\"")

heights <- heights.all$Height
means <- c(49,90)
norm <- function(x,mean) {
  norm <- sqrt(sum((x-mean)^2))
  return(norm)
}

k.means <- function(testdata, means) {
  data <- cbind(testdata, c(1:length(testdata)))
  test.init <- 1
  test.final <- 0
  mu1 <- means[1]
  mu2 <- means[2]
  means <- c(mu1,mu2)
  distances <- NULL
  
  while(test.init > test.final) {
  test.init <- sum(sapply(data[,1], norm, mean = data[,2])^2)
  distances <- c(distances,test.init)
  
  #assign x's based on mus
  data.1 <- sapply(data[,1], norm, mean = mu1 )
  data.2 <- sapply(data[,1], norm, mean = mu2 )
  data[,2] <- ifelse(data.1^2 < data.2^2, mu1, mu2)
   
  #calculate new means based on classes
  mu1 <- mean(subset(data[,1], data[,2] == mu1 ))
  mu2 <- mean(subset(data[,1], data[,2] == mu2 ))
  means <- rbind(means,c(mu1,mu2))
  
  #assign x's based on mus
  data.1 <- sapply(data[,1], norm, mean = mu1 )
  data.2 <- sapply(data[,1], norm, mean = mu2 )
  data[,2] <- ifelse(data.1^2 < data.2^2, mu1, mu2)
  
  test.final <- sum(sapply(data[,1], norm, mean = data[,2])^2)
  }
  distances <- c(distances,test.final)
  clusters <- ifelse(data[,2] == mu1, 1, 2)
  results <- list(means = c(mu1, mu2), clusters = as.vector(clusters), meanchanges = means, distances = distances)
  return(results)
}

mymeans <- k.means(heights, means)
rmeans <- kmeans(heights, centers=2)

#comparing to R's kmeans generator
cbind(mymeans$cluster,rmeans$cluster)

#comparing to actual data set's means
mu.female <- mean(subset(heights.all$Height, heights.all$Gender==1))
mu.male <- mean(subset(heights.all$Height, heights.all$Gender==2))

print(c(mu.female, mu.male))
mymeans$means

#graphing distances
plot(x=mymeans$distances, type ="l", ylab = "Value of f(x)", xlab = "Iteration")

#comparing to actual genders
diffs <- ifelse(mymeans$cluster == heights.all$Gender, 0, 1)
sum(diffs)
#classifies 87/100 correctly

heights <- hope$Height

#using hint in reading, I made the starting variance the variance of the entire sample, and the starting mu's two random values from the sample
params <- c(65, 74, 16, 17, .5) #mu1, mu2, sig1, sig2, p

EM.est <- function(x, params) {
  x <-heights
  mu1 <- params[1]
  mu2 <- params[2]
  sig1 <- params[3]
  sig2 <- params[4]
  p <- params[5]
  log <- loglike(x, mu1, mu2, sig1, sig2, p)
  path <- NULL

  #while (abs(path[i,1] - path[i-1,1]) > .01 | abs(path[i,2] - path[i-1,2]) > .01 | abs(path[i,3] - path[i-1,3]) > .01 | abs(path[i,4] - path[i-1,4]) > .01 ) 
  
  #the while loop runs until the path starts to repeat
  while (!(mu1 %in% path) | !(mu2 %in% path) | !(sig1 %in% path) | !(sig2 %in% path)  ) {
    path <- rbind(path,c(mu1, mu2, sig1, sig2, p, log))
    
  #expectation step
  r1 <- p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) / (p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)))
  r2 <- (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)) / (p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1)))
  
  #maximization step
  mu1 <- sum((1-r1) * x)/sum((1-r1))
  mu2 <- sum(r1 * x)/sum(r1)
  
  sig1 <- sum((1-r1) * (x-mu1)^2)/sum(1-r1)
  sig2 <- sum(r1 * (x-mu2)^2)/sum(r1)
  
  p <- 1/100 * sum(r1)
  p2 <- 1/100 * sum(r2)
    
  log <- loglike(x, mu1, mu2, sig1, sig2, p)
  
  }

  return(path)
}

path <- EM.est(heights, params)

#graph log likelihood
loglike <- function(x,mu.one,mu.two,sigma.one,sigma.two,p) {
  fxn <- sum( log(p*(1/sqrt(2*pi*sigma.one))*exp(-(x-mu.one)^2/(2*sigma.one)) + (1-p)*(1/sqrt(2*pi*sigma.two))*exp(-(x-mu.two)^2/(2*sigma.two))) )
  return(fxn)
}

plot(path[2:183,6], ylim=c(-281.6,-281.46))

#part c ii
mu1 <- path[183,1]
mu2 <- path[183,2]
sig1 <- path[183,3]
sig2 <- path[183,4]
p <- path[183,5]

r1 <- sum(p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1))) / sum((p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1))))
r2 <- sum((1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1))) / sum((p*(sqrt(2*pi*sig1)^(-1))*exp((-(x-mu1)^2)*(2*sig1)^(-1)) + (1-p)*(sqrt(2*pi*sig2)^(-1))*exp((-(x-mu2)^2)*(2*sig2)^(-1))))


