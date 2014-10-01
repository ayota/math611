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
means <- c(50,90)
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

#comparing to actual genders
cbind(mymeans$cluster,heights.all$Gender)

#graphing distances
plot(x=mymeans$distances, type ="l", ylab = "Value of f(x)", xlab = "Iteration")
