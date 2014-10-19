#2 a)
image(M)
M <- matrix(0,100,100)

#step = number of iterations
#

sampler <- function(step) {
  #initialize intial graph
  M <- matrix(0,100,100)
  
  sums <- 0
  
  for(i in 1:step) {
    i <- trunc(runif(1,1,101)) #choose an i value at random
    j <- trunc(runif(1,1,101)) #choose a j value at random
    
    move <- trunc(runif(1,0,2)) #determine 1 -- move or 0 -- don't move
    
    #If the move is 1, then we make sure we can flip our new point. If 0, then we pick another point.
    if(move) {
      #Now we will check to see if we can flip our new point
      flip <- TRUE
      if(i+1 <= 100) {
        if(M[i+1,j]==1) {flip <- FALSE} #Check right
      }
      if (i >=2 & flip) {
        if(M[i-1,j]==1) {flip <- FALSE} #Check left
      }
      if(j+1<=100 & flip) {
        if(M[i,j+1]==1) {flip<-FALSE} #Check upper
      }
      if(j>=2 & flip) {
        if(M[i,j-1] == 1) {flip<-FALSE} #Check lower
      }
      if(flip) {M[i,j]<-1}
    } else {M[i,j]<-0}
  #sums <- rbind(sums,sum(M)/10000)
  }
  return(sum(M))
}

M2 <- sampler(100000)

cmpfun(sampler)
M3 <- sapply(1:1000, function(x) sampler(30000)) 
hist(M3)
image(M2)

plot(M2,type="l")
