## gamma simulation 
set.seed(1927351)
# set parameters here
N <- 200
t <- 100
t_seq <- 1:t

linkFunc <- function(v){
  res <- (1+exp(-v))^(-1)
  return(res)
}

# compute the first order derivative of link function 
linkDev <- function(v,dv){
  res <- (dv * exp(-v)) / (1+exp(-v))^2
  return(res)
}

# generate X
generateX <- function(N,mu,sigma,t){
  # t time points + 1 col of 1's
  X <- matrix(rep(NA,(t+1)*N),nrow = N)
  X[,1] <- rep(1,N)
  for (i in 1:t){
    X[,i+1] <- rnorm(N,mean = mu,sd= sigma)
  }
  return(X)
}

# generate gamma
gamma0 <- -15 + 20*sin(t_seq*pi/60)
gamma1 <- 4 - ((t_seq-20)/10)^2
gammas <- cbind(beta0,beta1)

# generate epsilon
epsilon <- 4*exp(-(t_seq-1))

# generate score
generateScore <- function(X,gammas,epsilon,N,t){
  # X*beta: 200 * 100
  score <- matrix(rep(NA,N*t),nrow = N,ncol = t)
  for (i in 1:t){
    score[,i] <- X[,c(1,i+1)] %*% matrix(gammas[i,]) + epsilon[i]
  }
  score <- linkFunc(score)
  return(score)
}


score <- generateScore(X,gammas,epsilon,200,100)



# generate D
# E[ D= 1 | X]
#D <- matrix(rep(1,N),nrow = N)

# simulate data
# generate warm up sample
X <- generateX(N,0,0.5,100)
epsilon <- 4*exp(-(t_seq-1))
#Y <- generateY(X,gammas,epsilon,200,100)


# initial estimate
X_init <- t(X[,1:2])
gamma_init <- X_init %*% (D - 1/2)

# update
X_new <-  t(X[,c(1,3)])

