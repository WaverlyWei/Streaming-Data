## Monte Carlo Simulation
# 30 time points
t <- seq(1,120,4)
# generate beta
beta <- matrix(-5 + (30-t)^3/5000,nrow = length(t))
# 200 subjects
N <- 200
# simulate X,Y, [2 time points ]
generateX <- function(N,mu,sigma,t){
  X <- matrix(rep(NA,t*N),ncol = 200)
  for (i in 1:t){
    X[i,] <- rnorm(N,mean = mu,sd= sigma)
  }
  return(X)
}
X <- generateX(N,0,0.5,30)
# X: 200 * 30

# matrix epsilon 
epsilon <- matrix(c(rep(4*exp(-(t-1)),N)), 
                    ncol = length(t),nrow = 200, byrow = TRUE)
# generate Y
generateY <- function(X,beta,epsilon,N,t){
  Xbeta <- matrix(rep(NA,N*t),nrow = t,ncol = N)
  for (i in 1:t){
    Xbeta[i,] <- X[i,] * beta[i,]
  }
  Y <- t(Xbeta) + epsilon
  return(Y)
}

Y <- generateY(X,beta,epsilon,200,30)

# estimate beta
estimateBeta <- function(X,Y,lambda,t){
  est_beta <- c()
  XX_old <- 0 
  XY_old <- 0 
  for (i in 1:t){
    wi <- lambda*exp(-(i-1))
    XX <- solve(XX_old + wi*(matrix(X[i,],nrow = 1) %*% 
                                 t(matrix(X[i,],nrow = 1))))
    XY <- XY_old+wi*(matrix(X[i,],nrow=1) %*% Y[,i])
    beta_hat <- XX %*% XY
    est_beta <- c(est_beta,beta_hat)
  }
  return(est_beta)
}

simulation <- function(N,t,mu,sigma,iter,lambda){
  res <- matrix(rep(NA,t*iter),ncol = t,nrow = iter)
  bias <- c()
  for (i in 1: iter){
    # interval = 4
    t_seq <- seq(1,t*4,4)
    beta <- matrix(-5 + (30-t_seq)^3/5000,nrow = t)
    X <-  generateX(N,mu,sigma,t)
    epsilon <- matrix(c(rep(4*exp(-(t_seq-1)),N)), 
                      ncol = t,nrow = N, byrow = TRUE)
    Y <-  generateY(X,beta,epsilon,N,t)
    beta_hat <- estimateBeta(X,Y,lambda,t)
    res[i,] <- beta_hat 
  }
  return(res)
}

res <- simulation(200,t = 30,0,0.5,1000,0.7)
var(res[,1])

microbenchmark::microbenchmark(simulation(200,t = 30,0,0.5,1000,0.7),
                               times = 2)
