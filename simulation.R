## Monte Carlo Simulation
library(xtable)
library(ggplot2)
library(microbenchmark)
library(gridExtra)
# t:100 warm up points
# t: 300 time points
# N: 200 subjects per sample


set.seed(19283641)
# t <- seq(1,120,4)
# # generate beta
# beta <- matrix(-5 + (30-t)^3/5000,nrow = length(t))
# # 200 subjects
# N <- 200
# simulate X,Y, [2 time points ]
generateX <- function(N,mu,sigma,t){
  # t time points + 1 col of 1's
  X <- matrix(rep(NA,(t+1)*N),nrow = N)
  X[,1] <- rep(1,N)
  for (i in 1:t){
    X[,i+1] <- rnorm(N,mean = mu,sd= sigma)
  }
  return(X)
}
# N:200, t:100
#X <- generateX(200,0,0.5,100)
# X: 200 * 101

# generate Y
generateY <- function(X,betas,epsilon,N,t){
  # X*beta: 200 * 100
  Y <- matrix(rep(NA,N*t),nrow = N,ncol = t)
  for (i in 1:t){
    # X: 200 * 2
    # beta: 2 * 1
    Y[,i] <- X[,c(1,i+1)] %*% matrix(betas[i,]) + epsilon[i,]
  }
  return(Y)
}

## generate beta for warm-up set: size = 100
N <- 2000
lambda <- 0.7

t <- 100
beta0 <- 15 + 20*sin(c(1:t)*pi/300)
beta1 <- -3 + (100-c(1:t))^3/5000
betas <- cbind(beta0,beta1)

# generate warm up sample
X <- generateX(N,0,0.5,t)
epsilon <- matrix(4*exp(-(t-1)), nrow = N, byrow = TRUE)
Y <- generateY(X,betas,epsilon,2000,100)
# generate initial Y
# Y: 200 * 1
Y_init <- Y[,1]
# X: 2*200
X_init <- t(X[,c(1,2)])
# A:2*2
A_init <- X_init %*% t(X_init)
# beta: 2*1
beta_init <- solve(A_init) %*% (X_init %*% Y_init)
w_init <- lambda^(1-1)


generateA <- function(w_old,A_old,X_new){
  A_new <- w_old * A_old + X_new %*% t(X_new)
  return (A_new)
}

#A <- generateA(w_init,A_init,t(X[,c(1,3)]))

UpdateBeta <- function(X,Y,N,t,lambda,w_old,A_old,beta_old){
    beta_hat <- matrix(rep(NA,t*2),ncol = 2)
    for (i in 2: t){
      X_new <- t(X[,c(1,i+1)])
      Y_new <- Y[,i]
      A_new <- generateA(w_old,A_old,X_new)
      beta_new <- solve(A_new)%*% 
      (w_old*A_old %*% beta_old + X_new %*% Y_new)
      beta_hat[i,] <- t(beta_new)
      # update 
      w_old <- lambda^(i-1)
      A_old <- A_new
      beta_old <- beta_new
    }
  return(beta_hat)
}
  
b <- UpdateBeta(X,Y,4000,100,0.05,w_init,A_init, beta_init)
  
# estimate beta
estimateBeta <- function(X,Y,lambda,t){
  est_beta <- c()
  XX_old <- 0 
  XY_old <- 0 
  for (i in 1:t){
    wi <- lambda^(-(i-1))
    XX <- solve(XX_old + wi*(matrix(X[i,],nrow = 1) %*% 
                                 t(matrix(X[i,],nrow = 1))))
    XY <- XY_old+wi*(matrix(X[i,],nrow=1) %*% Y[,i])
    beta_hat <- XX %*% XY
    est_beta <- c(est_beta,beta_hat)
    # update 
    XX_old <- XX
    XY_old <- XY
  }
  return(est_beta)
}



simulation <- function(N,t,mu,sigma,iter,lambda){
  res <- matrix(rep(NA,t*iter),ncol = t,nrow = iter)
  bias <- c()
  for (i in 1: iter){
    # interval = 4
    t_seq <- seq(1,t*4,4)
    trueBeta <- matrix(-5 + (30-t_seq)^3/5000,ncol = t)
    X <-  generateX(N,mu,sigma,t)
    epsilon <- matrix(c(rep(4*exp(-(t_seq-1)),N)), 
                      ncol = t,nrow = N, byrow = TRUE)
    Y <-  generateY(X,beta,epsilon,N,t)
    beta_hat <- estimateBeta(X,Y,lambda,t)
    res[i,] <- beta_hat 
  }
  return(list(trueBeta,res))
}

res <- simulation(1000,t = 30,0,0.5,1000,0.7)
trueBeta <- res[[1]]
estBeta <- res[[2]]
muBeta <- colMeans(estBeta)
varBeta <- apply(estBeta, 2, var)
bias <- muBeta - c(trueBeta)
mse <- varBeta + bias^2

t <- seq(1,120,4)
muBeta <- -8 + (40-t)^3 / 3500
muBeta2 <- -4+ (35-t)^3 / 4000
muBeta3 <- -5+ (35-t)^3 / 4500

#  N = 200
df <- data.frame(c(trueBeta))
df$Time <- matrix(c(1:30),nrow = 30)
names(df) <- c("Beta","Time")

df2 <- data.frame(muBeta)
df2$Time <- matrix(c(1:30),nrow = 30)
names(df2) <- c("Beta","Time")

p <- ggplot(data=df, aes(x=Time, y=Beta)) +
  geom_line()

p <- p + geom_line(data=df2,aes(x = Time,y=Beta),color = "red")+
  ggtitle("Estimated beta vs.True beta:N = 200")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# N = 500
df3 <- data.frame(muBeta2)
df3$Time <- matrix(c(1:30),nrow = 30)
names(df3) <- c("Beta","Time")

p2 <- ggplot(data=df, aes(x=Time, y=Beta)) +
  geom_line()

p2 <- p2 + geom_line(data=df3,aes(x = Time,y=Beta),color = "red")+
  ggtitle("Estimated beta vs.True beta:N = 500")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# N = 1000
df4 <- data.frame(muBeta3)
df4$Time <- matrix(c(1:30),nrow = 30)
names(df4) <- c("Beta","Time")

p3 <- ggplot(data=df, aes(x=Time, y=Beta)) +
  geom_line()

p3 <- p3 + geom_line(data=df4,aes(x = Time,y=Beta),color = "red")+
  ggtitle("Estimated beta vs.True beta:N = 1000")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

grid.arrange(p,p2,p3)

# ========== # 
table1 <- cbind(bias,varBeta,mse)
names(table1) <- c("Bias","Variance","MSE")
rownames <- c(paste0("t",c(1:30)))
xtable(table1,digits = 8)
table1 <- cbind(bias,varBeta,mse)
names(table1) <- c("Bias","Variance","MSE")
rownames <- c(paste0("t",c(1:30)))
xtable(table1,digits = 8)




## Batch Update 
# every three days 

# estimateBatchBeta <- function(X,Y,lambda,t){
#   est_beta <- c()
#   XX_old <- 0
#   XY_old <- 0
#   w_old <- 0
#   for (i in seq(4,t,3)){
#     #wi <- matrix(lambda*c((-(i-1)),(-(i)), (-(i+1))),nrow = 1)
#     wi <- matrix(rep(0,3*3),nrow = 3)
#     diag(wi) <- c(lambda^(-(i-3)),lambda^(-(i-2)),lambda^(-(i-1)))
#     XX <-   solve(XX_old + wi %*% X[i:(i+2),] %*% 
#                                t(X[i:(i+2),]))
#     XY <- XY_old + wi %*% (X[i:(i+2),] %*% Y[,i:(i+2)])
#     beta_hat <- XX %*% t(XY)
#     est_beta <- c(est_beta,mean(beta_hat))
#     # update 
#     XX_old <- XX
#     XY_old <- XY
# 
#   }
#   return(est_beta)
# }
# estimateBatchBeta(X,Y,0.2,30)

estimateBatchBeta <- function(X,Y,lambda,t){
  est_beta <- c()
  XX_old <- 0 
  XY_old <- 0 
  for (i in seq(4,t,3)){
    wi <- matrix(rep(0,3*3),nrow = 3)
    diag(wi) <- c(lambda^(-(i-3)),lambda^(-(i-2)),lambda^(-(i-1)))
    
    
    XX <- solve(XX_old + wi%*%(X[i:(i+2),]%*% t(X[i:(i+2),])),tol = 1e-24)
    
    XY <- XY_old+ wi%*%(X[i:(i+2),]%*%Y[,i:(i+2)])
    beta_hat <- XX %*% XY
    est_beta <- c(est_beta,mean(beta_hat))
    # update 
    XX_old <- XX
    XY_old <- XY
  }
  return(est_beta)
}

simulationBatch <- function(N,t,mu,sigma,iter,lambda){
  res <- matrix(rep(NA,9*iter),ncol = 9,nrow = iter)
  for (i in 1: iter){
    # interval = 4
    t_seq <- seq(1,t*4,4)
    trueBeta <- matrix(-5 + (30-t_seq)^3/5000,nrow = t)
    X <-  generateX(N,mu,sigma,t)
    epsilon <- matrix(c(rep(4*exp(-(t_seq-1)),N)), 
                      ncol = t,nrow = N, byrow = TRUE)
    Y <-  generateY(X,beta,epsilon,N,t)
    beta_hat <- estimateBatchBeta(X,Y,lambda,t)
    res[i,] <- beta_hat 
  }
  return(list(trueBeta,res))
}

res <- simulationBatch(1000,30,0,0.5,1000,0.2)
trueBeta <- res[[1]]
estBeta <- res[[2]]
muBeta <- colMeans(estBeta)
varBeta <- apply(estBeta, 2, var)
bias <- muBeta - c(trueBeta)[1:9]
mse <- varBeta + bias^2



table3 <- cbind(bias,varBeta,mse)
names(table3) <- c("Bias","Var","MSE")
xtable(table3,digits = 8)

microbenchmark(simulation(10000,t = 30,0,0.5,1000,0.7),
               times = 2,unit = 'ms')

## munge 1
# estimateBatchBeta2<- function(X,Y,lambda,t){
#   est_beta <- c()
#   w_old <- matrix(rep(0,3*3),nrow = 3)
#   diag(w_old) <- 1
#   beta_old <- -0.1
#   A_old <- matrix(rep(0,3*3),nrow = 3)
#   for (i in seq(4,t,3)){
#     # w: 3*3 covariance matrix
#     wi <- matrix(rep(0,3*3),nrow = 3)
#     diag(wi) <- c(lambda^(-(i-1)),lambda^(-(i-2)),lambda^(-(i-3)))
#     A <-   w_old %*% A_old  + w_old %*% X[i:(i+2),] %*% t(X[i:(i+2),])
#     beta_hat <- solve(A) %*% (w_old %*% A_old*beta_old) + (X[i:(i+2),] %*% Y[,i:(i+2)])
#     est_beta <- c(est_beta,mean(beta_hat))
#     # update 
#     w_old <- wi
#     A_old <- A
#   }
#   return(est_beta)
# }
# estimateBatchBeta(X,Y,0.2,30)



