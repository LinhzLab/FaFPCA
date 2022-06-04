########################
#### load libraries ####
########################

library(MASS)
library(fda)


sim.scenario1 <- function (sample.len, time.len, p, K, sigma) {
  
  # sample.len: sample sizes 
  # time.len: number of time points for each individual
  # p: dimension of predictors
  # K: number of eigenfunctions
  # sigma: sd matrix of error 
  
  
  # generate time points
  
  time <- list()
  for (i in 1 : sample.len) {
    set.seed(1)
    time[[i]] <- sort( runif(time.len[i], min = 0, max = 10) )
  }
  

  # generate eigenfunctions
  
  eigenfun <- list()
  for (i in 1 : sample.len) {
    eigenfun[[i]] <- matrix(0, K, time.len[i])
    for (k in 1 : K) eigenfun[[i]][k,] <-  sin( (2 * k - 1) * pi * time[[i]] / 10 ) * sqrt(2)
  }

  
  # generate FPCA score
  
  score.true <- array(0, c(p, K, sample.len))
  evar <- array(0, c(p, p, K))
  for (k in 1 : K) {
    
    for (i in 1 : p) {
      for (j in 1 : p) {
        evar[i,j,k] <- 0.5 ^ ( abs (i - j) )
      }
    }
    diag( evar[,,k] ) <- diag( evar[,,k] ) * (4 - k)
    
  }

  for (k in 1 : K) score.true[,K,] <- t( mvrnorm(sample.len, mu = rep(0,p) , evar[,,k]) )
  
  
  # generate functional data

  X <- list()
  for (i in 1 : sample.len) X[[i]] <- score.true[,,i] %*% eigenfun[[i]] + matrix( mvrnorm(n = time.len[i],
                                                                                rep(0,p),sigma), p, time.len[i], byrow = F ) 
  
  
  return( list( time = time, X = X) )
  
}



###################
#### load codes ###
###################
source("FaFPCA.R")


###################################
#### generation of the example ####
###################################

sample.len <- 100
time.len <- rep(20, 100)
p <- 100
K.setting1 <- 2
K.setting2 <- 10
sigma.setting1 <- diag(p)
set.seed(6)
sigma.setting2 <- diag(p)
for (i in 1 : p) {
  for (j in 1 : p) sigma.setting2[i,j] <- 0.3
  sigma.setting2[i,i] <- 1 
}
q <- 60 # 80
K <- 4 # 2
spline.break <- c(0, 2.5, 5, 7.5, 10)
num.spline <- 5


run.sim <- sim.scenario1(sample.len, time.len, p, K.setting1, sigma.setting1)
X <- run.sim$X
time <- run.sim$time

Bt <- Bt.spline(sample.len, time, spline.break, num.spline)
BT <- as.matrix( Bt[[1]] )
for (i in 2 : sample.len)  BT <- rbind(BT, as.matrix( Bt[[i]] ))
Bt.qr <- qr(BT)
BT <- qr.Q(Bt.qr)
Bt[[1]] <- BT[1 : time.len[1], ]
for (i in 2 : sample.len) Bt[[i]] <- BT[ ( sum( time.len[1 : (i - 1)] ) + 1 ) :  sum( time.len[1 : i] ), ]

runcode <- FaFPCA(X, time, sample.len, time.len, p, q, K, num.spline, Bt)


