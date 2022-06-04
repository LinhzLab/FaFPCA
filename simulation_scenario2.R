########################
#### load libraries ####
########################

library(MASS)
library(fda)

sim.scenario2 <- function (sample.len, time.len, p, q, K, sigma) {
  
  # sample.len: sample sizes 
  # time.len: number of time points for each individual
  # p: dimension of predictors
  # q: dimension of loadings
  # K: number of eigenfunctions
  # sigma: sd matrix of error 
  
  # generate time points
  
  time <- list()
  for (i in 1 : sample.len) {
    set.seed(1)
    time[[i]] <- sort( runif(time.len[i], min = 0, max = 10) )
  }
  
  
  # generate loading
  
  ideqr.B <- function(X) {
    B.qr <- qr(X)
    B <- sqrt(p) * qr.Q(B.qr)
    for (i in 1 : q) B[,i] <- sign( B[1,i] ) * B[,i]
    return(B)
  }
  
  set.seed(5)
  
  b.var <- matrix(0, p, p)
  for (i in 1 : p) {
    for (j in 1 : p) {
      b.var[i,j] <- 0.5 ^ ( abs(i - j)) 
    }
  }
  b.ran <- mvrnorm(sample.len, rep(0,p), b.var)
  Bn.true <- b.ran %*% t(b.ran) 
  Bn.qr <- qr(Bn.true)
  Bn.true <- t(b.ran) %*% qr.Q(Bn.qr)[,1 : q] 
  B.true <- sqrt(p) * qr.Q( qr(Bn.true) )
  B.true  <- ideqr.B(B.true)
  
  
  # generate eigenfunctions

  eigenfun <- list()
  for (i in 1 : sample.len) {
    eigenfun[[i]] <- array(0,c(time.len[i], K, q))
    
    for (l in 1 : q)  {
      for (k in 1 : (K/2)) {
        eigenfun[[i]][,(2 * k) - 1,l] <- sin(2 * k * pi * time[[i]] / 10) * sqrt(2)
        eigenfun[[i]][,2 * k,l] <- cos(2 * k * pi * time[[i]] / 10) * sqrt(2)
      }
    }
    
    for (k in 1 : (K/2)) {
      eigenfun[[i]][,(2 * k) - 1,2] <- sin(k * pi * time[[i]] / 10) * sqrt(2)
      eigenfun[[i]][,2 * k,2] <- sin((k+2) * pi * time[[i]] / 10) * sqrt(2)
    }
    
  }
  
  for (l in 1 : q) {
    for (k in 1 : K) {
      eigenkl <- c(eigenfun[[1]][,k,l])
      for (i in 2 : sample.len) eigenkl <- c( eigenkl, eigenfun[[i]][,k,l] )
      for (i in 1 : sample.len) eigenfun[[i]][,k,l] <- eigenfun[[i]][,k,l] / norm(eigenkl, type="2") * sqrt( sample.len )
    }
  }
   
  # generate score
  
  idesvd.score <- function(x) {
    score <- matrix(0, sample.len, q * K)
    score.n <- array(NA, c(q, K, sample.len))
    score.n <- x
    for (i in 1 : sample.len) score[i,] <- as.vector( t(score.n[,,i]) )
    score.svd <- svd(score)
    score <- score.svd$u %*% diag(score.svd$d)
    for (j in 1 : ( q * K )) score[,j] <- sign( score[1,j] ) * score[,j]
    for (i in 1 : sample.len) score.n[,,i] <- t(matrix(score[i,], K, q))
    result <- list(score.n = score.n, score = score)
    return(result) 
  }
  
  set.seed(3)
  score.var <- matrix(0, K * q, K * q)
  for (i in 1 : (K * q)) score.var[i,i] <- 1 / i * K * q
  score.true <- mvrnorm(n = sample.len, rep(0, K * q), score.var)
  score.true.ran <- score.true[,1 : (K * q)]
  score.n.true.ran <- array(NA,c(q, K, sample.len))
  for (i in 1 : sample.len) score.n.true.ran[,,i] <- t( matrix(score.true.ran[i,], K, q) )
  score.n.true <- idesvd.score(score.n.true.ran) $ score.n
  score.true <- idesvd.score(score.n.true.ran) $ score
  
  
  # generate functional data
  
  Ht.true <- list()
  for (i in 1 : sample.len) {
    
    Ht.true[[i]] <- matrix(0, q, time.len[i])
    for (j in 1 : q) {
      for (t in 1 : time.len[i]) {
        Ht.true[[i]][j,t] <- t( score.n.true[j,,i] ) %*% eigenfun[[i]][t,,j]
      }
    }
    
  }
  
  X <- list()
  for (i in 1 : sample.len) X[[i]] <- B.true %*% Ht.true[[i]] + matrix( mvrnorm(n = time.len[i], rep(0, p), sigma), 
                                                                        p, time.len[i], byrow = F ) 
  
  
  return( list( time = time, X = X, B = B.true, score = score.true, eigenfun = eigenfun) )
  
}



###################
#### load codes ###
###################
source("FaFPCA.R")


###################################
#### generation of the example ####
###################################

sample.len <- 100 # 500
tt <- 20 # 5, 6, 7, 8, 9, 10, 15, 20, 30 ,40 ,50
time.len <- rep(tt, 100) 
p <- 100 # 500 
q <- 5
K <- 2
sigma <- diag(p)
spline.break <- c(0, 2.5, 5, 7.5, 10)
num.spline <- 5


run.sim <- sim.scenario2(sample.len, time.len, p, q, K, sigma)
X <- run.sim$X
time <- run.sim$time
Bt <- Bt.spline(sample.len, time, spline.break, num.spline)

runcode <- FaFPCA(X, time, sample.len, time.len, p, q, K, num.spline, Bt)

