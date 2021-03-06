rm(list = ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
iter <- as.numeric(repenv)


library(ks, lib.loc = 'packages')
library(mvtnorm, lib.loc = 'packages')
source('vbetelIVaux.R')
source('reparamVB.R')
source('metropolisHastings.R')
library(Rcpp, lib.loc = 'packages')
library(RcppArmadillo, lib.loc = 'packages')
sourceCpp('vbetelMatrixCalculations.cpp')


Nseq <- c(250, 1000, 2000)
N <- max(Nseq)
set.seed(iter)

alpha <- 1
beta <- 0.5
delta <- 0.7

z1 <- rnorm(N, 0.5, 1)
z2 <- rnorm(N, 0.5, 1)
w <- runif(N)

errors <- rmvnorm(N, rep(0, 2), matrix(c(1, 0.7, 0.7, 1), 2))
u <- errors[,2]

eps <- qMix(pnorm(errors[,1]), m1 = 0.5, m2 = -0.5, sd1 = 0.5, sd2 = 1.118, p = 0.5)

x <- z1 + z2 + w + u
y <- alpha + beta * x + delta * w + eps

dataFull <- cbind(y, x, z1, z2, w)

results <- data.frame()
for(N in Nseq){
  data <- dataFull[1:N,]
  draw <- c(alpha = 1, beta = 0.5, delta = 0.5)
  mcmcM1 <- metropolisHastings(start = draw, 
                               iter = 20000, 
                               model = densityM1,
                               data = data,
                               stepsize = 0.1)
 
  mcmcM2 <- metropolisHastings(start = draw, 
                               iter = 20000, 
                               model = densityM2,
                               data = data,
                               stepsize = 0.1)
  
  postDensM1 <- kde(as.matrix(mcmcM1[seq(15001, 20000, 10),1:3]))
  postDensM2 <- kde(as.matrix(mcmcM2[seq(15001, 20000, 10),1:3]))
  
  indexM1 <- which(postDensM1$estimate > 0, arr.ind = TRUE)
  estimM1 <- rep(0, 100)
  for(i in 1:100){
    dim <- sample(nrow(indexM1), 1)
    thetaM1 <- sapply(1:3, function(x) postDensM1$eval.points[[x]][indexM1[dim, x]])
    logJointM1 <- densityM1(data, thetaM1)
    estimM1[i] <- logJointM1 - log(postDensM1$estimate[indexM1[dim,1], indexM1[dim, 2], indexM1[dim, 3]])
  }
  
  indexM2 <- which(postDensM2$estimate > 0, arr.ind = TRUE)
  estimM2 <- rep(0, 100)
  for(i in 1:100){
    dim <- sample(nrow(indexM2), 1)
    thetaM2 <- sapply(1:3, function(x) postDensM2$eval.points[[x]][indexM2[dim, x]])
    logJointM2 <- densityM2(data, thetaM2)
    estimM2[i] <- logJointM2 - log(postDensM2$estimate[indexM2[dim,1], indexM2[dim, 2], indexM2[dim, 3]])
  }
  
  mcmcDiff <- mean(estimM1) - mean(estimM2)
  
  lambda <- c(draw, diag(0.1, 3))
  
  fitM1 <- reparamVB(data = data,
                     lambda = lambda, 
                     model = reparamDerivM1, 
                     alpha = 0.03,
                     S = 10, 
                     dimTheta = 3,
                     maxIter = 1000,
                     threshold = 0.0001 * N,
                     RQMC = FALSE)  
  
  UM1 <- matrix(fitM1$lambda[4:12], 3)
  SigmaM1 <- t(UM1) %*% UM1
   
  fitM2 <- reparamVB(data = data,
                     lambda = lambda, 
                     model = reparamDerivM2, 
                     alpha = 0.03,
                     S = 10, 
                     dimTheta = 3,
                     maxIter = 1000,
                     threshold = 0.0001 * N,
                     RQMC = FALSE)  
  
  UM2 <- matrix(fitM2$lambda[4:12], 3)
  SigmaM2 <- t(UM2) %*% UM2
     
  tM1 <- rmvnorm(100, fitM1$lambda[1:3], SigmaM1)
  tM2 <- rmvnorm(100, fitM2$lambda[1:3], SigmaM2)
  
  vbDiff <- mean(sapply(1:100, function(i){
    logJoint1 <- densityM1(data, tM1[i,]);
    logJoint2 <- densityM2(data, tM2[i,]);
    logPost1 <- dmvnorm(tM1[i, ], fitM1$lambda[1:3], SigmaM1, log = TRUE);
    logPost2 <- dmvnorm(tM2[i, ], fitM2$lambda[1:3], SigmaM2, log = TRUE);
    logJoint1 - logPost1 - (logJoint2 - logPost2)
  }))
  
  results <- rbind(results,
                   data.frame(
                     id = iter,
                     N = N, 
                     mcmcDiff = mcmcDiff,
                     vbDiff = vbDiff))
  
}

write.csv(results, paste0('vbetel/id', iter, '.csv'), row.names = FALSE)


