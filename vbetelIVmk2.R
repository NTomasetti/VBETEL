rm(list = ls())
repenv <- Sys.getenv("SLURM_ARRAY_TASK_ID")
iter <- as.numeric(repenv)

library(mvtnorm, lib.loc = 'packages')
source('vbetelIVaux.R')
source('reparamVB.R')
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
  draw <- c(alpha = 1, beta = 0.5, delta = 0.5, nu = 1)
  mode <- nlm(negLogETEL, draw, data = data, hessian = TRUE)
  covar <- solve(mode$hessian)
  
  mcmcFull <- betelMCMC(start = draw, 
                        iter = 20000, 
                        burnIn = 10000,
                        thin = 10,
                        model = densityNested,
                        mode = mode$estimate,
                        covar = covar,
                        data = data)
  
  meanFull <- colMeans(mcmcFull[,1:4])
  
  
  draw <- c(alpha = 1, beta = 0.5, delta = 0.5)
  modeRes <- nlm(negLogETELres, draw, data = data, hessian = TRUE)
  covarRes <- solve(modeRes$hessian)
  mcmcRestricted <- betelMCMC(start = draw, 
                              iter = 20000, 
                              burnIn = 10000,
                              thin = 10,
                              model = densityRestricted,
                              mode = modeRes$estimate,
                              covar = covarRes,
                              data = data)
  
  meanRes <- colMeans(mcmcRestricted[,1:3])

  
  logjoint <- densityNested(data, meanFull)
  propFull <- proposalDens(meanFull, mode$estimate, solve(covar), 5)
  
  index <- sample(1:nrow(mcmcFull), 100)
  E1 <- 0
  E2 <- 0
  for(i in index){
    theta <- as.numeric(mcmcFull[i, 1:4])
    tDens <- densityNested(data, theta)
    tProp <- proposalDens(theta, mode$estimate, solve(covar), 5)
    alpha <- exp(logjoint - tDens + tProp - propFull)
    E1 <- E1 + 1/100 * alpha * exp(propFull)
    
    qDraw <- mode$estimate + t(chol(covar)) %*% rt(4, 5)
    qDens <- densityNested(data, qDraw)
    qProp <- proposalDens(qDraw, mode$estimate, solve(covar), 5)
    alphaQ <- exp(qDens - logjoint + propFull - qProp)
    E2 <- E2 + 1/100 * alphaQ
  }
  logMargFull <- logjoint - log(E1/E2)
  
  
  logjointRes <- densityRestricted(data, meanRes)
  propRes <- proposalDens(meanRes, modeRes$estimate, solve(covarRes), 5)
  
  index <- sample(1:nrow(mcmcRestricted), 100)
  E1 <- 0
  E2 <- 0
  for(i in index){
    theta <- as.numeric(mcmcRestricted[i, 1:3])
    tDens <- densityRestricted(data, theta)
    tProp <- proposalDens(theta, modeRes$estimate, solve(covarRes), 5)
    alpha <- exp(logjointRes - tDens + tProp - propRes)
    E1 <- E1 + 1/100 * alpha * exp(propRes)
    
    qDraw <- modeRes$estimate + t(chol(covarRes)) %*% rt(3, 5)
    qDens <- densityRestricted(data, qDraw)
    qProp <- proposalDens(qDraw, modeRes$estimate, solve(covarRes), 5)
    alphaQ <- exp(qDens - logjointRes + propRes - qProp)
    E2 <- E2 + 1/100 * alphaQ
  }
  logMargRes <- logjointRes - log(E1/E2)
  
  mcmcDiff <- logMargRes - logMargFull
  
  
  lambda <- c(1, 0.5, 0.5, 1, diag(0.1, 4))
  
  fitFull <- reparamVB(data = data,
                       lambda = lambda, 
                       model = reparamDerivFull, 
                       alpha = 0.03,
                       S = 10, 
                       dimTheta = 4,
                       maxIter = 1000,
                       threshold = 0.0001 * N,
                       RQMC = FALSE)  
  
  UF <- matrix(fitFull$lambda[5:20], 4)
  SigmaF <- t(UF) %*% UF
  
  lambda <- c(1, 0.5, 0.5, diag(0.1, 3))
  fitRes <- reparamVB(data = data,
                     lambda = lambda, 
                     model = reparamDerivRestricted,
                     alpha = 0.03,
                     S = 10, 
                     dimTheta = 3,
                     maxIter = 1000,
                     threshold = 0.0001 * N,
                     RQMC = FALSE)  
  
  UR <- matrix(fitRes$lambda[4:12], 3)
  SigmaR <- t(UR) %*% UR
  
  tFull <- rmvnorm(100, fitFull$lambda[1:4], SigmaF)
  tRes <- rmvnorm(100, fitRes$lambda[1:3], SigmaR)
  
  vbDiff <- mean(sapply(1:100, function(i){
    logJoint1 <- densityNested(data, tFull[i,]);
    logJoint2 <- densityRestricted(data, tRes[i,]);
    logPost1 <- dmvnorm(tFull[i, ], fitFull$lambda[1:4], SigmaF, log = TRUE);
    logPost2 <- dmvnorm(tRes[i, ], fitRes$lambda[1:3], SigmaR, log = TRUE);
    logJoint2 - logPost2 - (logJoint1 - logPost1)
  }))
  
  results <- rbind(results,
                   data.frame(
                     id = iter,
                     N = N, 
                     mcmcDiff = mcmcDiff,
                     vbDiff = vbDiff))
  
}

write.csv(results, paste0('vbetel/id', iter, '.csv'), row.names = FALSE)


