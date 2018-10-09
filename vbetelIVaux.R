hLinReg <- function(lam, g){
  mean(exp(lam %*% t(g)))
}

gradScaleT <- function(scale, df, x){
  value <- 0
  grad <- rep(0, length(x))
  for(i in 1:length(x)){
    value <- value + log(1 / scale[i] * dt(x[i] / scale[i], df[i]))
    grad[i] <- -x[i] * (df[i] + 1) / (df[i] * scale[i]^2 + x[i]^2)
  }
  list(val = value, grad = grad)
}

betelMCMC <- function (start, iter, model, mode, covar, burnIn = 0, thin = 1, suppressProgress = FALSE, ...) {
  dim <- length(start)
  nSave <- floor((iter - burnIn)/thin)
  saveDraws <- matrix(0, nSave, dim)
  draw <- start
  invCov <- solve(covar)
  oldDens <- model(draw, ...)
  for (i in 1:iter) {
    if (i == 50) {
      startTime <- Sys.time()
    }
    else if (i == 150) {
      timePerIter <- (Sys.time() - startTime)/100
      class(timePerIter) <- "numeric"
      print(paste0("Estimated Finishing Time: ", Sys.time() + 
                     timePerIter * (iter - 150)))
      if (attr(timePerIter, "units") == "mins") {
        attr(timePerIter, "units") = "secs"
        timePerIter <- timePerIter * 60
      }
    }
    if (!suppressProgress & i%%1000 == 0) {
      mins <- (iter - i) * timePerIter[1]/60
      if (mins > 180) {
        print(paste0("Iteration: ", i, ". Est. Time Remaining: ", 
                     round(mins/60, 2), " hours."))
      }
      else if (mins > 1) {
        print(paste0("Iteration: ", i, ". Est. Time Remaining: ", 
                     round(mins, 1), " minutes."))
      }
      else {
        print(paste0("Iteration: ", i, ". Est. Time Remaining: ", 
                     ceiling(mins * 60), " seconds."))
      }
    }
    
    candidate <- mode + t(chol(covar)) %*% rt(length(draw), 5)
    
    canDens <- model(candidate, ...)
    oldPropDens <- proposalDens(draw, mode, invCov, 5)
    canPropDens <- proposalDens(candidate, mode, invCov, 5)
    ratio <- exp(canDens - oldDens + oldPropDens - canPropDens)
    if (runif(1) < ratio) {
      draw <- candidate
      oldDens <- canDens
    }
    if (i > burnIn & i%%thin == 0) {
      saveDraws[(i - burnIn)/thin, ] <- draw
    }
  }
  draws <- data.frame(saveDraws)
  if (is.null(names(start))) {
    colnames(draws) <- paste0("V", 1:dim)
  }
  else {
    colnames(draws) <- names(start)
  }
  draws$iter <- seq(burnIn + thin, iter, thin)
  draws
}

densityRestricted <- function(data, theta){
  y <- data[,1]
  x <- data[,2]
  z1 <- data[,3]
  z2 <- data[,4]
  w <- data[,5]
  
  alpha <- theta[1]
  beta <- theta[2] 
  delta <- theta[3]
  
  prior <- log(0.2 * dt(0.2 * alpha, 2.5)) + log(0.2 * dt(0.2 * beta, 2.5)) + log(0.2 * dt(0.2 * delta, 2.5))
  
  g1 <- y - alpha - x * beta - w * delta
  g2 <- z1 * g1
  g3 <- z2 * g1
  g4 <- w * g1
  g <- cbind(g1, g2, g3, g4)
  
  lambdaHat <- nlm(hLinReg, rep(0, 4), g = g)$estimate
  
  exponent <- c(exp(lambdaHat %*% t(g)))
  prior + sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent))
}

densityNested <- function(data, theta){
  y <- data[,1]
  x <- data[,2]
  z1 <- data[,3]
  z2 <- data[,4]
  w <- data[,5]
  
  alpha <- theta[1]
  beta <- theta[2] 
  delta <- theta[3]
  nu <- theta[4]
  
  prior <- gradScaleT(rep(5, 4), rep(2.5, 4), theta)$val
  
  g1 <- y - alpha - x * beta - w * delta
  g2 <- z1 * g1
  g3 <- z2 * g1 - nu
  g4 <- w * g1
  g <- cbind(g1, g2, g3, g4)
  
  lambdaHat <- nlm(hLinReg, rep(0, 4), g = g)$estimate
  
  exponent <- c(exp(lambdaHat %*% t(g)))
  prior + sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent))
}

proposalDens <- function(draw, mean, InvCov, df){
  p <- length(draw)
  lgamma((df + p)/2) - lgamma(df/2) - p/2 * log(df * 3.14159) + 0.5 * log(det(InvCov)) -
    -(df + p)/2 * log(1 + 1/df * t(draw - mean) %*% InvCov %*% (draw - mean))
}



negLogETEL <- function(theta, data){
  y <- data[,1]
  x <- data[,2]
  z1 <- data[,3]
  z2 <- data[,4]
  w <- data[,5]
  
  alpha <- theta[1]
  beta <- theta[2] 
  delta <- theta[3]
  nu <- theta[4]
  
  g1 <- y - alpha - x * beta - w * delta
  g2 <- z1 * g1
  g3 <- z2 * g1 - nu
  g4 <- w * g1
  g <- cbind(g1, g2, g3, g4)
  
  lambdaHat <- nlm(hLinReg, rep(0, 4), g = g)$estimate
  
  exponent <- c(exp(lambdaHat %*% t(g)))
  -(sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent)))
}

negLogETELres <- function(theta, data){
  y <- data[,1]
  x <- data[,2]
  z1 <- data[,3]
  z2 <- data[,4]
  w <- data[,5]
  
  alpha <- theta[1]
  beta <- theta[2] 
  delta <- theta[3]
  nu <- theta[4]
  
  g1 <- y - alpha - x * beta - w * delta
  g2 <- z1 * g1
  g3 <- z2 * g1
  g4 <- w * g1
  g <- cbind(g1, g2, g3, g4)
  
  lambdaHat <- nlm(hLinReg, rep(0, 4), g = g)$estimate
  
  exponent <- c(exp(lambdaHat %*% t(g)))
  -(sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent)))
}

reparamDerivFull <- function(data, lambda, epsilon){
  n <- nrow(data)
  y <- data[,1]
  x <- data[,2]
  z1 <- data[,3]
  z2 <- data[,4]
  w <- data[,5]
 
  
  d <- length(epsilon)
  # Transform
  U <- matrix(lambda[(d+1):length(lambda)], d)
  theta <- lambda[1:d] + t(U) %*% epsilon
  
  alpha <- theta[1]
  beta <- theta[2] 
  delta <- theta[3]
  nu <- theta[4]
  
  # Set of transform derivatives
  dtdl <- matrix(0, length(lambda), d)
  diag(dtdl[1:d, 1:d]) <- 1
  
  for(i in 1:d){
    for(j in 1:i){
      dtdl[i*d+j,i] <- epsilon[j]
    }
  }
  
  djdl <- rep(0, length(lambda))
  djdl[seq(d+1, d*(d+1), d+1)] <- 1 / lambda[seq(d+1, d*(d+1), d+1)]
  
  # Calculate dpdt components
  g1 <- y - alpha - x * beta - w * delta
  g2 <- z1 * g1
  g3 <- z2 * g1 - nu
  g4 <- w * g1
  g <- cbind(g1, g2, g3, g4)
  
  minima <- nlm(hLinReg, rep(0, 4), g = g, hessian = TRUE)
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  dgdt <- vapply(1:n,
                 function(i) matrix(c(-1, -x[i], -w[i], 0, 
                                      -z1[i], -z1[i] * x[i], -z1[i] * w[i], 0,
                                      -z2[i], -z2[i] * x[i], -z2[i] * w[i], -1,
                                      -w[i], -w[i] * x[i], -w[i]^2, 0),
                                    nrow = 4,
                                    byrow = TRUE),
                 matrix(runif(16), nrow = 4))
  
  # export numerical work to C++
  gradients <- vbetelMatrixCalculations(g, dh2dlam2, lambdaHat, exponent, dgdt)
  prior <- gradScaleT(rep(5, d), rep(2.5, d), theta)
  
  dpdt <- gradients$grad + prior$grad
  logp <- gradients$val + prior$val
  
  # Calculate the derivative of the  ELBO
  dELBO <- dtdl %*% dpdt + djdl
  
  list(grad = dELBO, val = logp)
}

reparamDerivRestricted <- function(data, lambda, epsilon){
  n <- nrow(data)
  y <- data[,1]
  x <- data[,2]
  z1 <- data[,3]
  z2 <- data[,4]
  w <- data[,5]
  
  
  d <- length(epsilon)
  # Transform
  U <- matrix(lambda[(d+1):length(lambda)], d)
  theta <- lambda[1:d] + t(U) %*% epsilon
  
  alpha <- theta[1]
  beta <- theta[2] 
  delta <- theta[3]
  
  # Set of transform derivatives
  dtdl <- matrix(0, length(lambda), d)
  diag(dtdl[1:d, 1:d]) <- 1
  
  for(i in 1:d){
    for(j in 1:i){
      dtdl[i*d+j,i] <- epsilon[j]
    }
  }
  
  djdl <- rep(0, length(lambda))
  djdl[seq(d+1, d*(d+1), d+1)] <- 1 / lambda[seq(d+1, d*(d+1), d+1)]
  
  # Calculate dpdt components
  g1 <- y - alpha - x * beta - w * delta
  g2 <- z1 * g1
  g3 <- z2 * g1
  g4 <- w * g1
  g <- cbind(g1, g2, g3, g4)
  
  minima <- nlm(hLinReg, rep(0, 4), g = g, hessian = TRUE)
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  dgdt <- vapply(1:n,
                 function(i) matrix(c(-1, -x[i], -w[i],
                                      -z1[i], -z1[i] * x[i], -z1[i] * w[i],
                                      -z2[i], -z2[i] * x[i], -z2[i] * w[i],
                                      -w[i], -w[i] * x[i], -w[i]^2),
                                    nrow = 4,
                                    byrow = TRUE),
                 matrix(runif(12), nrow = 4))
  
  # export numerical work to C++
  gradients <- vbetelMatrixCalculations(g, dh2dlam2, lambdaHat, exponent, dgdt)
  prior <- gradScaleT(rep(5, d), rep(2.5, d), theta)
  dpdt <- gradients$grad + prior$grad
  logp <- gradients$val + prior$val
  
  # Calculate the derivative of the  ELBO
  dELBO <- dtdl %*% dpdt + djdl
  
  list(grad = dELBO, val = logp)
}

qMix <- function(q, m1, m2, sd1, sd2, p){
  support <- seq(-4, 3, 0.001)
  dens1 <- dnorm(support, m1, sd1)
  dens2 <- dnorm(support, m2, sd2)
  densTotal <- p * dens1 + (1 - p) * dens2
  CDF <- cumsum(densTotal) / sum(densTotal)
  quantiles <- numeric(length(q))
  for(k in 1:length(q)){
    x0 <- CDF[max(which(CDF < q[k]))]
    x1 <- CDF[min(which(CDF > q[k]))]
    y0 <- support[max(which(CDF < q[k]))]
    y1 <- support[min(which(CDF > q[k]))]
    quantiles[k] <- y0 + (q[k] - x0) * (y1 - y0) / (x1 - x0)
  }
  quantiles
}
  
  
