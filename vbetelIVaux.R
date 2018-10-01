hLinReg <- function(lam, g){
  mean(exp(lam %*% t(g)))
}

densityM1 <- function(data, theta){
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

densityM2 <- function(data, theta){
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
  g3 <- w * g1
  g <- cbind(g1, g2, g3)
  
  lambdaHat <- nlm(hLinReg, rep(0, 3), g = g)$estimate
  
  exponent <- c(exp(lambdaHat %*% t(g)))
  prior + sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent))
}

reparamDerivM1 <- function(data, lambda, epsilon){
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
  gradients <- VBfuns::vbetelMatrixCalculations(g, dh2dlam2, lambdaHat, exponent, dgdt)
  prior <- gradScaleT(rep(5, 3), rep(2.5, 3), theta)
  
  dpdt <- gradients$grad + prior$grad
  logp <- gradients$val + prior$val
  
  # Calculate the derivative of the  ELBO
  dELBO <- dtdl %*% dpdt + djdl
  
  list(grad = dELBO, val = logp)
}

reparamDerivM2 <- function(data, lambda, epsilon){
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
  g3 <- w * g1
  g <- cbind(g1, g2, g3)
  
  minima <- nlm(hLinReg, rep(0, 3), g = g, hessian = TRUE)
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  dgdt <- vapply(1:n,
                 function(i) matrix(c(-1, -x[i], -w[i],
                                      -z1[i], -z1[i] * x[i], -z1[i] * w[i],
                                      -w[i], -w[i] * x[i], -w[i]^2),
                                    nrow = 3,
                                    byrow = TRUE),
                 matrix(runif(9), nrow = 3))
  
  # export numerical work to C++
  gradients <- VBfuns::vbetelMatrixCalculations(g, dh2dlam2, lambdaHat, exponent, dgdt)
  dpdt <- gradients$grad
  logp <- gradients$val
  
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
  
  