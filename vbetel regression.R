library(tidyverse)
library(VBfuns)

N <- 250
z <- rnorm(N, 0.5, 1)
alpha <- 0
beta <- 1
y <- alpha + beta * z
for(i in 1:N){
  u <- runif(1)
  if(u < 0.5){
    y[i] <- y[i] + rnorm(1, 0.75, 0.75)
  } else {
    y[i] <- y[i] + rnorm(1, -0.75, 1.25)
  }
}
nu <- 0.5 * (3 * 0.75 * 0.75^2 + 0.75^3) + 0.5 * (3 * (-0.75) * 1.25^2 + (-0.75)^3)
var <- 0.5 * (0.75^2 + 0.75^2) + 0.5 * ((-0.75)^2 + 1.25^2)
lam <- c(1, 2, 3)

hLinReg <- function(lam, g){
  mean(exp(lam %*% t(g)))
}

density <- function(data, theta){
  y <- data[,1]
  z <- data[,2]
  alpha <- theta[1]
  beta <- theta[2] 
  nu <- theta[3]
  sig <- theta[4]
  
  g1 <- y - alpha - beta * z
  g2 <- z * g1
  g3 <- g1^3 - nu
  g4 <- g1^2 - sig
  g <- cbind(g1, g2, g3, g4)
  
  lambdaHat <- nlm(hLinReg, rep(0, 4), g = g)$estimate
  

  exponent <- c(exp(lambdaHat %*% t(g)))
  sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent))
  
#  list(dens = sum(lambdaHat %*% t(g)) - length(y) * log(sum(exponent)),
 #      lambdaHat = lambdaHat)
}

draw <- c(alpha = 0, beta = 1, nu = -1, `sigma^{2}` = 2)
betelMCMC <- metropolisHastings(start = draw, 
                                iter = 10000, 
                                model = density,
                                data = cbind(y, z))


betelMCMC %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free', labeller = label_parsed)
  

betelMCMC %>%
  filter(iter > 2500) %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(draw), stat = 'density') + facet_wrap(~var, scales = 'free', labeller = label_parsed)


trueDensity <- function(data, theta){
  y <- data[,1]
  z <- data[,2]
  dens <- 0.5 * 1 / sqrt(2 * pi * 0.75^2) * exp(-(y - theta[1] - theta[2] * z - 0.75)^2 / (2 * 0.75^2)) +
    0.5 * 1 / sqrt(2 * pi * 1.25^2) *  exp(-(y - theta[1] - theta[2] * z + 0.75)^2 / (2 * 1.75^2))
  sum(log(dens))
}

trueMCMC <- metropolisHastings(start = c(alpha = 0, beta = 1),
                               iter = 10000, 
                               model = trueDensity,
                               data = cbind(y, z))


trueMCMC %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(iter, draw)) + facet_wrap(~var, scales = 'free')

trueMCMC %>%
  filter(iter > 2500) %>%
  gather(var, draw, -iter) %>%
  ggplot() + geom_line(aes(draw), stat = 'density') + facet_wrap(~var, scales = 'free')

MCMCdraws <-
  rbind(trueMCMC %>%
          filter(iter > 2500) %>%
          gather(var, draw, -iter) %>%
          mutate(method = 'MCMC-true'),
        betelMCMC %>%
          filter(iter > 2500) %>%
          gather(var, draw, -iter) %>%
          mutate(method = 'MCMC-empirical'))

ggplot(MCMCdraws) + geom_line(aes(draw, colour = method), stat = 'density') + 
  geom_vline(data = data.frame(true = c(nu, var), method = 'MCMC-true', var = c('nu', 'sigma^{2}')), aes(xintercept = true, colour = method)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed)

reparamDeriv <- function(data, lambda, epsilon){
  y <- data[,1]
  z <- data[,2]
  d <- length(epsilon)
  # Transform
  U <- matrix(lambda[(d+1):length(lambda)], d)
  theta <- lambda[1:d] + t(U) %*% epsilon
  
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
  n <- length(y)
  alpha <- theta[1]
  beta <- theta[2] 
  nu <- theta[3]
  sigma2 <- theta[4]
  
  g1 <- y - alpha - beta * z
  g2 <- z * g1
  g3 <- g1^3 - nu
  g4 <- g1^2 - sigma2
  g <- cbind(g1, g2, g3, g4)
  
  minima <- nlm(hLinReg, rep(0, 4), g = g, hessian = TRUE)
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  dgdt <- vapply(1:n,
                 function(x) matrix(c(-1, z[n], 0, 0,
                                      -z[n], -z[n]^2, 0, 0,
                                      -3 * g1[n]^2, -3 * z[n] * g1[n]^2, -1, 0,
                                      -2 *g1[n], -2 * z[n] * g1[n], 0, -1),
                                    nrow = 4,
                                    byrow = TRUE),
                 matrix(runif(16), 4))
  
  # export numerical work to C++
  gradients <- vbetelMatrixCalculations(g, dh2dlam2, lambdaHat, exponent, dgdt)
  dpdt <- gradients$grad
  logp <- gradients$val
  
  # Calculate the derivative of the  ELBO
  dELBO <- dtdl %*% dpdt + djdl
  
  list(grad = dELBO, val = logp)
}

lambda <- c(0, 0, 0, 1, diag(0.1, 4))

fit <- gaussianVB(data = cbind(y, z), 
                  lambda = lambda, 
                  model = reparamDeriv, 
                  alpha = 0.03,
                  S = 5, 
                  dimEpsilon = 4,
                  maxIter = 1000,
                  RQMC = FALSE)


qplot(1:fit$iter, fit$LB[1:fit$iter], geom = 'line')


U <- matrix(fit$lambda[5:20], 4)
Sigma <- t(U) %*% U

vbDens <- gaussianDensity(mu = fit$lambda[1:4], 
                          sigma = sqrt(diag(Sigma)),
                          transform = rep('identity', 4),
                          names = c('alpha', 'beta', 'nu', 'sigma^{2}'))
vbDens$method <- 'VB-Reparam'

ggplot(MCMCdraws) + geom_line(aes(draw, colour = method), stat = 'density') + 
  geom_vline(data = data.frame(true = c(nu, var), method = 'MCMC-true', var = c('nu', 'sigma^{2}')), aes(xintercept = true, colour = method)) + 
  geom_line(data = vbDens, aes(support, density, colour = method)) + 
  facet_wrap(~var, scales = 'free', labeller = label_parsed)


### Stein
steinDeriv <- function(y, z, theta){
  
  # Calculate dpdt components
  n <- length(y)
  alpha <- theta[1]
  beta <- theta[2] 
  nu <- theta[3]
  sigma2 <- theta[4]
  
  g1 <- y - alpha - beta * z
  g2 <- z * g1
  g3 <- g1^3 - nu
  g4 <- g1^2 - sigma2
  g <- cbind(g1, g2, g3, g4)
  minima <- nlm(hLinReg, rep(0, 4), g = g, hessian = TRUE)
  
  lambdaHat <- minima$estimate
  dh2dlam2 <- minima$hessian
  exponent <- c(exp(lambdaHat %*% t(g)))
  
  # export numerical work to C++
  gradients <- matrixCalculations(y, z, alpha, beta, g, dh2dlam2, lambdaHat, exponent)
  dpdt <- gradients$grad
  logp <- gradients$val
  
  list(grad = dpdt, val = logp)
}

N <- 75
particles <- betelMCMC[sample(5001:10000, N), ] + matrix(rnorm(4 * N, 0, 0.35), ncol = 4)
initial <- data.frame(particle = c(particles), var = rep(c('alpha', 'beta', 'nu', 'sigma^{2}'), rep(N, 4)))

ggplot(MCMCdraws) + geom_line(aes(draw, colour = method), stat = 'density') + 
  geom_vline(data = data.frame(true = c(nu, var), method = 'MCMC-true', var = c('nu', 'sigma^{2}')), aes(xintercept = true, colour = method)) + 
  geom_line(data = initial, aes(particle), stat = 'density') + 
  labs(colour = 'MCMC Likelihood') +
  facet_wrap(~var, scales = 'free', labeller = label_parsed)

h <- 0.1

fitStein <- steinVB(particles = particles,
                                 model = steinDeriv,
                                 threshold = 0.01,
                                 alpha = 0.025,
                                 y = y,
                                 z = z)

steinDf <- data.frame(particle = c(fitStein), var = rep(c('alpha', 'beta', 'nu', 'sigma^{2}'), rep(N, 4)), method = 'VB-Stein')

cbPalette <- c("#000000", "#E69F00", "#2684B9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
MCMCdraws %>%
 #filter(var != 'sigma^{2}') %>%
  mutate(method = ifelse(method == 'MCMC-true', 'MCMC-parametric', method)) %>%
  ggplot() + geom_line(aes(draw, colour = method), stat = 'density') +
  geom_vline(data = data.frame(value = c(nu, var), method = 'MCMC-parametric', var = c('nu', 'sigma^{2}')), aes(xintercept = value, colour = method)) + 
 # geom_vline(data = data.frame(value = nu, method = 'MCMC-parametric', var = 'nu'), aes(xintercept = value, colour = method)) + 
  geom_line(data = steinDf , aes(particle, colour = method), stat = 'density') + 
  geom_line(data = vbDens, aes(support, density, colour = method)) +
  theme_bw() + 
  facet_wrap(~var, scales = 'free', label = label_parsed) + 
  scale_colour_manual(values=cbPalette) + 
  theme(legend.position = 'bottom') + 
  labs(x = 'support')
