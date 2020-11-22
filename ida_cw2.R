##############################2(b)##############################
install.packages("maxLik")
library(maxLik)

# Load data
load("dataex2.Rdata")
x = dataex2$X
r = dataex2$R

# Define the log likelihood function
log_like_func <- function(r, x, theta){
  sum(r*dnorm(x, mean = theta, sd = 1.5, log = TRUE) + 
        (1-r)*pnorm(x, mean = theta, sd = 1.5, log = TRUE))
}

# Do the optimisation of log likelihood automatically by using "maxLik" package
mle <- maxLik(logLik = log_like_func, x = x, r=r, start = c(15))
print(mle)

##############################4##############################
# Define the log likelihood function
log_like_func <- function(data, param) {
  x = data[[1]]$X
  y = data[[1]]$Y
  y_obs = y[!is.na(y)]
  y_mis = y[is.na(y)]
  x_o = x[!is.na(y)]
  x_m = x[is.na(y)]
  p <- data[[2]]
  beta_0 <- param[1]
  beta_1 <- param[2]
  exp = exp(beta_0+x_m*beta_1)
  sum(y_obs*(beta_0+x_o*beta_1) - log(1+exp(beta_0+x_o*beta_1))) +
    sum(p*(beta_0+x_m*beta_1) - log(1+exp(beta_0+x_m*beta_1)))
}

# Define the EM algorithm
em.bernoulli <- function(data, beta0, eps){
  
  x = data$X
  y = data$Y
  x_m = x[is.na(y)]
  
  beta <- beta0
  beta_0 <- beta[1]; beta_1 <- beta[2]
  
  diff <- 1
  while(diff > eps){
    beta.old <- beta
    #E-step
    ptilde <- exp(beta_0+x_m*beta_1)/(1+exp(beta_0+x_m*beta_1))
    #M-step
    mle <- maxLik(logLik = log_like_func, data = list(data, ptilde), start = c(10, 10))
    beta_0 <- mle$estimate[1]
    beta_1 <- mle$estimate[2]
    
    beta <- c(beta_0, beta_1)
    diff <- sum(abs(beta - beta.old))
  }
  return(beta)
}

# Load data
load("dataex4.Rdata")

res <- em.bernoulli(dataex4, c(1, 2), 0.00001)
beta_0 <- res[1]
beta_1 <- res[2]
beta_0; beta_1


##############################5##############################
# Define the EM algorithm
em.mixture.lognormal.exponential <- function(y, theta0, eps){
  n <- length(y)
  theta <- theta0
  p <- theta[1]
  mu <- theta[2]; sigma_square <- theta[3]
  lambda <- theta[4]
  diff <- 1
  while(diff > eps){
    theta.old <- theta
    #E-step
    ptilde1 <- p*dlnorm(y, meanlog = mu, sdlog = sqrt(sigma_square))
    ptilde2 <- (1 - p)*dexp(y, lambda)
    ptilde <- ptilde1/(ptilde1 + ptilde2)
    #M-step
    p <- mean(ptilde)
    mu <- sum(log(y)*ptilde)/sum(ptilde)
    sigma_square <- sum(((log(y) - mu)^2)*ptilde)/sum(ptilde)
    lambda <- sum(1-ptilde)/sum((1-ptilde)*y)
    theta <- c(p, mu, sigma_square, lambda)
    diff <- sum(abs(theta - theta.old))
  }
  return(theta)
}

# Load data
load("dataex5.Rdata")
y = dataex5

res <- em.mixture.lognormal.exponential(y = y, c(0.1, 1, 0.5^2, 2), 0.00001)
p <- res[1]
mu <- res[2]
sigma_square <- res[3]
lambda <- res[4]
p; mu; sigma_square; lambda

# Draw the density function
hist(y, main = "The Mixture distribution",
     xlab = "Value",
     ylab="Density", freq = F, cex.main = 1.5,
     cex.lab = 1.5, cex.axis = 1.4, ylim = c(0, 0.2))
curve(p*dlnorm(x, mu, sigma_square) + (1 - p)*dexp(x, lambda),
      add = TRUE, lwd = 2, col = "blue2")
