# Author(s): Dean Hansen
# Poisson Gradients
fi2_poisson_grad <- function(params) {
  L1 <- params[1]
  L2 <- params[2]
  L3 <- params[3]
  gL1 <- 0
  gL2 <- 0
  gL3 <- 0
  output <- matrix(c(gL1,gL2,gL3),nrow=1)
  return(output)
}

gdi_poisson_grad <- function(params) {
  L1 <- params[1]
  L2 <- params[2]
  L3 <- params[3]
  mu1 <- L1+L3
  mu2 <- L2+L3
  D <- mu1^2 + mu2^2
  D_sq <- D^2
  gL1 <- (2*mu2*L3*(mu2^2 - mu1^2))
  gL2 <- (2*mu1*L3*(mu1^2 - mu2^2))
  gL3 <- (2*mu1*mu2*D + 2*L3*(mu1 + mu2)^2*(mu1 - mu2))
  output <- matrix(c(gL1,gL2,gL3),nrow=1) / D_sq
  return(output)
}

mdi_poisson_grad <- function(params) {
  L1 <- params[1]
  L2 <- params[2]
  L3 <- params[3]
  gL1 <- 0
  gL2 <- 0
  gL3 <- 0
  output <- matrix(c(gL1,gL2,gL3),nrow=1)
  return(output)
}

# Author(s): Dean Hansen
# Poisson Variance
fi2_poisson_delta_var <- function(params,varcov) {
  grad <- fi2_poisson_grad(params)
  output <- grad %*% varcov %*% t(grad)
  output <- as.numeric(output)
  return(output)
}

gdi_poisson_delta_var <- function(params,varcov) {
  grad <- gdi_poisson_grad(params)
  output <- grad %*% varcov %*% t(grad)
  output <- as.numeric(output)
  return(output)
}

mdi_poisson_delta_var <- function(params,varcov) {
  grad <- mdi_poisson_grad(params)
  output <- grad %*% varcov %*% t(grad)
  output <- as.numeric(output)
  return(output)
}

# Author(s): Dean Hansen
# Negative Binomial Gradient
fi2_negative_binomial_grad <- function(params) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]
  b1 <- params[4]
  b2 <- params[5]
  mu1 <- a0+a1
  mu2 <- a0+a2
  b1_bar <- b1+1
  b2_bar <- b2+1
  N <- (mu1*mu2*b1_bar*b2_bar*(b1_bar + b2_bar) - 2*a0^2*b1*b2*sqrt(b1_bar*b2_bar))
  D <- mu1*mu2*b1_bar*b2_bar - (a0^2)*b1*b2
  D_sq <- D^2
  ga0 <- (2*mu1*mu2*(b1*b2)^2*(mu1*(b2-b1) + mu2*(b1-b2)))
  ga1 <- (mu2*b1_bar*b2_bar*a0^2*b1*b2*(2*sqrt(b1_bar*b2_bar) - (b1_bar + b2_bar)))
  ga2 <- (mu1*b1_bar*b2_bar*a0^2*b1*b2*(2*sqrt(b1_bar*b2_bar) - (b1_bar + b2_bar)))
  gb1 <- (mu1*mu2*b2_bar*(2*b1_bar + b2_bar) - a0^2*b2*sqrt(b1_bar*b2_bar)*(1-b1*(1/b1_bar)))*D - (N * (mu1*mu2*b2_bar - a0^2*b2))
  gb2 <- (mu1*mu2*b1_bar*(2*b2_bar + b1_bar) - a0^2*b1*sqrt(b1_bar*b2_bar)*(1-b2*(1/b2_bar)))*D - (N * (mu1*mu2*b1_bar - a0^2*b1))
  output <- matrix(c(ga0,ga1,ga2,gb1,gb2),nrow=1) / D_sq
  return(output)
}

gdi_negative_binomial_grad <- function(params) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]
  b1 <- params[4]
  b2 <- params[5]
  mu1 <- a0+a1
  mu2 <- a0+a2
  N <- mu1^2*b1^3 + mu2^2*b2^3 + 2*a0*b1^(3/2)*b2^(3/2)*sqrt(mu1*mu2)
  D <- mu1^2*b1^2 + mu2^2*b2^2
  D_sq <- D^2
  ga0 <- ((2*mu1^(3/2)*sqrt(mu1)*b1^3 + 2*sqrt(mu1)*mu2^(3/2)*b2^3 + b1^(3/2)*b2^(3/2)*(2*mu1*mu2 + a0*(mu1 + mu2)))*D) / (sqrt(mu1*mu2)) - (2*sqrt(mu1*mu2)*(mu1*b1^2 + mu2*b2^2)*N)
  ga1 <- ((2*mu1^(3/2)*b1^3 + a0*b1^(3/2)*b2^(3/2)*sqrt(mu2))*D) / sqrt(mu1) - (2*mu1*b1^2*N)
  ga2 <- ((2*mu2^(3/2)*b2^3 + a0*b1^(3/2)*b2^(3/2)*sqrt(mu1))*D) / sqrt(mu1) - (2*mu2*b2^2*N)
  gb1 <- (3*(mu1^2*b1^2 + a0*sqrt(b1)*b2^(3/2)*sqrt(mu1*mu2)))*D -2*mu1^2*b1*N
  gb2 <- (3*(mu2^2*b2^2 + a0*b1^(3/2)*sqrt(b2)*sqrt(mu1*mu2)))*D -2*mu2^2*b2*N
  output <- matrix(c(ga0,ga1,ga2,gb1,gb2),nrow=1) / D_sq
  return(output)
}

mdi_negative_binomial_grad <- function(params) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]
  b1 <- params[4]
  b2 <- params[5]
  mu1 <- a0+a1
  mu2 <- a0+a2
  N <- mu1^2*b1^3 + mu2^2*b2^3
  D <- mu1^2*b1^2 + mu2^2*b2^2
  D_sq <- D^2
  ga0 <- (2*mu1*mu2*b1^2*b2^2*(mu1*(b2 - b1) + mu2*(b1 - b2)))
  ga1 <- (2*mu1*mu2^2*b1^2*b2^2*(b1 - b2))
  ga2 <- (2*mu1^2*mu2*b1^2*b2^2*(b2 - b1))
  gb1 <- (mu1^2*b1*(mu1^2*b1^3 + mu2^2*b2^2*(3*b1 - 2*b2)))
  gb2 <- (mu2^2*b2*(mu2^2*b2^3 + mu1^2*b1^2*(3*b2 - 2*b1)))
  output <- matrix(c(ga0,ga1,ga2,gb1,gb2),nrow=1) / D_sq
  return(output)
}

# Author(s): Dean Hansen
# Negative Binomial Variance
fi2_negative_binomial_delta_var <- function(params,varcov) {
  grad <- fi2_negative_binomial_grad(params)
  output <- grad %*% varcov %*% t(grad)
  output <- as.numeric(output)
  return(output)
}

gdi_negative_binomial_delta_var <- function(params,varcov) {
  grad <- gdi_negative_binomial_grad(params)
  output <- grad %*% varcov %*% t(grad)
  output <- as.numeric(output)
  return(output)
}

mdi_negative_binomial_delta_var <- function(params,varcov) {
  grad <- mdi_negative_binomial_grad(params)
  output <- grad %*% varcov %*% t(grad)
  output <- as.numeric(output)
  return(output)
}

# Author(s): Dean Hansen
# Polya-Aeppli Gradient
fi2_polya_aeppli_grad <- function(params) {
  L1<-params[1]
  L2<-params[2]
  L3<-params[3]
  rhoo<-params[4]
  gL1 <- 0
  gL2 <- 0
  gL3 <- 0
  gRho <- 4/(1-rhoo)^2
  output <- matrix(c(gL1,gL2,gL3,gRho),nrow=1)
  return(output)
}

gdi_polya_aeppli_grad <- function(params) {
  L1 <- params[1]
  L2 <- params[2]
  L3 <- params[3]
  rhoo <- params[4]
  mu1 <- L1+L3
  mu2 <- L2+L3
  D <- (1-rhoo)*(mu1^2+mu2^2)
  D_sq <- D^2
  gL1 <- sqrt(mu1/mu2)*(2*(1-rhoo)*(mu1^2 - 3*mu2^2))
  gL2 <- sqrt(mu2/mu1)*(2*(1-rhoo)*(mu2^2 - 3*mu1^2))
  gL3 <- (2*mu1*mu2*D + L3*(1-rhoo)*(mu1 - mu2)^2*(mu1 + mu2)) / (sqrt(mu1*mu2))
  gRho <- (2*(mu1^2 + mu2^2) - 2*L3*sqrt(mu1*mu2)) / ((1-rhoo))
  output <- matrix(c(gL1,gL2,gL3,gRho),nrow=1) / D_sq
  return(output)
}

mdi_polya_aeppli_grad <- function(params) {
  L1 <- params[1]
  L2 <- params[2]
  L3 <- params[3]
  rhoo <- params[4]
  gL1 <- 0
  gL2 <- 0
  gL3 <- 0
  gRho <- 2/(1-rhoo)^2
  output <- matrix(c(gL1,gL2,gL3,gRho),nrow=1)
  return(output)
}

# Author(s): Dean Hansen
# Polya Aeppli Variance
fi2_polya_aeppli_delta_var <- function(params,varcov) {
  grad <- fi2_polya_aeppli_grad(params)
  output <- as.numeric(grad %*% varcov %*% t(grad))
  return(output)
}

gdi_polya_aeppli_delta_var <- function(params,varcov) {
  grad <- gdi_polya_aeppli_grad(params)
  output <- as.numeric(grad %*% varcov %*% t(grad))
  return(output)
}

mdi_polya_aeppli_delta_var <- function(params,varcov) {
  grad <- mdi_polya_aeppli_grad(params)
  output <- as.numeric(grad %*% varcov %*% t(grad))
  return(output)
}

# Author(s): Dean Hansen
# Wrapper for Delta Method
mle_indices_var <- function(params, varcov, dist) {
  if (dist == "POIS") {
    output <- c(fi2_poisson_delta_var(params,varcov), gdi_poisson_delta_var(params,varcov), mdi_poisson_delta_var(params,varcov))
  }
  if (dist == "NB") {
    output <- c(fi2_negative_binomial_delta_var(params,varcov), gdi_negative_binomial_delta_var(params,varcov), mdi_negative_binomial_delta_var(params,varcov))
  }
  if (dist == "PA") {
    output <- c(fi2_polya_aeppli_delta_var(params,varcov), gdi_polya_aeppli_delta_var(params,varcov), mdi_polya_aeppli_delta_var(params,varcov))
  }
  return(output)
}
