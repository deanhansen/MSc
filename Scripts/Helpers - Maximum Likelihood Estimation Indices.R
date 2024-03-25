# Author(s): Dean Hansen

# Bivariate Bernoulli
fi2_bernoulli <- function(params) {
  p00 <- params[1]
  p10 <- params[2]
  p01 <- params[3]
  p11 <- 1-sum(p00,p10,p01)
  mu1 <- sum(p00,p10)
  mu2 <- sum(p01,p11)
  a <- p00*p11 - p01*p10
  output <- (mu1*mu2-2*(a^2)*sqrt(mu1*mu2)) / (mu1*mu2 - (a^2)) / 2
  return(output)
}

gdi_bernoulli <- function(params) {
  p00 <- params[1]
  p10 <- params[2]
  p01 <- params[3]
  p11 <- 1-sum(p00,p10,p01)
  mu1 <- sum(p00,p10)
  mu2 <- sum(p01,p11)
  a <- p00*p11 - p01*p10
  output <- 1 + (2*a*sqrt(mu1*mu2) - (mu1^3 + mu2^3)) / (mu1^2 + mu2^2)
  return(output)
}

mdi_bernoulli <- function(params) {
  p00 <- params[1]
  p10 <- params[2]
  p01 <- params[3]
  p11 <- 1-sum(p00,p10,p01)
  mu1 <- sum(p00,p10)
  mu2 <- sum(p01,p11)
  a <- 0
  output <- 1 + (2*a*sqrt(mu1*mu2) - (mu1^3 + mu2^3)) / (mu1^2 + mu2^2)
  return(output)
}

# Author(s): Dean Hansen
# Bivariate Poisson
fi2_poisson <- function(params) {
  l1 <- params[1]
  l2 <- params[2]
  l3 <- params[3]
  sigma <- matrix(c(l1 + l3, l3, l3, l2 + l3), byrow=TRUE, nrow=2)
  sigma_det <- (l1 + l3) * (l2 + l3) - l3^2
  sigma_inv <- (1/sigma_det) * matrix(c(l2 + l3, -l3, -l3, l1 + l3), byrow=TRUE,nrow=2)
  fi <- matrix(c(1,1))
  output <- (t(fi) %*% (sigma * sigma_inv) %*% fi) / 2
  return(output)
}

gdi_poisson <- function(params) {
  l1 <- params[1]
  l2 <- params[2]
  l3 <- params[3]
  mu1 <- l1 + l3
  mu2 <- l2 + l3
  output <- 1 + (2*mu1*mu2*l3) / (mu1^2 + mu2^2)
  return(output)
}

mdi_poisson <- function(params) {
  l1 <- params[1]
  l2 <- params[2]
  l3 <- params[3]
  mu1 <- l1 + l3
  mu2 <- l2 + l3
  output <- (mu1^2 + mu2^2) / (mu1^2 + mu2^2)
  return(output)
}

# Author(s): Dean Hansen
# Bivariate Negative Binomial
fi2_negative_binomial <- function(params) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]
  b1 <- params[4]
  b2 <- params[5]
  w <- (a0+a1)*(a0+a2)*(b1+1)*(b2+1)
  output <- ((w*(2+b1+b2)-2*(a0^2)*b1*b2*sqrt(b1+1)*sqrt(b2+1)) / (w-(a0^2)*b1*b2)) * (1/2)
  return(output)
}

gdi_negative_binomial <- function(params) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]
  b1 <- params[4]
  b2 <- params[5]
  mu1 <- a0 + a1
  mu2 <- a0 + a2
  output <- 1 + ((mu1^2)*b1 + (mu2)^2*b2 + 2*a0*b1*b2*sqrt(mu1)*sqrt(mu2)) / (mu1^2 + mu2^2)
  return(output)
}

mdi_negative_binomial <- function(params) {
  a0 <- params[1]
  a1 <- params[2]
  a2 <- params[3]
  b1 <- params[4]
  b2 <- params[5]
  mu1 <- a0 + a1
  mu2 <- a0 + a2
  output <- 1 + ((mu1^2)*b1 + (mu2)^2*b2) / (mu1^2 + mu2^2)
  return(output)
}

# Author(s): Dean Hansen
# Bivariate Pólya-Aeppli
fi2_polya_aeppli <- function(params) {
  l1 <- params[1]
  l2 <- params[2]
  l3 <- params[3]
  rho <- params[4]
  output <- (1+rho)/(1-rho)
  return(output)
}

gdi_polya_aeppli <- function(params) {
  l1 <- params[1]
  l2 <- params[2]
  l3 <- params[3]
  rho <- params[4]
  mu1 <- (l1 + l3) / (1-rho)
  mu2 <- (l2 + l3)  / (1-rho)
  a <- l3 / ((1-rho)^2)
  output <- 1 + 2 * (rho*(mu1^2+mu2^2) + a * sqrt((l1+l3)*(l2+l3))) / ((1-rho)*(mu1^2+mu2^2))
  return(output)
}

mdi_polya_aeppli <- function(params) {
  l1 <- params[1]
  l2 <- params[2]
  l3 <- params[3]
  rho <- params[4]
  output <- 1 + ((2*rho)/ (1-rho))
  return(output)
}

# Author(s): Dean Hansen
# Wrapper for MLEs
mle_indices <- function(params, dist) {
  if (dist == "BERN") {
    output <- c(fi2_bernoulli(params), gdi_bernoulli(params), mdi_bernoulli(params))
  }
  if (dist == "POIS") {
    output <- c(fi2_poisson(params), gdi_poisson(params), mdi_poisson(params))
  }
  if (dist == "NB") {
    output <- c(fi2_negative_binomial(params), gdi_negative_binomial(params), mdi_negative_binomial(params))
  }
  if (dist == "PA") {
    output <- c(fi2_polya_aeppli(params), gdi_polya_aeppli(params), mdi_polya_aeppli(params))
  }
  return(output)
}

