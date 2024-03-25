# Author(s): Dean Hansen

library('extraDistr')
library('VGAM')
source("Main/Scripts/Probability Mass Functions.R")
palette(viridis::viridis(100))
a <- 1600

# Bivariate Poisson ---------------------------------------------------------
gen_biv_poisson <- function(mu1, mu2, mu3, t) {
x1 <- 0:t
x2 <- 0:t
z <- function(x1,x2) {z <- extraDistr::dbvpois(x1, x2, a = mu1, b = mu2, c = mu3)}
f <- t(outer(x1,x2,z))
return(list(x1=x1, x2=x2,f=f))
}

p1 <- gen_biv_poisson(2,3,1,15)
p2 <- gen_biv_poisson(2,3,2,15)
p3 <- gen_biv_poisson(2,3,3,15)
p4 <- gen_biv_poisson(2,3,4,15)

png("main/plots/Bivariate Densities/bivariate_poisson.png", pointsize=10, width=a, height=a, res=300) 
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1))
image(p1$f, xlim = c(0, 15), ylim = c(0,15), x = p1$x1, y = p1$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 1)), cex.main = 1.25, col.main = "black")
image(p2$f, xlim = c(0, 15), ylim = c(0,15), x = p2$x1, y = p2$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 2)), cex.main = 1.25, col.main = "black")
image(p3$f, xlim = c(0, 15), ylim = c(0,15), x = p3$x1, y = p3$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 3)), cex.main = 1.25, col.main = "black")
image(p4$f, xlim = c(0, 15), ylim = c(0,15), x = p4$x1, y = p4$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 4)), cex.main = 1.25, col.main = "black")
dev.off()

jpeg("main/plots/Bivariate Densities/bivariate_poisson.jpeg", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1))
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1))
image(p1$f, xlim = c(0, 15), ylim = c(0,15), x = p1$x1, y = p1$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 1)), cex.main = 1.25, col.main = "black")
image(p2$f, xlim = c(0, 15), ylim = c(0,15), x = p2$x1, y = p2$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 2)), cex.main = 1.25, col.main = "black")
image(p3$f, xlim = c(0, 15), ylim = c(0,15), x = p3$x1, y = p3$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 3)), cex.main = 1.25, col.main = "black")
image(p4$f, xlim = c(0, 15), ylim = c(0,15), x = p4$x1, y = p4$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 4)), cex.main = 1.25, col.main = "black")
dev.off()

# Bivariate Negative Binomial ---------------------------------------------------------
gen_biv_negative_binomial <- function(a0, a1, a2, b1, b2, t) {
  x1 <- 0:t
  x2 <- 0:t
  z <- function(x1,x2) {
    output <- matrix(NA,nrow=length(x1),ncol=length(x2))
    for (k in x1) {
      for (m in x2) {
        output[k+1,m+1] <- dbvnb(k,m,a0=a0,a1=a1,a2=a2,b1=b1,b2=b2)
      }
    }
    return(output)
  }
  f <- t(z(x1,x2))
  return(list(x1=x1,x2=x2,f=f))
}

p1 <- gen_biv_negative_binomial(a0=5,a1=5,a2=8,b1=0.25,b2=0.50,t=15)
p2 <- gen_biv_negative_binomial(a0=5,a1=5,a2=2,b1=0.50,b2=0.25,t=15)
p3 <- gen_biv_negative_binomial(a0=5,a1=8,a2=8,b1=0.25,b2=0.25,t=15)
p4 <- gen_biv_negative_binomial(a0=5,a1=8,a2=2,b1=0.25,b2=0.50,t=15)

png("main/plots/Bivariate Densities/bivariate_negative_binomial.png", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1), cex=0.8)
image(p1$f, xlim = c(0, 15), ylim = c(0,15), x = p1$x1, y = p1$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 5, 8)", ", ", bold(beta) == "(0.25, 0.50)")), cex.main = 1.25, col.main = "black")
image(p2$f, xlim = c(0, 15), ylim = c(0,15), x = p2$x1, y = p2$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 5, 2)", ", ", bold(beta) == "(0.50, 0.25)")), cex.main = 1.25, col.main = "black")
image(p3$f, xlim = c(0, 15), ylim = c(0,15), x = p3$x1, y = p3$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 8, 8)", ", ", bold(beta) == "(0.25, 0.25)")), cex.main = 1.25, col.main = "black")
image(p4$f, xlim = c(0, 15), ylim = c(0,15), x = p4$x1, y = p4$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 8, 2)", ", ", bold(beta) == "(0.25, 0.50)")), cex.main = 1.25, col.main = "black")
dev.off()

jpeg("main/plots/Bivariate Densities/bivariate_negative_binomial.jpeg", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1), cex=0.8)
image(p1$f, xlim = c(0, 15), ylim = c(0,15), x = p1$x1, y = p1$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 5, 8)", ", ", bold(beta) == "(0.25, 0.50)")), cex.main = 1.25, col.main = "black")
image(p2$f, xlim = c(0, 15), ylim = c(0,15), x = p2$x1, y = p2$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 5, 2)", ", ", bold(beta) == "(0.50, 0.25)")), cex.main = 1.25, col.main = "black")
image(p3$f, xlim = c(0, 15), ylim = c(0,15), x = p3$x1, y = p3$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 8, 8)", ", ", bold(beta) == "(0.25, 0.25)")), cex.main = 1.25, col.main = "black")
image(p4$f, xlim = c(0, 15), ylim = c(0,15), x = p4$x1, y = p4$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(bold(alpha) == "(5, 8, 2)", ", ", bold(beta) == "(0.25, 0.50)")), cex.main = 1.25, col.main = "black")
dev.off()

# Bivariate Polya-Aeppli ---------------------------------------------------------
gen_biv_polya_aeppli <- function(L1, L2, L3, rhoo, t) {
  x1 <- 0:t
  x2 <- 0:t
  z <- pmf.M(t,t,L1,L2,L3,rhoo)
  return(list(x1=x1, x2=x2,f=z))
}

p1 <- gen_biv_polya_aeppli(L1=2, L2=1, L3=1, rhoo=0.25, t=15)
p2 <- gen_biv_polya_aeppli(L1=2, L2=3, L3=1, rhoo=0.50, t=15)
p3 <- gen_biv_polya_aeppli(L1=2, L2=1, L3=3, rhoo=0.25, t=15)
p4 <- gen_biv_polya_aeppli(L1=2, L2=3, L3=3, rhoo=0.50, t=15)

png("main/plots/Bivariate Densities/bivariate_polya_aeppli.png", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1))
image(p1$f, xlim = c(0, 15), ylim = c(0,15), x = p1$x1, y = p1$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 1, ", ", lambda[3] == 1, ", ", rho == 0.25)), cex.main = 1.25, col.main = "black")
image(p2$f, xlim = c(0, 15), ylim = c(0,15), x = p2$x1, y = p2$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 1, ", ", rho == 0.50)), cex.main = 1.25, col.main = "black")
image(p3$f, xlim = c(0, 15), ylim = c(0,15), x = p3$x1, y = p3$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 1, ", ", lambda[3] == 3, ", ", rho == 0.25)), cex.main = 1.25, col.main = "black")
image(p4$f, xlim = c(0, 15), ylim = c(0,15), x = p4$x1, y = p4$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 3, ", ", rho == 0.50)), cex.main = 1.25, col.main = "black")
dev.off()

png("main/plots/Bivariate Densities/bivariate_polya_aeppli.jpeg", pointsize=10, width=a, height=a, res=300)
par(mfrow = c(2, 2), mar=c(3, 3, 2, 2), oma = c(1, 1, 1, 1))
image(p1$f, xlim = c(0, 15), ylim = c(0,15), x = p1$x1, y = p1$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 1, ", ", lambda[3] == 1, ", ", rho == 0.25)), cex.main = 1.25, col.main = "black")
image(p2$f, xlim = c(0, 15), ylim = c(0,15), x = p2$x1, y = p2$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 1, ", ", rho == 0.50)), cex.main = 1.25, col.main = "black")
image(p3$f, xlim = c(0, 15), ylim = c(0,15), x = p3$x1, y = p3$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 1, ", ", lambda[3] == 3, ", ", rho == 0.25)), cex.main = 1.25, col.main = "black")
image(p4$f, xlim = c(0, 15), ylim = c(0,15), x = p4$x1, y = p4$x2, xlab = "", ylab = "", col = palette())
title(main = expression(paste(lambda[1] == 2, ", ", lambda[2] == 3, ", ", lambda[3] == 3, ", ", rho == 0.50)), cex.main = 1.25, col.main = "black")
dev.off()
