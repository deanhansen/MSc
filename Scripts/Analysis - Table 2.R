# Author(s): Dean Hansen
# Table 2 from Kokonendji, Puig (2018)
mu <- c(11965.52, 44455.59, 11654.23)
vars <- c(96707437, 2550265224, 109833650)
d <- sqrt(vars/mu)
mu_sqrt <- sqrt(mu)
var_matrix <- matrix(0,nrow=3,ncol=3)
corr_matrix <- matrix(c(1,0.381,0.565,0.381,1,0.388,0.565,0.388,1),nrow=3)
diag(var_matrix) <- vars
Sigma <- sqrt(var_matrix) %*% corr_matrix %*% sqrt(var_matrix)
Sigma_inv <- solve(Sigma)
fi3 <- (t(d) %*% (Sigma * Sigma_inv) %*% d) / 3
gdi <- (mu_sqrt %*% Sigma %*% mu_sqrt) / (mu %*% mu)
mdi <- (mu_sqrt %*% diag(diag(Sigma)) %*% mu_sqrt) / (mu %*% mu)
c(fi3, gdi, mdi)
