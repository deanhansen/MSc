# Author(s): Dean Hansen
# Function to calculate sample indices
sample_indices <- function(x) {
  m <- apply(x, MARGIN = 2, mean)
  v <- apply(x, MARGIN = 2, var)
  fi <- matrix(sqrt(v / m))
  gdi <- (t(sqrt(m)) %*% var(x) %*% sqrt(m)) / (t(m) %*% m)
  mdi <- (t(sqrt(m)) %*% diag(diag(var(x))) %*% sqrt(m)) / (t(m) %*% m)
  fi_p <- (t(fi) %*% (var(x) * solve(var(x))) %*% fi) / length(m)
  output <- list(fi=as.numeric(fi_p), gdi=as.numeric(gdi), mdi=as.numeric(mdi))
  return(output)
}
