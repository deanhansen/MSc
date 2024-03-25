# Author(s): Dean Hansen
# Cross-tabulate bivariate count data
crosstabulate <- function(data) {
  x <- data[,1]
  y <- data[,2]
  xmax <- max(x)
  ymax <- max(y)
  xlevels <- as.character(0:xmax)
  ylevels <- as.character(0:ymax)
  xfactor <- factor(x, levels = xlevels)
  yfactor <- factor(y, levels = ylevels)
  xytable <- table(xfactor,yfactor,dnn = c("N1","N2"))
  return(xytable)
}

# Author(s): Dean Hansen
# Summary statistics of the simulated dataset
statistics <- function(data) {
  sim_range <- apply(data,2,range)
  sim_mean <- apply(data,2,mean) |> round(2)
  sim_var <- apply(data,2,var) |> round(2)
  sim_di <- (sim_var / sim_mean) |> round(2)
  output <- list(
    sim_range=sim_range,
    sim_mean=sim_mean,
    sim_var=sim_var,
    sim_di=sim_di
  )
  return(output)
}

# Author(s): Ye Yang
# Generate random sample from Bivariate Pólya-Aeppli Distribution
rbivpa <- function(params,m) {
  BivPA <- matrix(NA,ncol=2,nrow=m)
  for (i in 1:m) {
    Z1 <- rpois(1,params[1])
    Z2 <- rpois(1,params[2])
    Z3 <- rpois(1,params[3])
    n1 <- rgeom(Z1+Z3,1-params[4])
    n2 <- rgeom(Z2+Z3,1-params[4])
    N1 <- sum(n1)+Z1+Z3
    N2 <- sum(n2)+Z2+Z3
    BivPA[i,1] <- N1
    BivPA[i,2] <- N2
  }
  output <- BivPA
  return(output)
}

# Author(s): Ye Yang, Dean Hansen
# Log-likelihood of Bivariate Pólya-Aeppli
evalLL <- function(params,data) {
  x <- nrow(data)
  y <- ncol(data)
  output <- sum(data*log(pmf.M(x=y-1,y=x-1,params[1],params[2],params[3],params[4])))
  return(output)
}

# Author(s): Dean Hansen
# Confidence Intervals for Indices
confint_indices <- function(true_params, true_indices, indices, varcov, dist) {
  q <- qnorm(0.975)
  true_indices_var <- mle_indices_var(true_params,varcov,dist)
  true_indices_sd <- sqrt(true_indices_var)
  confint_fi2 <- true_indices[1] + q*c(-1,1)*true_indices_sd[1]
  coverage_prob_fi2 <- ((indices[1] >= confint_fi2[1]) && (indices[1] <= confint_fi2[2]))
  confint_gdi <- true_indices[2] + q*c(-1,1)*true_indices_sd[2]
  coverage_prob_gdi <- ((indices[2] >= confint_gdi[1]) && (indices[2] <= confint_gdi[2]))
  confint_mdi <- true_indices[3] + q*c(-1,1)*true_indices_sd[3]
  coverage_prob_mdi <- ((indices[3] >= confint_mdi[1]) && (indices[3] <= confint_mdi[2]))
  output <- list(
    confint_fi2=confint_fi2,
    coverage_prob_fi2=coverage_prob_fi2,
    confint_gdi=confint_gdi,
    coverage_prob_gdi=coverage_prob_gdi,
    confint_mdi=confint_mdi,
    coverage_prob_mdi=coverage_prob_mdi
    )
  return(output)
}

# Author(s): Dean Hansen
# Wrapper for fitting MLEs
mle <- function(data, params=NULL, dist, show=FALSE) {
  if (dist == "POIS") {
    output <- bivpois::bp.mle(data)
  }
  if (dist == "PA" && (is.null(params) == FALSE)) {
    output <- NR.F(params, data, show)
  }
  if (dist == "PA" && (is.null(params) == TRUE)) {
    params <- MofM(data)
    output <- NR.F(params, data, show)
  }
  if (dist == "NB") {
    output <- bzinb::bnb(data[,1], data[,2], maxiter=1e5, tol=1e-4, showFlag=show)
  }
  return(output)
}

# Author(s): Dean Hansen
# Construct varcov for Bivariate Poisson distribution
varcov_POIS <- function(params) {
  params <- as.vector(params)
  L11 <- params[1] + params[3]
  L22 <- params[2] + params[3]
  L12 <- params[3]
  output <- matrix(c(L11,L12,L12,L22),nrow=2,ncol=2)
  return(output)
}

# Author(s): Dean Hansen
# Function to run the simulation study given data
simulation_study <- function(B,n,data) {
  output <- matrix(NA,nrow=B,ncol=32)
  successes <- 0
  failures <- 0
  N <- nrow(data)
  
  # run the bootstrap loop B times
  while (successes < B) {
    
    # create bootstrap sample
    rows <- sample(N,n,replace=TRUE)
    bsample <- data[rows,]
    bsample_tab <- crosstabulate(bsample)
    
    # find the MLEs for negative binomial and polya-aeppli
    warning_message <- tryCatch({
      boot_PA <- mle(data = bsample_tab, dist = "PA", show = FALSE)
      boot_NB <- mle(data = bsample, dist = "NB", show = FALSE)
      NULL
    }, warning = function(w) {
      return(conditionMessage(w))
    })
    
    # check PA
    if (is.null(warning_message)) {
      
      # start PA
      boot_mle_PA <- as.vector(boot_PA[,1])
      boot_varcov_PA <- (-1)*solve(hL(boot_mle_PA,bsample_tab))
      boot_varcov_diag_PA <- as.vector(diag(boot_varcov_PA))

      # check varcov is not negative
      if (sum(boot_varcov_diag_PA > 0,na.rm=TRUE) == 4) {
        
        # PA
        boot_loglik_PA <- as.vector(evalLL(boot_mle_PA,bsample_tab))
        boot_indices_PA <- mle_indices(boot_mle_PA,"PA")
        boot_indices_var_PA <- mle_indices_var(boot_mle_PA,boot_varcov_PA,"PA")
        
        # NB
        boot_mle_NB <- as.vector(boot_NB$coefficients[,1])
        boot_varcov_NB <- boot_NB$vcov
        boot_varcov_diag_NB <- as.vector(diag(boot_varcov_NB))
        boot_loglik_NB <- as.vector(boot_NB$lik)
        boot_indices_NB <- mle_indices(boot_mle_NB,"NB")
        boot_indices_var_NB <- mle_indices_var(boot_mle_NB,boot_varcov_NB,"NB")
        
        output[successes+1,] <- c(
          boot_mle_PA,
          boot_varcov_diag_PA,
          boot_loglik_PA,
          boot_indices_PA,
          boot_indices_var_PA,
          boot_mle_NB,
          boot_varcov_diag_NB,
          boot_loglik_NB,
          boot_indices_NB,
          boot_indices_var_NB
        )
        successes <- successes + 1
        cat("Success! Number: ", successes, "\n")
      } 
      else {
        failures <- failures + 1
        cat("Failure! Retrying...Failure:", failures, "\n")
      }
    }
    else {
      failures <- failures + 1
      cat("Failure! Retrying...Failure:", failures, "\n")
    }
  }
  return(output)
}
