#' Hessian for Bivariate Pólya-Aeppli
#' Author(s): Ye-Yang, Dean Hansen
#'
#' @param x 
#' @param y 
#' @param params 
#' @param data cross-tabulated dataset
#'
#' @return Hessian matrix over entire dataset
#' @export
#'
#' @examples
#' data <- bivpois::rbp(500,c(1,2,3))
#' data <- crosstabulate(data)
#' x <- ncol(data)-1
#' y <- nrow(data)-1
#' params <- c(1,2,3,0.25)
#' hL(params, data)
hL <- function(params, data) {
  x <- ncol(data) - 1
  y <- nrow(data) - 1
  m <- length(params)
  # Initialize empty hessian matrix
  output <- matrix(NA,ncol=m,nrow=m)
  # Probability mass function
  pmf <- pmf.M(x,y,params[1],params[2],params[3],params[4])
  pmf_reciprocated <- 1 / pmf
  pmf_reciprocated_squared <- 1 / pmf^2
  # First derivatives
  dL1 <- D.L1(x,y,params[1],params[2],params[3],params[4])
  dL2 <- D.L2(x,y,params[1],params[2],params[3],params[4])
  dL3 <- D.L3(x,y,params[1],params[2],params[3],params[4])
  dRho <- D.Rho(x,y,params[1],params[2],params[3],params[4])
  # Second Derivatives
  ddL11 <- DD.L1(x,y,params[1],params[2],params[3],params[4])
  ddL22 <- DD.L2(x,y,params[1],params[2],params[3],params[4])
  ddL33 <- DD.L3(x,y,params[1],params[2],params[3],params[4])
  ddRho <- DD.Rho(x,y,params[1],params[2],params[3],params[4])
  # Mixed partial derivative
  ddL12 <- DD.L12(x,y,params[1],params[2],params[3],params[4])
  ddL13 <- DD.L13(x,y,params[1],params[2],params[3],params[4])
  ddL23 <- DD.L23(x,y,params[1],params[2],params[3],params[4])
  ddL1R <- DD.L1R(x,y,params[1],params[2],params[3],params[4])
  ddL2R <- DD.L2R(x,y,params[1],params[2],params[3],params[4])
  ddL3R <- DD.L3R(x,y,params[1],params[2],params[3],params[4])
  # Calculate the hessian element-wise
  df11 <- ddL11 * pmf_reciprocated - pmf_reciprocated_squared * dL1^2
  df12 <- ddL12 * pmf_reciprocated - pmf_reciprocated_squared * dL1*dL2
  df13 <- ddL13 * pmf_reciprocated - pmf_reciprocated_squared * dL1*dL3
  df1r <- ddL1R * pmf_reciprocated - pmf_reciprocated_squared * dL1*dRho
  df22 <- ddL22 * pmf_reciprocated - pmf_reciprocated_squared * dL2^2
  df23 <- ddL23 * pmf_reciprocated - pmf_reciprocated_squared * dL2*dL3
  df2r <- ddL2R * pmf_reciprocated - pmf_reciprocated_squared * dL2*dRho
  df33 <- ddL33 * pmf_reciprocated - pmf_reciprocated_squared * dL3^2
  df3r <- ddL3R * pmf_reciprocated - pmf_reciprocated_squared * dL3*dRho
  dfrr <- ddRho * pmf_reciprocated - pmf_reciprocated_squared * dRho^2
  # 'data' is cross-tabulated so everything aligns perfectly
  # in the element-wise product. What 'data' does is weighs the derivatives
  # according to the number of observations for the particular (N_1,N_2) value
  D11 <- sum(data*df11)
  D12 <- sum(data*df12)
  D13 <- sum(data*df13)
  D1r <- sum(data*df1r)
  D22 <- sum(data*df22)
  D23 <- sum(data*df23)
  D2r <- sum(data*df2r)
  D33 <- sum(data*df33)
  D3r <- sum(data*df3r)
  Drr <- sum(data*dfrr)
  # Populate the hessian matrix
  output[1,] <- c(D11,D12,D13,D1r)
  output[2,] <- c(D12,D22,D23,D2r)
  output[3,] <- c(D13,D23,D33,D3r)
  output[4,] <- c(D1r,D2r,D3r,Drr)
  return(output)
}
