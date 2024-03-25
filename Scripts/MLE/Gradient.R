#' Gradient for Bivariate Pólya-Aeppli
#' Author(s): Ye-Yang, Dean Hansen
#'
#' The dataset must be cross-tabulated for all values of
#' of c(0,max(x)) and c(0,max(y)). If there is no observation
#' present, the element-wise product will set the f'/f term 
#' to zero in the final gradient calculation so this will not
#' cause an issue in the calculation.
#' 
#' @param x number of columns in (cross-tabulated) dataset minus one 
#' @param y number of rows in (cross-tabulated) dataset minus one
#' @param params vector of parameters for BivPA density
#' @param data (y+1) by (x+1) cross-tabulated dataset 
#'
#' @return gradient vector evaluated over entire dataset
#'
#' @examples 
#' data <- bivpois::rbp(500,c(1,2,4))
#' data <- crosstabulate(data)
#' params <- c(1,2,4,0.34)
#' dL(params, data)
dL <- function(params, data) {
  x <- ncol(data) - 1 
  y <- nrow(data) - 1
  pmf_reciprocated <- 1 / pmf.M(x,y,params[1],params[2],params[3],params[4])
  f1f <- D.L1(x,y,params[1],params[2],params[3],params[4]) * pmf_reciprocated
  f2f <- D.L2(x,y,params[1],params[2],params[3],params[4]) * pmf_reciprocated
  f3f <- D.L3(x,y,params[1],params[2],params[3],params[4]) * pmf_reciprocated
  frf <- D.Rho(x,y,params[1],params[2],params[3],params[4]) * pmf_reciprocated
  dL1 <- sum(f1f*data)
  dL2 <- sum(f2f*data)
  dL3 <- sum(f3f*data)
  dRho <- sum(frf*data)
  output <- c(dL1,dL2,dL3,dRho)
  return(output)
}
