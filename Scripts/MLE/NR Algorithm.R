#' Newton-Rapshon Algorithm
#' Author(s): Ye-Yang, Dean Hansen
#'
#' @param params 
#' @param data cross-tabulated dataset
#'
#' @return
#' @export
#'
#' @examples
#' NR.F(params=c(1,2,3,0.25),data)
NR.F <- function(params, data, show=TRUE) {
  m <- length(params)
  max <- 10
  err <- 1e-4
  ii <- 0
  er <- rep(max,m)
  Pa <- params
  loop <- 0
  output <- matrix(NA, ncol=2, nrow=m)
  while ((ii <= max) && all(er >= err)) {
    Pb <- Pa - solve(hL(Pa,data)) %*% dL(Pa,data)
    e1 <- abs(Pb-Pa)
    Pa <- Pb
    if ((min(Pa) <= 0) || is.na(Pa[1])) {break;}
    er <- max(e1)
    ii <- ii+1
    output[,1] <- Pa
    output[,2] <- e1
    if (show) {print(loop)}
    loop <- loop+1
  }
  output <- data.frame(output,row.names=c("L1","L2","L3","Rho"))
  names(output) <- c("MLE", "Abs. Err")
  return(output)
}
