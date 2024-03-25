#' Method of Moments for Bivariate Polya-Aeppli Distribution
#' Author(s): Ye-Yang, Dean Hansen
#'
#' @param data cross-tabulated dataset
#'
#' @return
#' @export
#'
#' @examples
MofM <- function(data) {
  i1<-nrow(data)
  j1<-ncol(data)
  i2 <- 0
  j2 <- 0
  for(i in 1:i1) {
    i2 <- i2 + sum(data[i,])*(i-1)
  }
  for(j in 1:j1) {
    j2 <- j2 + sum(data[,j])*(j-1)
  }
  sa1 <- 0
  sa2 <- 0
  sb <- 0
  n1.b <- i2/sum(data)
  n2.b <- j2/sum(data)
  for (i in 1:i1) {
    sa1 <- sa1 + sum(data[i,])*((i-1)-n1.b)^2
  }
  for(j in 1:j1) {
    sa2 <- sa2+sum(data[,j])*((j-1)-n2.b)^2
  }
  for (i in 1:i1) {
    for (j in 1:j1) {
      sb <- sb+data[i,j]*((i-1)-n1.b)*((j-1)-n2.b)
    }
  }
  s.1 <- sa1/(sum(data)-1)
  s.2 <- sa2/(sum(data)-1)
  s.12 <- sb/(sum(data)-1)
  rho.h <- (s.1+s.2-n1.b-n2.b)/(s.1+s.2+n1.b+n2.b)
  L1.h <- (1-rho.h)*n1.b-(1-rho.h)^2*s.12
  L2.h <- (1-rho.h)*n2.b-(1-rho.h)^2*s.12
  L3.h <- (1-rho.h)^2*s.12
  rho.h <- (s.1+s.2-n1.b-n2.b)/(s.1+s.2+n1.b+n2.b)
  output <- c(L1=L1.h, L2=L2.h,L3=L3.h, rhoo=rho.h)
  return(output)
}

#' Method of Moments with Bootstrap CIs
#' Author(s): Ye-Yang, Dean Hansen
#'
#' @param data 
#' @param P 
#' @param B 
#'
#' @return
#' @export
#'
#' @examples
MofMF <- function(data,P,B) {
  MofM_estimates <- MofM(data)
  m <- length(P)
  if (min(MofM_estimates) <= 0) {return(NA)}
  b <- matrix(NA,ncol=m,nrow=B)
  N <- sum(data)
  for (i in 1:B) {
    Da <- BivPA(MofM_estimates,N)
    b[i,] <- MofM(Da)
  }
  s <- apply(b,2,sort)
  output <- apply(s,2,function(x){quantile(x,c(0.025,0.975))})
  output <- data.frame(t(output), row.names = c("L1","L2","L3","Rho"),check.names=F)
  return(output)
}
