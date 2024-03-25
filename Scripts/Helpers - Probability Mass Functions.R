#' Bivariate Polya-Aeppli Probability Mass Function
#' Author(s): Ye Yang, Dean Hansen
#'
#' If you want to calculate f(N_1,N_2) it should enter the
#' function as pmf.M(N_2,N_1,L1,L2,L3,rhoo)
#'
#' @param x 
#' @param y 
#' @param L1 
#' @param L2 
#' @param L3 
#' @param rhoo 
#'
#' @return (y+1) by (x+1) grid where pmf.M(x,y) -> A[y+1,x+1] -> f(y,x)
#' where the returned value is the opposite of (x,y) it is (y,x). If you 
#' want to calculate the PMF f(i,j) at (x=5,y=3), you need to input it
#' in the reverse as pmf.M(x=3,y=5) -> A[6,4] -> f(5,3)
#'
#' @examples
#' x<-5; y<-3; L1<-2; L2<-4; L3<-3; rhoo<-0.25;
#' pmf.M(x,y,L1,L2,L3,rhoo)
pmf.M <- function(x,y,L1,L2,L3,rhoo) {
  
  a <- function(i) {2*rhoo + ((1-rhoo)*L1-2*rhoo)/i}
  b <- function(i) {rhoo^2*(1-(2/i))}
  c <- function(j) {2*rhoo + ((1-rhoo)*L2-2*rhoo)/j}
  v <- function(i) {2*rhoo^2-((1-rhoo)*(L3-rhoo*(L1+L3))+2*rhoo^2)/i}
  w <- function(j) {2*rhoo^2-((1-rhoo)*(L3-rhoo*(L2+L3))+2*rhoo^2)/j}
  
  # f(i,j) == A[j+1,i+1]
  A <- matrix(NA, nrow=(y+1), ncol=(x+1))
  
  # f(0,0)
  counter<-0
  i<-counter; im<-i+1
  j<-counter; jm<-j+1
  A[im,jm] <- exp(-(L1+L2+L3))
  if (x==0 && y==0) {return(A)}
  
  # f(0,1)
  j<-counter+1; jm<-j+1
  if (x != 0) {
    A[im,jm] <- c(j)*A[im,jm-1]
  }
  # f(0,j)
  # j=2,3,...;
  for (l in ((j+1):x)) {
    if (x <= 1) {break}
    jm <- l+1
    A[im,jm] <- c(l)*A[im,jm-1] - b(l)*A[im,jm-2]
  }
  
  # f(1,0)
  i<-counter+1; im<-i+1
  j<-counter; jm<-j+1
  if (y!=0) {
    A[im,jm] <- a(i)*A[im-1,jm]
  }
  # f(i,0)
  # i=2,3,...;
  for (k in ((i+1):y)) {
    if (y<=1) {break}
    im <- k+1
    A[im,jm] <- a(k)*A[im-1,jm] - b(k)*A[im-2,jm]
  }
  
  # if x=1 or y=1 stop
  if (min(x,y) <= counter) {return(A)}
  counter <- 1
  
  # f(i,j)
  # i=1;
  # j=1,2,...,x;
  i<-counter; im<-i+1
  j<-counter; jm<-j+1
  for (l in (j:x)) {
    jm <- l+1
    A[im,jm] <- rhoo*A[im,jm-1] + a(i)*A[im-1,jm] - v(i)*A[im-1,jm-1]
  }
  
  # f(i,j)
  # i=1,2,3,...,y;
  # j=1;
  j<-counter; jm<-j+1
  for (k in ((i+1):y)) {
    if (y < 2) {break}
    im <- k+1
    A[im,jm] <- rhoo*A[im-1,jm] + c(j)*A[im,jm-1] - w(j)*A[im-1,jm-1]
  }
  #recursion loop
  counter <- 2
  while (counter <= min(x,y)) {
    # f(i,j)
    # i=2,3,...;
    # j=1,2,...;
    i<-counter; im<-i+1
    j<-counter; jm<-j+1
    for (l in (j:x)) {
      jm <- l+1
      A[im,jm] <- rhoo*A[im,jm-1] + a(i)*A[im-1,jm] - b(i)*(A[im-2,jm] - rhoo*A[im-2,jm-1]) - v(i)*A[im-1,jm-1]
    }
    if ((y <= x) && (counter==y)) {break;}
    # f(i,j)
    # i=1,2,...;
    # j=2,3,...;
    j<-counter; jm<-j+1
    for (k in ((i+1):y)) {
      im <- k+1
      A[im,jm] <- rhoo*A[im-1,jm] + c(j)*A[im,jm-1] - b(j)*(A[im,jm-2] - rhoo*A[im-1,jm-2]) - w(j)*A[im-1,jm-1]
    }
    # move down one row and right one column
    counter <- counter + 1
  }
  return(A) 
}

#' Bivariate Negative Binomial Probability Mass Function
#' Author(s): Dean Hansen
#'
#' @param x1 
#' @param x2 
#' @param a0 
#' @param a1 
#' @param a2 
#' @param b1 
#' @param b2 
#'
#' @return PMF at (x1,x2)
#'
#' @examples
#' dbvnb(1,4,2,3,4,0.25,0.5)
#' dbvnb(10,0,3,4,5,0.2,0.5)
dbvnb <- function(x1, x2, a0, a1, a2, b1, b2) {
  output <- matrix(NA,nrow=(x1+1),ncol=(x2+1))
  for (k in 0:x1) {
    for (m in 0:x2) {
      output[k+1,m+1] <- choose(a0+x1+x2-k-m-1,a0+x2-m-1) * choose(a0+x2-m-1,a0-1) * choose(a1+k-1,a1-1) * choose(a2+m-1,a2-1) * ((b1^x1)*(b2^x2)*(b1+b2+1)^(k+m-x1-x2-a0)/((b1+1)^(k+a1) * (b2+1)^(m+a2)))
      }
    }
  output <- sum(output)
  return(output)
}
