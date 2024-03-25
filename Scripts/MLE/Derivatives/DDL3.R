#' Second derivative of PMF with respect to $\lambda_3$
#' Author(s): Dean Hansen
#'
#' If you want to calculate f''(N_1,N_2) it should enter the
#' function as DD.L3(N_2,N_1,L1,L2,L3,rhoo)#'
#' @param x 
#' @param y 
#' @param L1 
#' @param L2 
#' @param L3 
#' @param rhoo 
#'
#' @return
#' @export
#'
#' @examples
DD.L3 <- function(x,y,L1,L2,L3,rhoo) {
  
  # helper functions
  a <- function(i) {2*rhoo+((1-rhoo)*L1-2*rhoo)/i}
  b <- function(i) {rhoo^2*(1-(2/i))}
  c <- function(j) {2*rhoo+((1-rhoo)*L2-2*rhoo)/j}
  v <- function(i) {2*rhoo^2-((1-rhoo)*(L3-rhoo*(L1+L3))+2*rhoo^2)/i}
  w <- function(j) {2*rhoo^2-((1-rhoo)*(L3-rhoo*(L2+L3))+2*rhoo^2)/j}
  
  # derivatives of helper functions w.r.t $\lambda_3$
  a_prime <- function(i) {0}
  c_prime <- function(j) {0}
  b_prime <- function(i) {0}
  v_prime <- function(i) {((-1)*(1-rhoo)^2)/i}
  w_prime <- function(j) {((-1)*(rhoo-1)^2)/j}
  
  # derivatives of helper functions w.r.t $\lambda_3$
  a_prime_prime <- function(i) {0}
  c_prime_prime <- function(j) {0}
  b_prime_prime <- function(i) {0}
  v_prime_prime <- function(i) {0}
  w_prime_prime <- function(j) {0}
  
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
    A[im,jm] <- rhoo*A[im,jm-1] + a(i)*A[im-1,jm] - 2*v_prime(i)*D.L3(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v(i)*A[im-1,jm-1]
  }
  
  # f(i,j)
  # i=1,2,3,...,y;
  # j=1;
  j<-counter; jm<-j+1
  for (k in ((i+1):y)) {
    if (y < 2) {break}
    im <- k+1
    A[im,jm] <- rhoo*A[im-1,jm] + c(j)*A[im,jm-1] - 2*w_prime(j)*D.L3(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w(j)*A[im-1,jm-1]
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
      A[im,jm] <- rhoo*A[im,jm-1] + a(i)*A[im-1,jm] - b(i)*(A[im-2,jm] - rhoo*A[im-2,jm-1]) - 2*v_prime(i)*D.L3(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v(i)*A[im-1,jm-1]
    }
    if ((y <= x) && (counter==y)) {break;}
    # f(i,j)
    # i=1,2,...;
    # j=2,3,...;
    j<-counter; jm<-j+1
    for (k in ((i+1):y)) {
      im <- k+1
      A[im,jm] <- rhoo*A[im-1,jm] + c(j)*A[im,jm-1] - b(j)*(A[im,jm-2] - rhoo*A[im-1,jm-2]) - 2*w_prime(j)*D.L3(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w(j)*A[im-1,jm-1]
    }
    # move down one row and right one column
    counter <- counter + 1
  }
  return(A) 
}
