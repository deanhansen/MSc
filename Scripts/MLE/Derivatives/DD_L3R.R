#' Mixed derivatives w.r.t $\lambda_3$ and $\rho$
#' Author(s): Dean Hansen
#'
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
DD.L3R <- function(x,y,L1,L2,L3,rhoo) {
  
  # helper functions
  a <- function(i) {2*rhoo+((1-rhoo)*L1-2*rhoo)/i}
  b <- function(i) {rhoo^2*(1-(2/i))}
  c <- function(j) {2*rhoo+((1-rhoo)*L2-2*rhoo)/j}
  v <- function(i) {2*rhoo^2-((1-rhoo)*(L3-rhoo*(L1+L3))+2*rhoo^2)/i}
  w <- function(j) {2*rhoo^2-((1-rhoo)*(L3-rhoo*(L2+L3))+2*rhoo^2)/j}
  
  # derivatives of helper functions w.r.t $\lambda_3$
  a_prime_3 <- function(i) {0}
  c_prime_3 <- function(j) {0}
  b_prime_3 <- function(i) {0}
  v_prime_3 <- function(i) {((-1)*(1-rhoo)^2)/i}
  w_prime_3 <- function(j) {((-1)*(rhoo-1)^2)/j}
  
  # derivatives of helper functions w.r.t $\rho$
  a_prime_r <- function(i) {2 + (-L1-2)/i}
  c_prime_r <- function(j) {2 + (-L2-2)/j}
  b_prime_r <- function(i) {2*rhoo*(1-(2/i))}
  v_prime_r <- function(i) {4*rhoo - (L1*(2*rhoo-1)+2*L3*(rhoo-1)+4*rhoo)/i}
  w_prime_r <- function(j) {4*rhoo - (L2*(2*rhoo-1)+2*L3*(rhoo-1)+4*rhoo)/j}
  
  # derivatives of helper functions w.r.t $\lambda_3$ and $\rho$
  a_prime_3_r <- function(i) {0}
  c_prime_3_r <- function(j) {0}
  b_prime_3_r <- function(i) {0}
  v_prime_3_r <- function(i) {(2*(1-rhoo))/i}
  w_prime_3_r <- function(j) {(2*(1-rhoo))/j}
  
  # f(i,j) == A[j+1,i+1]
  A <- matrix(NA, nrow=(y+1), ncol=(x+1))
  
  # f(0,0)
  counter<-0
  i<-counter; im<-i+1
  j<-counter; jm<-j+1
  A[im,jm] <- 0
  if (x==0 && y==0) {return(A)}
  
  # f(0,1)
  j<-counter+1; jm<-j+1
  if (x != 0) {
    A[im,jm] <- c_prime_r(j)*D.L3(j-1,i,L1,L2,L3,rhoo)[i+1,j] + c(j)*A[im,jm-1]
  }
  # f(0,j)
  # j=2,3,...;
  for (l in ((j+1):x)) {
    if (x <= 1) {break}
    jm <- l+1
    # l=j
    A[im,jm] <- c_prime_r(l)*D.L3(l-1,i,L1,L2,L3,rhoo)[i+1,l] + c(l)*A[im,jm-1] - b_prime_r(l)*D.L3(l-2,i,L1,L2,L3,rhoo)[i+1,l-1] - b(l)*A[im,jm-2]
  }
  
  # f(1,0)
  i<-counter+1; im<-i+1
  j<-counter; jm<-j+1
  if (y!=0) {
    A[im,jm] <- a_prime_r(i)*D.L3(j,i-1,L1,L2,L3,rhoo)[i,j+1] + a(i)*A[im-1,jm]
  }
  # f(i,0)
  # i=2,3,...;
  for (k in ((i+1):y)) {
    if (y<=1) {break}
    im <- k+1
    # k=i
    A[im,jm] <- a_prime_r(k)*D.L3(j,k-1,L1,L2,L3,rhoo)[k,j+1] + a(k)*A[im-1,jm] - b_prime_r(k)*D.L3(j,k-2,L1,L2,L3,rhoo)[k-1,j+1] - b(k)*A[im-2,jm]
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
    # l=j
    A[im,jm] <- D.L3(l-1,i,L1,L2,L3,rhoo)[i+1,l] + rhoo*A[im,jm-1] + a_prime_r(i)*D.L3(l,i-1,L1,L2,L3,rhoo)[i,l+1] + a(i)*A[im-1,jm] - v_prime_3_r(i)*pmf.M(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v_prime_r(i)*D.L3(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v_prime_3(i)*D.Rho(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v(i)*A[im-1,jm-1]
  }
  
  # f(i,j)
  # i=1,2,3,...,y;
  # j=1;
  j<-counter; jm<-j+1
  for (k in ((i+1):y)) {
    if (y < 2) {break}
    im <- k+1
    # k=i
    A[im,jm] <- D.L3(j,k-1,L1,L2,L3,rhoo)[k,j+1] + rhoo*A[im-1,jm] + c_prime_r(j)*D.L3(j-1,k,L1,L2,L3,rhoo)[k+1,j] + c(j)*A[im,jm-1] - w_prime_3_r(j)*pmf.M(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w_prime_r(j)*D.L3(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w_prime_3(j)*D.Rho(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w(j)*A[im-1,jm-1]
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
      A[im,jm] <- D.L3(l-1,i,L1,L2,L3,rhoo)[i+1,l] + rhoo*A[im,jm-1] + a_prime_r(i)*D.L3(l,i-1,L1,L2,L3,rhoo)[i,l+1] + a(i)*A[im-1,jm] - b_prime_r(i)*(D.L3(l,i-2,L1,L2,L3,rhoo)[i-1,l+1] - rhoo*D.L3(l-1,i-2,L1,L2,L3,rhoo)[i-1,l]) - b(i)*(A[im-2,jm] - D.L3(l,i-2,L1,L2,L3,rhoo)[i-1,l+1] - rhoo*A[im-2,jm-1]) - v_prime_3_r(i)*pmf.M(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v_prime_r(i)*D.L3(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v_prime_3(i)*D.Rho(l-1,i-1,L1,L2,L3,rhoo)[i,l] - v(i)*A[im-1,jm-1]
    }
    if ((y <= x) && (counter==y)) {break;}
    # f(i,j)
    # i=1,2,...;
    # j=2,3,...;
    j<-counter; jm<-j+1
    for (k in ((i+1):y)) {
      im <- k+1
      # k=i
      A[im,jm] <- D.L3(j,k-1,L1,L2,L3,rhoo)[k,j+1] + rhoo*A[im-1,jm] + c_prime_r(j)*D.L3(j-1,k,L1,L2,L3,rhoo)[k+1,j] + c(j)*A[im,jm-1] - b_prime_r(j)*(D.L3(j-2,k,L1,L2,L3,rhoo)[k+1,j-1] - rhoo*D.L3(j-2,k-1,L1,L2,L3,rhoo)[k,j-1]) - b(j)*(A[im,jm-2] - D.L3(j-2,k,L1,L2,L3,rhoo)[k+1,j-1] - rhoo*A[im-1,jm-2]) - w_prime_3_r(j)*pmf.M(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w_prime_r(j)*D.L3(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w_prime_3(j)*D.Rho(j-1,k-1,L1,L2,L3,rhoo)[k,j] - w(j)*A[im-1,jm-1]
    }
    # move down one row and right one column
    counter <- counter + 1
  }
  return(A) 
}
