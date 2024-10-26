#' Kernel-based Independence Dual Screening (KIDS)
#'
#' This function performs feature screening on ultrahigh dimensional right-censored data 
#' using the Kernel-based Independence Dual Screening (KIDS) procedure. 
#' The KIDS method selects important features by evaluating two kernel-based R-squared measures:
#' (1) a marginal measure: kernel-based R-squared between the feature and the censoring status, and
#' (2) a conditional measure: kernel-based partial R-squared between the feature and the observed survival time, adjusted for the censoring status.
#'
#' @param x A data matrix where each row corresponds to an observation and each column corresponds to a feature.
#' @param y A numeric vector representing the observed survival times.
#' @param delta A binary vector of the same length as `y`, indicating the censoring status (1 = event occurred, 0 = censored).
#' @param d An integer representing the number of features to be selected. If not supplied (NULL), the function returns utility measures for all features instead of selected features.
#' @param ir A logical value indicating whether to apply inverse regression for estimating the conditional measure. Set to `TRUE` to enable inverse regression, or `FALSE` to disable it.
#' 
#' @return If `d` is provided, the function returns a vector of indices corresponding to the `d` selected features. 
#' If `d` is not provided, the function returns utility measures for each feature, which can be used for feature ranking or selection.
#' 
#' @examples 
#' ### Simulate data from AFT model
#' n <- 200
#' p <- 5000
#' rho <- 0.5
#' Sigma <- rho^as.matrix(dist(1:p,diag = T,upper = T))
#' x <- rnorm_Choleski(n, rep(0,p), Sigma)
#' t <- exp(x[,1]+x[,2]+1.5*x[,7]^2+rnorm(n))
#' c <- rexp(n,1/5/exp(x[,1]))
#' y <- pmin(t,c)
#' delta <- as.numeric(t<c)
#' 
#' ### perform KIDS
#' ix = KIDS(x, y, delta, d = floor(n/log(n)), ir = F)
#' ix
#' 
#' @export
#' 
KIDS <- function(x, y, delta, d, ir){
  
  n <- nrow(x)
  p <- ncol(x)
  id0 <- delta==0
  n0 <- sum(id0)
  n1 <- n-n0
  w0 <- n0/n
  w1 <- 1-w0
  ohat <- matrix(0,p,2)
  
  if(ir==T){
    
    x0 <- x[id0,]
    x1 <- x[!id0,]
    bwy <- 1.06*sd(y)*n^(-1/5)
    tGGy0 <- NW_RBF(y[id0], n0, bwy)
    tGGy1 <- NW_RBF(y[!id0], n1, bwy)
    
    for(j in 1:p){
      
      distu = dist(x[,j], diag=T, upper=T)^2
      sigma2u = 0.5*median(distu)
      if(sigma2u==0){sigma2u = 0.001}
      K = Gram_RBF(distu, sigma2u, T)
      K0 = Gram_RBF(x0[,j], sigma2u)
      K1 = Gram_RBF(x1[,j], sigma2u)      
      mK = mean(K)
      mK0 = mean(K0)
      mK1 = mean(K1)
      wmK = mK0*w0 + mK1*w1
      ohat[j,1] = (wmK-mK)/(1-mK)
      
      Hy0 = sum(tGGy0*K0)/n0 - mK0
      Hy1 = sum(tGGy1*K1)/n1 - mK1
      ohat[j,2] = (Hy0*w0 + Hy1*w1)/(1-wmK)
      
    }
    
  }else{
    
    x0 <- x[id0,]
    x1 <- x[!id0,]
    disty = dist(y, diag=T, upper=T)^2
    sigma2y = 0.5*median(disty)
    Ky0 <- Gram_RBF(y[id0], sigma2y)
    Ky1 <- Gram_RBF(y[!id0], sigma2y)
    mKy0 = mean(Ky0)
    mKy1 = mean(Ky1)
    wmKy = mKy0*w0 + mKy1*w1
    
    for(j in 1:p){
      
      distu = dist(x[,j], diag=T, upper=T)^2
      sigma2u = 0.5*median(distu)
      if(sigma2u==0){sigma2u = 0.001}
      K = Gram_RBF(distu, sigma2u, T)
      K0 = Gram_RBF(x0[,j], sigma2u)
      K1 = Gram_RBF(x1[,j], sigma2u)      
      mK = mean(K)
      mK0 = mean(K0)
      mK1 = mean(K1)
      wmK = mK0*w0 + mK1*w1
      ohat[j,1] = (wmK-mK)/(1-mK)
      
      bwu <- 1.06*sd(x[,j])*n^(-1/5)
      tGGx0 = NW_RBF(x0[,j], n0, bwu)
      tGGx1 = NW_RBF(x1[,j], n1, bwu)
      Hy0 = sum(tGGx0*Ky0)/n0 - mKy0
      Hy1 = sum(tGGx1*Ky1)/n1 - mKy1
      ohat[j,2] = (Hy0*w0 + Hy1*w1)/(1-wmKy)
      
    }
    
  }
  
  if(is.null(d)){return(ohat)}else{
    r1 <- rank(-ohat[,1])
    r2 <- rank(-ohat[,2])
    ow <- order(pmin(r1,r2),pmax(r1,r2))
    return(ow[1:d])
  }
  
}