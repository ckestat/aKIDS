# Gaussian kernel Gram matrix
Gram_RBF <- function(x, sigma2=NULL, is.dist=F){
  if(!is.dist){distx = dist(x, diag=T, upper=T)^2}else{distx = x}
  if(is.null(sigma2)){sigma2 = 0.5*median(distx)}
  return(exp(-as.matrix(distx)/2/sigma2))
}

# Nadaraya-Watson Estimation
NW_RBF <- function(y, n, bw=NULL){
  if(is.null(bw)){bw = 1.06*sd(y)*n^(-1/5)}
  dist = dist(y,diag=T,upper=T)^2
  G = exp(-as.matrix(dist)/2/bw^2)
  Gstar = G/rowSums(G)
  tGG = t(Gstar)%*%Gstar
  return(tGG)
}

# generate multivariate data
rnorm_Choleski <- function(n, mu, covMatrix) {
  
  p <- dim(covMatrix)[1]
  R <- chol(covMatrix)
  matrix(rnorm(p * n), n, p) %*% R + matrix(mu, n, p, byrow = TRUE)
  
}