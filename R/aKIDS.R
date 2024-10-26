#' alpha-controlled Kernel-based Independence Dual Screening (alpha-KIDS)
#'
#' This function performs feature screening with FDR control on ultrahigh dimensional right-censored data.
#' using the alpha-controlled Kernel-based Independence Dual Screening (alpha-KIDS) procedure
#' The alpha-KIDS method is a two-stage procedure that initially selects a set of potentially important features using the KIDS method,
#' followed by identifying a refined set under FDR control through a unified knockoﬀ procedure.
#'
#' @param x A data matrix where each row corresponds to an observation and each column corresponds to a feature.
#' @param y A numeric vector representing the observed survival times.
#' @param delta A binary vector of the same length as `y`, indicating the censoring status (1 = event occurred, 0 = censored).
#' @param n1.ratio A numeric value representing the proportion of observations to be used for KIDS screening. 
#' @param d An integer representing the number of features to be selected in the KIDS step. 
#' @param ir A logical value indicating whether to apply inverse regression for estimating the conditional measure. Set to `TRUE` to enable inverse regression, or `FALSE` to disable it.
#' @param modelx A logical value indicating the type of knockoff features to generate. If `TRUE`, the function generates approximate second-order modelX knockoff features; if `FALSE`, it generates fixed knockoff features.
#' @param offset A numeric value representing the offset to be applied in the FDR threshold. The default value is 1.
#' @param alpha A numeric vector specifying one or more false discovery rate (FDR) control levels.
#' 
#' @return The function returns:
#' \describe{
#'   \item{ix1}{A vector of indices corresponding to the features selected by the KIDS screening step.}
#'   \item{alpha}{The vector of FDR control levels that were applied during the selection process.}
#'   \item{ix}{A list of vectors, where each vector contains the indices of features selected by the alpha-KIDS procedure with FDR control for each value of `alpha`.}
#' }
#' 
#' @examples 
#' ### Simulate data from AFT model
#' n <- 2000
#' p <- 5000
#' rho <- 0.3
#' Sigma <- rho^as.matrix(dist(1:p,diag = T,upper = T))
#' x <- rnorm_Choleski(n, rep(0,p), Sigma)
#' t <- exp(x[,1]+x[,2]+1.5*x[,3]^2+x[,4]+x[,5]/log(1+abs(x[,5]))+x[,6]+1.5*sin(1/2/x[,7])+x[,8]+x[,9]+x[,10]+rnorm(n)/2)
#' c <- rexp(n,1/7/exp(x[,1]+x[,2]-x[,10]))
#' y <- pmin(t,c)
#' delta <- as.numeric(t<c)
#' ### Perform alpha-KIDS
#' ix <- aKIDS(x, y, delta, n1.ratio=0.25, d=100, ir=F, alpha=c(0.1,0.2))
#' ix
#' 
#' @export
#'
aKIDS <- function(x, y, delta, n1.ratio, d, ir, modelx = T, offset = 1, alpha = 0.1){
  
  ### KIDS
  ind1 <- 1:round(length(y)*n1.ratio)
  ix1 <- KIDS(x[ind1,], y[ind1], delta[ind1], d, ir)
  
  ### Knockoff
  # knockoff feature
  xk <- x[,ix1]
  if(modelx==F){
    xnorm <- apply(x[-ind1,ix1], 2, function(xj) sqrt(sum(xj^2)))
    sx2 <- sweep(x[-ind1,ix1], MARGIN = 2, STATS = xnorm, FUN = '/')
    sxk2 <- knockoff::create.fixed(sx2, randomize = T)$Xk
    xk[-ind1,] <- sweep(sxk2, MARGIN = 2, STATS = xnorm, FUN = '*')
  }else{
    xk[-ind1,] <- knockoff::create.second_order(x[-ind1,ix1], method="sdp")
  }
  # knockoff filters
  W <- KIDS(x[,ix1], y, delta, NULL, ir) - KIDS(xk, y, delta, NULL, ir)
  # thresholds and false discovery proportions
  ts <- rbind(rep(0,2),
              cbind(rep(sort(c(abs(W[,1]),Inf)),each=d+1),
                    rep(sort(c(abs(W[,2]),Inf)),d+1)))
  fdp <- sapply(1:nrow(ts), function(i) (offset + sum(W[,1] <= -ts[i,1] | W[,2] <= -ts[i,2])) / sum(W[,1] >= ts[i,1] | W[,2] >= ts[i,2]))
  ix <- vector("list", length(alpha))
  
  for(i in 1:length(alpha)){
    ok <- which(fdp <= alpha[i])
    if(length(ok)==0){
      ix[[i]] <- numeric(0)
    }else{
      # minimal elements
      ts.ok <- ts[ok,,drop=F]
      ts.ok1 <- do.call(rbind, by(ts.ok, list(ts.ok[,1]), head, n=1))
      ts.ok1s <- ts.ok1[order(ts.ok1[,2],ts.ok1[,1]),]
      ts.ok2 <- do.call(rbind, by(ts.ok1s, list(ts.ok1s[,2]), head, n=1))
      mes <- ts.ok2
      if(nrow(ts.ok2)>1){
        rm <- NULL
        for(j in 2:nrow(ts.ok2)){
          if(ts.ok2[j,1]>min(ts.ok2[1:(j-1),1])){rm = c(rm,j)}
        }
        if(!is.null(rm)){mes = mes[-rm,]}
      }
      row.names(mes) <- NULL
      # best minimal element that maximizes average utilities
      dt.tmp <- apply(mes, 1,
                      function(ti){
                        ix.tmp <- which(W[,1]>=ti[1]|W[,2]>=ti[2])
                        mean(W[ix.tmp,])
                      })
      th <- mes[which.max(dt.tmp),]
      ix.tmp <- unique(as.vector(ix1[which(W[,1]>=th$V1|W[,2]>=th$V2)]))
      if(i>1){
        if(length(ix.tmp)<length(ix[[i-1]])){ix[[i]]<-ix[[i-1]]}else{ix[[i]]<-ix.tmp}
      }else{
        ix[[i]]<-ix.tmp
      }
    }
  }
  return(list(ix1=ix1, alpha=alpha, ix=ix))
  
}
