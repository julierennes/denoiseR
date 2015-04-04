#' Singular Values Soft Thresholding
#'
#' This function estimates a low-rank signal from Gaussian noisy data using Soft Thresholding of the singular values. 
#' More precisely, any singular value smaller than a quantity lambda is set to zero (dl  = dl * max(1-lambda,0). This estimator minimizes the least-squares (Frobenius norm) penalized by the nuclear norm (sum of the singular values).  
#' The tuning parameter lambda is selected by minimizing an unbiased estimate of the risk namely a Stein Unbiased Risk Estimate (SURE). This requires to know the variance of the noise sigma.
#'  
#' @param X a data frame or a matrix with numeric entries
#' @param sigma integer, standard deviation of the Gaussian noise. By default sigma is estimated using the estim_sigma function with the MAD option 
#' @param center boolean, to center the data. By default "TRUE"
#' @param lambda0 integer, the initial value for lambda used to optimize SURE, by default the median of the singular values (must be in log scale)
#' @return mu.hat the estimator of the signal
#' @return nb.eigen the number of non-zero singular values
#' @return singval the singular values of the estimator
#' @return lambda the optimal lambda selected by minimizing SURE
#' @return low.rank the results of the SVD of the estimator
#' @details The SURE function is optimized in lambda using the optim function of the package stats (?optim). The initial lambda can be modified in the argument lambda0.  When sigma is not known, it can be estimated using the function estim_sigma. SVST automatically estimates the rank of the signal. Its value is given in the output nb.eigen corresponding to the number of non-zero eigenvalues. 
#' The estimated low rank matrix is given in the output mu.hat.
#' @references Candès, E. J., Sing-Long C. A. and Trzasko, J. D (2012). Unbiased risk estimates for singular value thresholding and spectral estimators. IEEE Transactions on Signal Processing 61(19), 4643-4657
#' @seealso \code{\link{estim_sigma}}
#' @seealso \code{\link{LRsim}}
#' @examples 
#' Xsim <- LRsim(200, 500, 100, 0.5)
#' truesigma <- 1/(0.5*sqrt(200*500))
#' svst.sure <-  SVST(Xsim$X, sigma = truesigma)
#' svst.sure$nb.eigen
#' svst.sure$lambda
#' svst.sure$singval
#' 
#'  Xsim <- LRsim(200, 500, 10, 4)
#'  sig <- estim_sigma(Xsim$X)
#'  svst.sure <- SVST(Xsim$X, sigma=sig)


SVST <- function(X,
                sigma,
                center = "TRUE",
                lambda0 = NA){

    ## Definition of the SURE(lambda) cost function
    SURE.optimisation <- function(X, sigma, lambda, svdXd){

      lambda <- exp(lambda)
      n <- nrow(X)
      p <- ncol(X)
      DD  <- svdXd
      DD2 <- DD^2
      lDD <- length(DD)
      D <- matrix(0, lDD, lDD)
      dhat <- DD * pmax( 1-(lambda/DD) , 0) 
      temp <- DD * dhat
      for (i in 1:lDD){
        DD2i <- DD2[i]
        diff2i <- DD2[i]-DD2
        diff2i[i] <- Inf
        D[i, ] <- temp[i]/diff2i
      }
      gradd <- (svdXd >= lambda)
      DIV <- sum(gradd + abs(n-p)*dhat/svdXd) + 2*sum(D)
      MSE.SVST <- -n*p*sigma^2 + sum((dhat-svdXd)^2)  + 2*sigma^2*DIV

      return(MSE.SVST)
}
    
    # housekeeping
    
    if(class(X) == "data.frame"){
      X <- as.matrix(X)
    }
    
    if(sum(sapply(X, is.numeric)) < ncol(X)){
      stop("all the variables are not numeric")
    }
    
    if(!is.na(sigma) & (sigma <= 0)){
      stop("sigma must be positive")
    }
   
    ###############
    ## Begin
        
    if(center == "TRUE"){
      moy <- apply(X, 2, mean)
      X <- scale(X, scale=F)
    }
    
    svdX <- svd(X) 
    svdXd <- svdX$d
    
    # infer unspecified choices
    
    if(is.na(sigma)){   
      sigma <- estim_sigma(X, method = "MAD", center = center)
      print(paste("sigma = ", round(sigma,6)))   
    }
    
    if(is.na(lambda0)){
      lambda0 <- log(median(svdXd))
    }  
    
    # find lambda minimizing SURE
    
    lambda.o <- optim(lambda0, SURE.optimisation, X = X, sigma = sigma, svdXd = svdXd)$par
    lambda.o <- exp(lambda.o)

    dhat <- svdXd * pmax( 1-(lambda.o/svdXd), 0) 
    temp <- t(t(svdX$u)*dhat)
    X.SVST <- temp %*% t(svdX$v)
    nbeig <- sum(svdXd > lambda.o)
    
    if(center == "TRUE"){ 
      X.SVST <-  X.SVST + matrix(moy, nrow(X), ncol(X), byrow = T)
    }
    
    svst.svd <- svdX
    svst.svd$d <- dhat   
    row.names(X.SVST) <- row.names(X)
    colnames(X.SVST) <- colnames(X)

    return(list(mu.hat = X.SVST, nb.eigen = nbeig, singval = dhat, lambda = lambda.o,  low.rank = svst.svd))
 }
