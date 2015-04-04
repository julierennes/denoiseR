#' Adaptive Shrinkage
#'
#' This function estimates a low-rank signal from Gaussian noisy data using the Adaptive Shrinker of the singular values. More precisely, the singular values are transformed using a function indexed by two parameters lambda and gamma as
#' dl  = dl * max(1-(lambda/dl)^gamma,0). This estimator is very flexible and adapts to the data whatever the noise regime. It can be seen as a compromise between hard and soft thresholding.  The parameters lambda and gamma are estimated by minimizing a Stein unbiased risk estimate (SURE) when the variance sigma^2 of the noise is known or a generalized SURE (GSURE) otherwise.
#' A method using an universal threshold for lambda is also available. 
#' 
#' @param X a data frame or a matrix with numeric entries
#' @param sigma integer, standard deviation of the Gaussian noise. By default sigma is estimated using the estim_sigma function with the MAD option 
#' @param method to select the two tunning parameters lambda and gamma. By default by minimizing GSURE
#' @param gamma.seq a vector for the sequence of gamma. (not used when method is UNIV). The values must be greater than 1
#' @param method.optim the method used in the optim function. By default BFGS
#' @param nbsim integer, number of replications used to calculate the universal threshold lambda when method is UNIV  
#' @param center boolean, to center the data. By default "TRUE"
#' @param lambda0 integer, the initial value for lambda used to optimize SURE and GSURE. By default the median of the singular values (must be in log scale)
#' @return mu.hat the estimator of the signal
#' @return nb.eigen the number of non-zero singular values
#' @return gamma the optimal gamma selected by minimizing SURE or GSURE
#' @return lambda the optimal lambda selected by minimizing SURE or GSURE
#' @return singval the singular values of the estimator
#' @return low.rank the results of the SVD of the estimator
#' @details When sigma is known, lambda and gamma can be estimated by minimizing SURE. To do this, a grid for gamma is defined in gamma.seq (gammas must be greater than 1) and the SURE function is optimized on lambda using the optim function of the package stats (?optim) with the optimization method by default sets to "BFGS". The initial lambda can be modified in the argument lambda0.  
#' A value for sigma has to be provided. When sigma is not known, it can be estimated using the function estim_sigma. An alternative which does not require to know or estimate sigma is estimate the two tuning parameters by minimizing GSURE.
#' UNIV consists in generating nbsim matrices of size n * p of Gaussian random variables with mean 0 and variance sigma^2 and computing the first singular value on each matrix. Then, the universal threshold lambda is calculated as the 1-alpha quantile of the null distribution (alpha is here sqrt(log(max(n,p)))). Then, gamma is estimated by minimizing a 1-dim SURE. This method is recommended when one is particularly interested in estimating the rank of the signal.        
#' The estimated low rank matrix is given in the output mu.hat. adashrink automatically estimates the rank of the signal. Its value is given in the output nb.eigen corresponding to the number of non-zero eigenvalues. 
#' @references Josse, J. & Sardy, S. (2015). Adaptive shrinkage of singular values. Statistics and Computing. 
#' @seealso \code{\link{estim_sigma}}
#' @seealso \code{\link{LRsim}}
#' @examples 
#' Xsim <- LRsim(200, 500, 100, 1)
#' ada.gsure <- adashrink(Xsim$X, method = "GSURE")
#' ada.gsure$nb.eigen
#' ada.gsure$singval
#' ada.gsure$lambda
#' ada.gsure$gamma
#' 
#' Xsim <- LRsim(200, 500, 10, 4)
#' sig <- estim_sigma(Xsim$X)
#' ada.sure <- adashrink(Xsim$X, method = "SURE", sigma = sig)

adashrink <- function(X,
                     sigma = NA,
                     method = c("GSURE", "UNIV", "SURE"), 
                     gamma.seq = seq(1, 5, by=.1),
                     nbsim = 500,
                     method.optim = "BFGS",
                     center = "TRUE",
                     lambda0 = NA){
  

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
  
  if(sum(gamma.seq<1)>0){
    stop("the gammas may be greater than 1")
  }
 
     
 method <- match.arg(method, c("GSURE", "Gsure", "gsure", "gSure","GSure",
                              "UNIV", "Univ", "univ", "SURE", "sure", "Sure"), several.ok = T)[1]
 method <- tolower(method)
  
  
  ## Definition of the SURE(lambda, gamma) cost function
  SURE.cost <- function(X,
                        method.ada,
                        sigma,
                        gamma,
                        lambda,
                        svdXd){
    n <- nrow(X)
    p <- ncol(X)
    if(method.ada == "univ"){
      gamma <- exp(gamma) + 1
    } else { 
    lambda <- exp(lambda)
    }  
    
    DD  <- svdXd
    DD2 <- DD^2
    lDD <- length(DD)
    D <-  matrix(0, lDD, lDD)
    dhat <- DD * pmax(1-(lambda/DD)^gamma, 0) 
    temp <- DD * dhat
    for (i in 1:lDD){
        DD2i <- DD2[i]
        diff2i <- DD2[i]-DD2
        diff2i[i] <- Inf
        D[i, ] <- temp[i]/diff2i
    }

    gradd <- (1+(gamma-1)*(lambda/svdXd)^gamma) * (svdXd >= lambda)
    DIV <- sum(gradd + abs(n-p)*dhat/svdXd) + 2*sum(D)
	
    if(method.ada == "gsure") {
      mse.ada <- sum((dhat-svdXd)^2)/(1-DIV/n/p)^2
    } else{       
     mse.ada <- -n*p*sigma^2 + sum((dhat-svdXd)^2)  + 2*sigma^2 *DIV
    }
    return(mse.ada)
}
	
  ###############
  ## Begin
 
 if(center == "TRUE"){
   moy <- apply(X, 2, mean)
   X <- scale(X, scale=F)
  }
 
 # infer unspecified choices
 
 if(is.na(sigma)&(method != "gsure")){   
   sigma <- estim_sigma(X, method = "MAD", center = center)
   warning(paste("sigma estimated by MAD:", round(sigma,6)))
 }
 
  val.optglob <- Inf
  svdX <- svd(X)
  n <- nrow(X)
  p <- ncol(X)
  svdXd <- svdX$d	
  
  if(method == "univ"){ 
   nbsim <- nbsim
   maxd <- rep(1, nbsim)
   for(i in 1:nbsim){
     maxd[i] <- max(svd(matrix(rnorm(n*p, sd = sigma), n, p))$d)
   }
   lambda.o <- quantile(maxd, 1-1/sqrt(log(max(n, p)))) 
   res.opti <- optim(0, SURE.cost, X = X, sigma = sigma, lambda = lambda.o, svdXd = svdXd, method.ada = method, method = method.optim)
   gamma.o <- exp(res.opti$par)+1
  } else {
   gamma.seq <- gamma.seq
   if(is.na(lambda0)){
     lambda0 <- log(median(svdXd))
   }  
   for (gamma in gamma.seq){
    res.opti <- optim(lambda0, SURE.cost, X = X, sigma = sigma, gamma = gamma, svdXd = svdXd,  method.ada = method, method = method.optim)
    lambda.temp <- exp(res.opti$par)
    val.opt <- res.opti$val
     if (val.opt < val.optglob) {	
       gamma.o <- gamma
       lambda.o <- lambda.temp
       val.optglob <- val.opt
    }
   }
  }
  
  dhat <- svdXd * pmax( 1-(lambda.o/svdXd)^gamma.o, 0) 
  temp <- t(t(svdX$u)*dhat)
  X.adaptive <- temp %*% t(svdX$v)
  adashrink.svd <- svdX
  adashrink.svd$d <- dhat
  nbeig <- sum(svdXd > lambda.o)
  if(center == "TRUE"){
    X.adaptive <-  X.adaptive + matrix(moy, nrow(X), ncol(X), byrow=T)
  }
  row.names(X.adaptive) <- row.names(X)
  colnames(X.adaptive) <- colnames(X)

  return(list(mu.hat = X.adaptive,
              nb.eigen = nbeig,
              gamma = gamma.o,
              lambda = lambda.o, 
              singval = dhat,
              low.rank = adashrink.svd))
}
