
#########################################################
### Section 5.1 - SIMU E.candes
#########################################################
rm(list=ls())
library("denoiseR")
library("FactoMineR")

## Functions
## SA 
SA <- function (X, lambda, S){   
  n <- nrow(X) ; p <- ncol(X)
  moyX <- apply(X, 2,mean)
  Xc <- scale(X, scale = F)
  recon <-  as.matrix(scale(X,scale=F))  %*%   solve(t(Xc)%*%(Xc) + diag(n*(p/min(p, (n-1)))*lambda, ncol(Xc)))   %*%  t(Xc)%*%(Xc)   
  recon <- reconst(PCA(recon, graph = F,scale = F, ncp = S),S)
  recon <- sweep(recon, 2, moyX, FUN="+")
  row.names(recon) = row.names(X)
  colnames(recon) = colnames(X)  
  return(recon)
}

TSVD_opt <- function(X, sigma = NULL){     
  n <- nrow(X)
  p <- ncol(X)
  n <- nrow(X) ; p <- ncol(X)
  moyX <- apply(X, 2,mean)
  X <- scale(X, scale = F)
  
  svdX <- svd(X)
  beta <- min(n, p)/max(n, p)
  lambdastar <- sqrt( 2*(beta+1) + 8*beta/((beta+1+(sqrt(beta^2+14*beta+1))))  )
  if (is.null(sigma)){
    wbstar <-  0.56*beta^3- 0.95*beta^2+1.82*beta+1.43
    ind <- which((svdX$d - wbstar*median(svdX$d))>0)  
  }else{
    ind <- which(svdX$d - lambdastar*sqrt(max(n, p)* sigma) >0)
  }      
  if (length(ind)==0){
    X.TSVD_opt <- matrix(0, n, p)
    X.TSVD_opt <- sweep(X.TSVD_opt, 2, moyX, FUN="+")
    nbeigTSVD_opt <- 0
  }else{  
    X.TSVD_opt <- svdX$u[, 1:length(ind), drop = F]%*% diag(svdX$d[1:length(ind)],length(ind),length(ind)) %*% t(svdX$v[, 1:length(ind), drop = F])
    X.TSVD_opt <- sweep(X.TSVD_opt, 2, moyX, FUN="+")
     nbeigTSVD_opt  <- length(ind)   
  }  
  return(list(recon = X.TSVD_opt, nbeig = nbeigTSVD_opt))}

## Code for simulation
nbsim <-50
lossFrob <- function(X.loss) sum((X.loss - MU)^2)

n <- 200
p <- 500
sigma <- 1/(SNR*sqrt(n*p))
k <- 100
SNR <- 1

results.list <- lapply(1:nbsim, function(iter){
  
  Xsim <- LRsim(n = n, p = p, k = k, SNR = SNR)
  X <- Xsim$X
  MU <- Xsim$mu
  
  SA.muhat <- SA(X, lambda = sigma^2, S = k)
  SA.loss <- lossFrob(SA.muhat)
  
  TSVD.muhat <- reconst(PCA(X, graph = F, scale = F, ncp = k), k)
  TSVD.loss <- lossFrob(TSVD.muhat)
  
  TSVD_opt.res <- TSVD_opt(X, sigma = sigma^2)
  TSVD_opt.loss <- lossFrob(TSVD_opt.res$recon)
  
  ISA.res  <-  ISA(X, sigma = sigma, maxiter = 50)
  ISA.loss  <- lossFrob(ISA.res$mu.hat)

  ASYMPT.res  <- optishrink(X, method = "ASYMPT", sigma = sigma)
  ASYMPT.loss  <- lossFrob(ASYMPT.res$mu.hat)
 
  SVST.res  <- adashrink(X, method = "SURE", sigma = sigma, gamma.seq = 1)
  SVST.loss <- lossFrob(SVST.res$mu.hat)
  
  LN.res <- optishrink(X, method = "LN", k = k, sigma = sigma) 
  LN.loss <- lossFrob(LN.res$mu.hat)

  loss.res  <- rbind(SA=SA.loss, ISA=ISA.loss, TSVD=TSVD.loss,  TSVD_opt=TSVD_opt.loss, ASYMPT=ASYMPT.loss, SVST=SVST.loss, LN=LN.loss)
  nbeig <- rbind(SA=k, ISA=ISA.res$nb.eigen, TSVD=k, TSVD_opt= TSVD_opt.res$nbeig, ASYMPT=ASYMPT.res$nb.eigen, SVST=SVST.res$nb.eigen, LN=k)
  cbind.data.frame(loss.res, nbeig)
})
results <- Reduce("+", results.list) / length(results.list)

#########################################################
### Section 5.2 - Poisson Noise
#########################################################

rm(list=ls())
library("denoiseR")
library("FactoMineR")

## Functions
TSVD  <- function(X, K = 2) {
  X.svd  <- svd(X, nu = K, nv = K)
    return(X.svd$u[, 1:K, drop = F] %*% diag(X.svd$d[1:K], K, K) %*% t(X.svd$v[, 1:K, drop = F]))
}

SA  <- function(X, K = 2, delta = 0.5) {
  XTX  <- t(X) %*% (X)
  S  <- diag(colSums(X))
  hatB  <- solve(XTX + delta / (1 - delta) * S)%*%XTX
  temp  <- eigen(XTX + delta / (1 - delta) * S)
  Matsquareroot  <- (temp$vect) %*% diag(sqrt(temp$val))%*%t(temp$vect)
  genesvd  <- svd(Matsquareroot%*% hatB)
  genesvdrec  <- (temp$vect) %*% diag(1/sqrt(temp$val))%*%t(temp$vect) %*% genesvd$u[,1:K]%*% diag(genesvd$d[1:K]) %*% t(genesvd$v[,1:K])
  X%*%genesvdrec
} 

TSVD_opt  <- function(X, sigma = NULL){     
  n  <- nrow(X)
  p  <- ncol(X)
  svdX  <- svd(X)
  beta  <-  min(n, p)/max(n, p)
  lambdastar  <- sqrt( 2*(beta+1) + 8*beta/((beta+1+(sqrt(beta^2+14*beta+1))))  )
  
  if (is.null(sigma)){
    wbstar  <-  0.56*beta^3- 0.95*beta^2+1.82*beta+1.43
    ind  <- which((svdX$d - wbstar*median(svdX$d))>0)  
  }else{
    ind  <- which(svdX$d - lambdastar*sqrt(max(n, p)* sigma) >0)
  }      
  if (length(ind)==0){
    X.TSVD_opt  <- matrix(0,n,p)
    nbeigTSVD_opt  <- 0
  }else{  
    X.TSVD_opt  <- svdX$u[,1:length(ind), drop = F]%*% diag(svdX$d[1:length(ind)],length(ind),length(ind)) %*% t(svdX$v[, 1:length(ind), drop = F])
    nbeigTSVD_opt  = length(ind)   
  }  
  return(list(recon=X.TSVD_opt, nbeig=nbeigTSVD_opt))}

## Code for simulation
m = 50
n = 20
MU1 = outer((1:m) / m, (1:n) / n) * 5
MU2 = outer(sin((1:m) / m * 2 * pi)^8, sin((1:n) / n * 5)^8) * 10
MU3 = outer(pmin(0.2, exp(-(1:m))), pmin(0.2, exp(-(1:n)))) * 800
MU = (MU1 + MU2 + MU3)
K = 3
TRUE.svd = svd(MU)
TRUE.eig = TRUE.svd$d[1:K]
TRUE.U = TRUE.svd$u[, 1:K]
TRUE.V = TRUE.svd$v[, 1:K]

Xl = unlist(c(MU))
N = sum(Xl)
proba = Xl/N
nbsim = 1000
loss = function(X.loss) sum((X.loss/sum(X.loss) - MU/sum(MU))^2)

n.obs.list = c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000)

results.all = sapply(n.obs.list, function(n.obs) {
  results.list = lapply(1:nbsim, function(iter){
    
  #  X= matrix(rpois(length(c(MU)), n.obs * c(MU)),nrow=nrow(MU))
 # X= matrix(rpois(length(c(MU)),c(MU)),nrow=nrow(MU))
  subsamplesize = n.obs
  X = rmultinom(1, size = subsamplesize, prob = proba)
  X = matrix(X, ncol = ncol(MU))
  while(any(apply(X,1,sum)==0)|| any(apply(X,2,sum)==0)) {X= rmultinom(1, size = subsamplesize, prob = proba); X=matrix(X,ncol=ncol(MU))}
  
  TSVD.muhat <- TSVD(X, K=K)
  TSVD.loss <- loss(TSVD.muhat)
  TSVD.svd <- svd(TSVD.muhat) 
  TSVD.eig <- TSVD.svd$d[1:K]
  TSVD.U <- TSVD.svd$u[, 1:K]
  TSVD.V <- TSVD.svd$v[, 1:K]
  
  LN.res <- optishrink(X, method = "LN", k = K, center = FALSE) 
  LN.loss <- loss(LN.res$mu.hat)
  LN.eig <- LN.res$low.rank$d[1:K]
  LN.U <- LN.res$low.rank$u[, 1:K]
  LN.V <- LN.res$low.rank$v[, 1:K]
  
  SA.muhat <- SA(X, K = K, delta = 0.5)
  SA.loss <- loss(SA.muhat)
  SA.svd <- svd(SA.muhat) 
  SA.eig <- SA.svd$d[1:K]
  SA.U <-SA.svd$u[, 1:K]
  SA.V <- SA.svd$v[, 1:K]
  
  ISA.res  <-  ISA(X, delta = 0.5, noise = "Binomial", maxiter = 70, center = FALSE)
  ISA.loss  <- loss(ISA.res$mu.hat)
  ISA.nbeig <- ISA.res$nb.eigen
  ISA.eig <- ISA.res$low.rank$d[1:K]
  ISA.U <- ISA.res$low.rank$u
  ISA.V <- ISA.res$low.rank$v
  
  
  ASYMPT.res  <- optishrink(X, method = "ASYMPT", center = FALSE)
  ASYMPT.loss  <- loss(ASYMPT.res$mu.hat)
  ASYMPT.nbeig <- ASYMPT.res$nb.eigen
  ASYMPT.eig <- ASYMPT.res$low.rank$d[1:K]
  ASYMPT.U <- ASYMPT.res$low.rank$u
  ASYMPT.V <- ASYMPT.res$low.rank$v
  
  TSVD_opt.muhat <- TSVD_opt(X, sigma=NULL)
  TSVD_opt.loss <- loss(TSVD_opt.muhat$recon)
  TSVD_opt.svd <- svd(TSVD_opt.muhat$recon) 
  TSVD_opt.eig <- TSVD_opt.svd$d[1:K]
  TSVD_opt.nbeig  <- TSVD_opt.muhat$nbeig
  TSVD_opt.U <- TSVD_opt.svd$u
  TSVD_opt.V <- TSVD_opt.svd$v
  
  loss.res = rbind(TSVD=TSVD.loss, TSVD_opt=TSVD_opt.loss,  ASYMPT=ASYMPT.loss,  LN=LN.loss, SA=SA.loss, ISA=ISA.loss)
  nbeig = rbind(TSVD=K,TSVD_opt=TSVD_opt.nbeig,   ASYMPT=ASYMPT.nbeig, LN=K, SA=K,   ISA=ISA.nbeig)
  eig = rbind(TSVD=TSVD.eig, TSVD_opt=TSVD_opt.eig,  ASYMPT= ASYMPT.eig, LN=LN.eig, SA=SA.eig, ISA=ISA.eig)
 
  if(ISA.nbeig >= K) {
    ISA.RV.U = coeffRV(TRUE.U[, 1:K, drop=F], ISA.U[, 1:K, drop=F])$rv
    ISA.RV.V = coeffRV(TRUE.V[, 1:K, drop=F], ISA.V[, 1:K, drop=F])$rv
  } else {
    ISA.RV.U = NA
    ISA.RV.V = NA
  }
  
  RVTRUEU = rbind(TSVD=coeffRV(TSVD.U,TRUE.U)$rv,  TSVD_opt=coeffRV(TRUE.U[,1:min(TSVD_opt.nbeig,K),drop=F],TSVD_opt.U[,1:min(TSVD_opt.nbeig,K),drop=F])$rv,  ASYMPT=coeffRV(TRUE.U[,1:min(ASYMPT.nbeig,K),drop=F],ASYMPT.U[,1:min(ASYMPT.nbeig,K),drop=F])$rv,LN=coeffRV(TRUE.U,LN.U)$rv, SA=coeffRV(TRUE.U,SA.U)$rv,ISA= ISA.RV.U)
  RVTRUEV = rbind(TSVD=coeffRV(TSVD.V,TRUE.V)$rv, TSVD_opt=coeffRV(TRUE.V[,1:min(TSVD_opt.nbeig,K),drop=F],TSVD_opt.V[,1:min(TSVD_opt.nbeig,K),drop=F])$rv, ASYMPT=coeffRV(TRUE.V[,1:min(ASYMPT.nbeig,K),drop=F],ASYMPT.V[,1:min(ASYMPT.nbeig,K),drop=F])$rv, LN=coeffRV(TRUE.V,LN.V)$rv, SA=coeffRV(TRUE.V,SA.V)$rv, ISA=ISA.RV.V)
  
  cbind.data.frame(loss.res,eig,RVTRUEU,RVTRUEV,nbeig)
  })
  
  results = Reduce("+", results.list) / length(results.list)
  results["ISA","RVTRUEU"] = mean(sapply(results.list, function(elem)elem["ISA","RVTRUEU"]), na.rm = TRUE)
  results["ISA","RVTRUEV"] = mean(sapply(results.list, function(elem)elem["ISA","RVTRUEV"]), na.rm = TRUE)
  return(results)
})

## Results
RV.U = Reduce(rbind, results.all[5,])[,4:6]
colnames(RV.V) = c("SVD", "SA", "ISA")

RV.V = Reduce(rbind, results.all[6,])[,4:6]
colnames(RV.V) = c("SVD", "SA", "ISA")

MSE = 1000 * Reduce(rbind, results.all[1,])
colnames(MSE) = c("TSVD", "TSVD_opt", "ASYMP", "LN", "SA", "ISA")

K.EST =  Reduce(rbind, results.all[7,])[,c(2, 3, 6)]
colnames(K.EST) = c("TSVD_opt", "ASYMP", "ISA")

## Graph
vals = sqrt(c(MU1, MU2, MU3))
breaks = (min(vals) + (max(vals) - min(vals)) * seq(0, 1, by = 0.01))^2
image(1:20, 1:50, t(MU1), breaks = breaks, xlab = "", ylab = "", col = rev(gray(seq(0, 1, by = 1/99))))
image(1:20, 1:50, t(MU2), breaks = breaks, xlab = "", ylab = "", col = rev(gray(seq(0, 1, by = 1/99))))
image(1:20, 1:50, t(MU3), breaks = breaks, xlab = "", ylab = "", col = rev(gray(seq(0, 1, by = 1/99))))

library(xtable)
xtab.mse = xtable(cbind(N=as.character(n.obs.list), round(MSE[,c(5, 6, 1, 2, 3, 4)], 2)))
print(xtab.mse, include.rownames = FALSE)
xtab.RV = xtable(cbind(N=as.character(n.obs.list), round(RV.U, 2), round(RV.V, 2)))
print(xtab.RV, include.rownames = FALSE)
xtab.K = xtable(cbind(N=as.character(n.obs.list), round(K.EST, 2)))
print(xtab.K, include.rownames = FALSE)


#############################
# Regularized CA - Results Section 6 
############################# 

rm(list = ls())
library(FactoMineR)
library(denoiseR)

## Functions

# Estimator LN 
LN <- function (X, K = 2){
  P <- as.matrix(X/sum(X))
  Rc <- apply(P, 2, sum)
  Rr <- apply(P, 1, sum)
  S = diag(1/sqrt(Rr))%*%(P-Rr%*%t(Rc))%*%diag(1/sqrt(Rc))
  svdRes = svd(S)
  n = nrow(X) ; p = ncol(X) 
  sigma2 = sum(svdRes$d[-c(1:K)]^2)/((n-1)*(p-1) - (n-1)*K - (p-1)*K + K^2) 
  lambda.shrinked = (svdRes$d[1:K]^2 - n*(p/min(p, (n-1)))*sigma2)/svdRes$d[1:K] 
  if(K == 1)
    recon = (svdRes$u[, 1] *lambda.shrinked)%*%t(svdRes$v[, 1])
  if(K > 1)
    recon = svdRes$u[, 1:K] %*%diag(lambda.shrinked)%*%t(svdRes$v[, 1:K])
  recon <- sum(X)*(sweep(sweep(recon,1,sqrt(Rr),FUN="*"),2,sqrt(Rc),FUN="*") + Rr%*%t(Rc))
  rownames(recon) = rownames(X)
  colnames(recon) = colnames(X) 
  U=diag(1/sqrt(Rr))%*%sweep(svdRes$u[, 1:K], 2, lambda.shrinked, FUN = "*")
  V=diag(1/sqrt(Rc))%*%sweep(svdRes$v[, 1:K], 2, lambda.shrinked, FUN = "*")
  res <- list(eig =  (lambda.shrinked^2)[1:K], row = U[,1:K,drop=FALSE], col = V[,1:K,drop=FALSE],recon=recon,svdRes=svdRes)
  return(res)  
}

#  Estimator SA
SA = function(X, K = 2, delta = 0.5) { 
  N = sum(X)
  r = rowSums(X)
  c = colSums(X)
  RC = matrix(r, length(r), 1) %*% matrix(c, 1, length(c))
  M = diag(1/sqrt(r)) %*% (X - RC / N) %*% diag(1/sqrt(c))
  M.svd = svd(M) 
  SM = diag(colSums(X/RC))*delta / (1 - delta)
  Bhat = solve(t(M) %*% M + SM, t(M) %*% M) 
  temp=eigen(t(M) %*% M + SM)
  Matsquareroot= (temp$vect) %*% diag(sqrt(temp$val))%*%t(temp$vect)
  genesvd=svd(Matsquareroot%*% Bhat)
  genesvdrec= (temp$vect) %*% diag(1/sqrt(temp$val))%*%t(temp$vect) %*% genesvd$u[,1:K]%*% diag(genesvd$d[1:K]) %*% t(genesvd$v[,1:K])
  hY=M %*% Bhat    
  recon <- (sweep(sweep(hY,1,sqrt(r),FUN="*"),2,sqrt(c),FUN="*") + r%*%t(c)/N)
  svdRes = svd(hY)
  U=diag(1/sqrt(r))%*%sweep(svdRes$u, 2, svdRes$d, FUN = "*")
  V=diag(1/sqrt(c))%*%sweep(svdRes$v, 2, svdRes$d, FUN = "*")
  res <- list(eig =  (svdRes$d^2)[1:K], row = U[,1:K,drop=FALSE], col = V[,1:K,drop=FALSE],recon=recon, svdRes=svdRes)
  return(res)
} 

## Code for simulation
don <-  read.table("http://factominer.free.fr/docs/perfume.txt",header=T,sep="\t",row.names=1)
X.TRUE=t(don)
colnames(X.TRUE)[4]="cinema"

# True population CA
K = 2
CA.TRUE=CA(X.TRUE, ncp=K,graph=F)
TRUE.eig=CA.TRUE$eig[1:K,1]
TRUE.row=CA.TRUE$row$coord[,1:K]
TRUE.col=CA.TRUE$col$coord[,1:K]

# Sample from the population data
Xl=unlist(c(X.TRUE))
N = sum(Xl)
proba = Xl/N
subsamplesize = 200
nbsim = 1000

results.list = lapply(1:nbsim, function(iter){
  X= rmultinom(1, size = subsamplesize, prob = proba)
  X=matrix(X,ncol=ncol(X.TRUE))
  while(any(apply(X,1,sum)==0)|| any(apply(X,2,sum)==0)) {X= rmultinom(1, size = subsamplesize, prob = proba);
  X=matrix(X,ncol=ncol(X.TRUE))}
  rownames(X)=rownames(X.TRUE)
  colnames(X)=colnames(X.TRUE)
  
  ## CA
  CA.X = CA(X,graph=F,ncp=K)
  CA.muhat = reconst(CA.X,K)
  CA.eig = CA.X$eig[1:K,1]
  CA.row = CA.X$row$coord
  CA.col = CA.X$col$coord
  
  ##LN
  LN.muhat = LN(X,K)
  LN.eig = LN.muhat$eig
  LN.row = LN.muhat$row
  LN.col = LN.muhat$col
  
  ## SA
  SA.muhat = SA(X,  K=K, delta=0.5)
  SA.eig = SA.muhat$eig
  SA.row = SA.muhat$row
  SA.col = SA.muhat$col
  
  ## ISA
  ISA.res = ISA(X, delta = 0.3, maxiter=300, transformation = "CA")
  ISA.eig = ISA.res$low.rank$d[1:K]
  ISA.nbeig = ISA.res$nb.eigen
  res.isa.ca <- CA(ISA.res$mu.hat, graph = FALSE)
  ISA.row = res.isa.ca$row$coord
  ISA.col = res.isa.ca$col$coord
  
  nbeig = rbind(CA=K,LN=K, SA=K,ISA=ISA.nbeig)
  eig = rbind(CA=CA.eig, LN=LN.eig, SA=SA.eig, ISA=ISA.eig)
  RVrow = rbind(CA=coeffRV(TRUE.row,CA.row)$rv, LN=coeffRV(TRUE.row,LN.row)$rv,  SA=coeffRV(TRUE.row,SA.row)$rv,  ISA=coeffRV(TRUE.row[,1:min(ISA.nbeig,K),drop=F],ISA.row[,1:min(ISA.nbeig,K),drop=F])$rv) 
  RVcol = rbind(CA=coeffRV(TRUE.col,CA.col)$rv, LN=coeffRV(TRUE.col,LN.col)$rv, SA=coeffRV(TRUE.col,SA.col)$rv, ISA=coeffRV(TRUE.col[,1:min(ISA.nbeig,K),drop=F],ISA.col[,1:min(ISA.nbeig,K),drop=F])$rv)
  cbind.data.frame(eig, nbeig, RVrow, RVcol)
})

results = Reduce("+", results.list) / length(results.list)

## Graphs to compare CA and ISA
X= rmultinom(1, size = subsamplesize, prob = proba)
X=matrix(X,ncol=ncol(X.TRUE))
while(any(apply(X,1,sum)==0)|| any(apply(X,2,sum)==0)) {X= rmultinom(1, size = subsamplesize, prob = proba);
X=matrix(X,ncol=ncol(X.TRUE))}
ISA.res = ISA(X, delta = 0.3, maxiter=300, transformation = "CA")
rownames(X)=rownames(ISA.res$mu.hat)=rownames(X.TRUE)
colnames(X)=colnames(ISA.res$mu.hat)=colnames(X.TRUE)

CA.res=CA(X,graph=F)
CA.res=CA(reconst(CA.res,2),graph=F)
plot(CA.res,title="CA",autoLab="no",cex=0.6,selectRow="contrib 20")
res.isa.ca <- CA(ISA.res$mu.hat, graph = FALSE)
plot(res.isa.ca,title="Regularized CA",autoLab="no",cex=0.6,selectRow="contrib 20")


 
