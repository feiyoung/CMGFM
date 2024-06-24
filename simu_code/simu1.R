rm(list=ls())
setwd("D:\\LearnFiles\\Research paper\\idea\\TodoProjects\\SpatialGFM\\Rcode\\")

setwd("Z:/ZXC/空间多模态数据分析")
source("selfdefinedFuncs.R")


sigma2vec <- c(0, 0.5, 1, 2)
ja <- 1 #每个sigma值需手动调整（1-4）
library(GFM)
q <- 6
sigmavec <- rep(sigma2vec[ja],2)
N <- 10
pvec <- c(50, 150)
methodNames <- c("CMGFM", "GFM", "MRRR", "LFR" , "GPCA", "COAP")
n_methods <- length(methodNames)
metricList <- list(F_tr = matrix(NA,N, n_methods), 
                   U_tr = matrix(NA,N, n_methods), 
                   Upsilon_tr = matrix(NA, N, n_methods),
                   beta_norm=matrix(NA, N, n_methods),
                   Lambda_norm=matrix(NA, N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
for(ii in seq_along(metricList)) colnames(metricList[[ii]]) <- methodNames
for(i in 1:N){
  #i <- 1
  
  datlist <- gendata_gauss(seed = i, n = 300, pvec = pvec, d = 3,
                                          q = q, rho = rep(1,length(pvec)), rho_z=1,
                                          sigmavec=sigmavec, sigma_eps=1)
  str(datlist)
  XList <- list(datlist$X)
  Z <- datlist$Z
  numvarmat <- datlist$numvarmat
  
  
  tic <- proc.time()
  rlist <- CMGFM(XList, Z, types="gaussian", q=q, numvarmat, init='random', add_IC_iter=F)
  toc <- proc.time()
  time_smgfm <- toc[3] - tic[3]
  Lambda0 <- rep(datlist$sigma_eps, ncol(XList[[1]]))
  hLambda <- 1/unlist(rlist$invLambdaf)
  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
  metricList$beta_norm[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- datlist$beta))
  metricList$Lambda_norm[i,1] <- normvec(Lambda0 - hLambda)
  metricList$Upsilon_tr[i,1] <- meanTr(hUplist, datlist$Uplist)
  metricList$F_tr[i,1] <- measurefun(rlist$M, datlist$F0)
  try({ metricList$U_tr[i,1] <-measurefun(Reduce(cbind,rlist$Xif), datlist$U0)}, silent=T)
  metricList$timeMat[i,1] <- time_smgfm
  
  ## GFM
  tic <- proc.time()
  res_gfm <- gfm(XList, types=c("gaussian"), q=q)
  toc <- proc.time()
  time_gfm <- toc[3] - tic[3]
  metricList$Upsilon_tr[i,2] <-meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), datlist$Uplist)
  metricList$F_tr[i,2] <- measurefun(res_gfm$hH, datlist$F0)
  metricList$timeMat[i,2] <- time_gfm


  ### MRRR
  familygroup <- rep(1, ncol(datlist$X))

  res_mrrr <- mrrr_run(datlist$X,  Z = Z,numvarmat, rank0=q, family=list(gaussian()),
                       familygroup = unlist(familygroup), maxIter=100) #

  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), datlist$Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- datlist$beta))
  metricList$timeMat[i,3] <- res_mrrr$time_use


  ## LFR
  XcList <- list(Reduce(cbind, XList))
  res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE)
  metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
  metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf- datlist$beta))
  metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, datlist$Uplist)
  metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
  metricList$timeMat[i,4] <- res_rf$time.use


  ## generalizedPCA
  library(generalizedPCA)
  t1 <- proc.time()
  res_gPCA <- generalizedPCA(datlist$X, k=q, family = "gaussian", quiet = F)
  t2 <- proc.time()
  time_gPCA <- t2[3]-t1[3]
  metricList$Upsilon_tr[i,5] <-meanTr(mat2list(cbind(res_gPCA$mu,res_gPCA$U), pvec), datlist$Uplist)
  metricList$F_tr[i,5] <- measurefun(res_gPCA$PCs, datlist$F0)
  metricList$timeMat[i,5] <- time_gPCA
}

save(metricList, file=paste0('metricList_simu1_gauss_sigma', sigma2vec[ja] ,'.rds'))
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)




# Poisson case ------------------------------------------------------------
sigma2vec <- c(0, 0.5, 1, 2)
ja <- 1 #每个sigma值需手动调整（1-4）
library(GFM)
q <- 6
sigmavec <- rep(sigma2vec[ja],2)
N <- 100
pvec <- c(50, 150)
methodNames <- c("CMGFM", "GFM", "MRRR", "LFR" , "GPCA", "COAP")
n_methods <- length(methodNames)
metricList <- list(F_tr = matrix(NA,N, n_methods), 
                   U_tr = matrix(NA,N, n_methods), 
                   Upsilon_tr = matrix(NA, N, n_methods),
                   beta_norm=matrix(NA, N, n_methods),
                   Lambda_norm=matrix(NA, N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
for(ii in seq_along(metricList)) colnames(metricList[[ii]]) <- methodNames
for(i in 1:N){
  #i <- 1
  
  datlist <- gendata_pois(seed = i, n = 300, pvec = pvec, d = 3,
                           q = q, rho = rep(2,length(pvec)), rho_z=0.2,
                           sigmavec=sigmavec, sigma_eps=1)
  str(datlist)
  XList <- list(datlist$X)
  max(datlist$X)
  Z <- datlist$Z
  numvarmat <- datlist$numvarmat
  
  library(CMGFM)
  tic <- proc.time()
  rlist <- CMGFM(XList, Z, types="poisson", q=q, numvarmat, init='LFM', epsELBO = 1e-8)
  toc <- proc.time()
  time_smgfm <- toc[3] - tic[3]
  Lambda0 <- rep(datlist$sigma_eps, ncol(XList[[1]]))
  hLambda <- 1/unlist(rlist$invLambdaf)
  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
  metricList$beta_norm[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- datlist$beta))
  metricList$Lambda_norm[i,1] <- normvec(Lambda0 - hLambda)
  metricList$Upsilon_tr[i,1] <- meanTr(hUplist, datlist$Uplist)
  metricList$F_tr[i,1] <- measurefun(rlist$M, datlist$F0)
  try({ metricList$U_tr[i,1] <-measurefun(Reduce(cbind,rlist$Xif), datlist$U0)}, silent=T)
  metricList$timeMat[i,1] <- time_smgfm
  
  ## GFM
  tic <- proc.time()
  res_gfm <- gfm(XList, types=c("poisson"), q=q)
  toc <- proc.time()
  time_gfm <- toc[3] - tic[3]
  metricList$Upsilon_tr[i,2] <-meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), datlist$Uplist)
  metricList$F_tr[i,2] <- measurefun(res_gfm$hH, datlist$F0)
  metricList$timeMat[i,2] <- time_gfm


  ### MRRR
  familygroup <- rep(1, ncol(datlist$X))

  res_mrrr <- mrrr_run(datlist$X,  Z = Z,numvarmat, rank0=q, family=list(poisson()),
                       familygroup = unlist(familygroup), maxIter=100) #

  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), datlist$Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- datlist$beta))
  metricList$timeMat[i,3] <- res_mrrr$time_use


  ## LFR
  XcList <- list(Reduce(cbind, XList))
  res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=TRUE)
  metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
  metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf)- datlist$beta)
  metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, datlist$Uplist)
  metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
  metricList$timeMat[i,4] <- res_rf$time.use


  ## generalizedPCA
  library(generalizedPCA)
  t1 <- proc.time()
  res_gPCA <- generalizedPCA(datlist$X, k=q, family = "poisson", quiet = F)
  t2 <- proc.time()
  time_gPCA <- t2[3]-t1[3]
  metricList$Upsilon_tr[i,5] <-meanTr(mat2list(cbind(res_gPCA$mu,res_gPCA$U), pvec), datlist$Uplist)
  metricList$F_tr[i,5] <- measurefun(res_gPCA$PCs, datlist$F0)
  metricList$timeMat[i,5] <- time_gPCA

  ### COAP
  tic <- proc.time()
  reslist <- COAP_run(datlist$X, Z = Z, numvarmat, q= q, rank = ncol(Z)+1, epsELBO = 1e-7)
  toc <- proc.time()
  time_apfactor <- toc[3] - tic[3]
  metricList$timeMat[i,6] <- time_apfactor
  metricList$F_tr[i,6] <- measurefun(reslist$H, datlist$F0)
  metricList$Upsilon_tr[i,6] <- meanTr(mat2list(cbind(reslist$bbeta[,1], reslist$B), pvec),datlist$Uplist)
  metricList$beta_norm[i, 6] <- normvec(as.vector(reslist$hbeta) - as.vector( datlist$beta))
  metricList$Lambda_norm[i,6] <- normvec(as.vector(1/reslist$invLambda) - hLambda)
  print(i)
  }

save(metricList, file=paste0('metricList_simu1_pois_sigma', sigma2vec[ja] ,'.rds'))
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)

