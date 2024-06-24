rm(list=ls())
setwd("D:\\LearnFiles\\Research paper\\idea\\TodoProjects\\SpatialGFM\\Rcode\\")

setwd("Z:/ZXC/空间多模态数据分析")
source("selfdefinedFuncs.R")




# Gaussian + Poisson ------------------------------------------------------
library(GFM)
pveclist <- list('gaussian'=c(50, 150),'poisson'=c(50, 150))
q <- 6
sigmavec <- rep(1,length(pveclist))
N <- 100
pvec <- unlist(pveclist)
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
  
  datlist <- gendata(pveclist = pveclist, seed = i, n = 300,d = 3,
                     q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
                     sigmavec=sigmavec, sigma_eps=1)
  str(datlist)
  XList <- datlist$XList
  Z <- datlist$Z
  numvarmat <- datlist$numvarmat
  Uplist <- list()
  k <- 1
  for(id in 1:length(datlist$Uplist)){
    for(im in 1:length(datlist$Uplist[[id]])){
      Uplist[[k]] <- datlist$Uplist[[id]][[im]]
      k <- k + 1
    }
  }
  types <- datlist$types
  
  tic <- proc.time()
  rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat, init='LFM', add_IC_iter=F)
  toc <- proc.time()
  time_smgfm <- toc[3] - tic[3]
  Lambda0 <- rep(datlist$sigma_eps, sum(sapply(XList, ncol)))
  hLambda <- 1/unlist(rlist$invLambdaf)
  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
  metricList$beta_norm[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- Reduce(cbind,datlist$beta0List)))
  metricList$Lambda_norm[i,1] <- normvec(Lambda0 - hLambda)
  metricList$Upsilon_tr[i,1] <- meanTr(hUplist, Uplist)
  metricList$F_tr[i,1] <- measurefun(rlist$M, datlist$F0)
  try({ metricList$U_tr[i,1] <-measurefun(Reduce(cbind,rlist$Xif), Reduce(cbind, datlist$U0))}, silent=T)
  metricList$timeMat[i,1] <- time_smgfm
  
  ## GFM
  tic <- proc.time()
  res_gfm <- gfm(XList, types=types, q=q)
  toc <- proc.time()
  time_gfm <- toc[3] - tic[3]
  metricList$Upsilon_tr[i,2] <-meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), Uplist)
  metricList$F_tr[i,2] <- measurefun(res_gfm$hH, datlist$F0)
  metricList$timeMat[i,2] <- time_gfm
  
  
  ### MRRR
  family_use <- list(gaussian(), poisson())
  familygroup <- lapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
  res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
                       familygroup = unlist(familygroup), maxIter=10) #100
  
  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- Reduce(cbind,datlist$beta0List)))
  metricList$timeMat[i,3] <- res_mrrr$time_use
  
  
  ## LFR
  XtmpList <- XList
  XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
  XcList <- list(Reduce(cbind, XtmpList))
  try({
    res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 500
    metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
    metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf-  Reduce(cbind,datlist$beta0List)))
    metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, Uplist)
    metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
    metricList$timeMat[i,4] <- res_rf$time.use
  }, silent=T)
  # ## generalizedPCA
  # library(generalizedPCA)
  # t1 <- proc.time()
  # res_gPCA <- generalizedPCA(datlist$X, k=q, family = "gaussian", quiet = F)
  # t2 <- proc.time()
  # time_gPCA <- t2[3]-t1[3]
  # metricList$Upsilon_tr[i,5] <-meanTr(mat2list(cbind(res_gPCA$mu,res_gPCA$U), pvec), datlist$Uplist)
  # metricList$F_tr[i,5] <- measurefun(res_gPCA$PCs, datlist$F0)
  # metricList$timeMat[i,5] <- time_gPCA
}
# install.packages("remotes")
save(metricList, file=paste0('metricList_simu3_np.rds'))
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)






# Gaussian + Binomial ------------------------------------------------------
library(GFM)
pveclist <- list('gaussian'=c(50, 150),'binomial'=c(50, 150))
q <- 6
sigmavec <- rep(1,length(pveclist))
N <- 5
pvec <- unlist(pveclist)
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
  
  datlist <- gendata(pveclist = pveclist, seed = i, n = 300,d = 3,
                     q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
                     sigmavec=sigmavec, sigma_eps=1)
  str(datlist)
  XList <- datlist$XList
  Z <- datlist$Z
  numvarmat <- datlist$numvarmat
  Uplist <- list()
  k <- 1
  for(id in 1:length(datlist$Uplist)){
    for(im in 1:length(datlist$Uplist[[id]])){
      Uplist[[k]] <- datlist$Uplist[[id]][[im]]
      k <- k + 1
    }
  }
  types <- datlist$types
  
  tic <- proc.time()
  rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat, init='LFM', add_IC_iter=F)
  toc <- proc.time()
  time_smgfm <- toc[3] - tic[3]
  Lambda0 <- rep(datlist$sigma_eps, sum(sapply(XList, ncol)))
  hLambda <- 1/unlist(rlist$invLambdaf)
  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
  metricList$beta_norm[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- Reduce(cbind,datlist$beta0List)))
  metricList$Lambda_norm[i,1] <- normvec(Lambda0 - hLambda)
  metricList$Upsilon_tr[i,1] <- meanTr(hUplist, Uplist)
  metricList$F_tr[i,1] <- measurefun(rlist$M, datlist$F0)
  try({ metricList$U_tr[i,1] <-measurefun(Reduce(cbind,rlist$Xif), Reduce(cbind, datlist$U0))}, silent=T)
  metricList$timeMat[i,1] <- time_smgfm
  
  ## GFM
  tic <- proc.time()
  res_gfm <- gfm(XList, types=types, q=q)
  toc <- proc.time()
  time_gfm <- toc[3] - tic[3]
  metricList$Upsilon_tr[i,2] <-meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), Uplist)
  metricList$F_tr[i,2] <- measurefun(res_gfm$hH, datlist$F0)
  metricList$timeMat[i,2] <- time_gfm
  
  
  ### MRRR
  family_use <- list(gaussian(), poisson())
  familygroup <- lapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
  res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
                       familygroup = unlist(familygroup), maxIter=10) #100
  
  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- Reduce(cbind,datlist$beta0List)))
  metricList$timeMat[i,3] <- res_mrrr$time_use
  
  
  ## LFR
  XtmpList <- XList
  XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
  XcList <- list(Reduce(cbind, XtmpList))
  try({
    res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 500
    metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
    metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf-  Reduce(cbind,datlist$beta0List)))
    metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, Uplist)
    metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
    metricList$timeMat[i,4] <- res_rf$time.use
  }, silent=T)
  # ## generalizedPCA
  # library(generalizedPCA)
  # t1 <- proc.time()
  # res_gPCA <- generalizedPCA(datlist$X, k=q, family = "gaussian", quiet = F)
  # t2 <- proc.time()
  # time_gPCA <- t2[3]-t1[3]
  # metricList$Upsilon_tr[i,5] <-meanTr(mat2list(cbind(res_gPCA$mu,res_gPCA$U), pvec), datlist$Uplist)
  # metricList$F_tr[i,5] <- measurefun(res_gPCA$PCs, datlist$F0)
  # metricList$timeMat[i,5] <- time_gPCA
}

save(metricList, file=paste0('metricList_simu3_nb.rds'))
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)






# Poisson + Binomial ------------------------------------------------------
library(GFM)
pveclist <- list('poisson'=c(50, 150),'binomial'=c(50, 150))
q <- 6
sigmavec <- rep(1,length(pveclist))
N <- 5
pvec <- unlist(pveclist)
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
  
  datlist <- gendata(pveclist = pveclist, seed = i, n = 300,d = 3,
                     q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
                     sigmavec=sigmavec, sigma_eps=1)
  str(datlist)
  XList <- datlist$XList
  Z <- datlist$Z
  numvarmat <- datlist$numvarmat
  Uplist <- list()
  k <- 1
  for(id in 1:length(datlist$Uplist)){
    for(im in 1:length(datlist$Uplist[[id]])){
      Uplist[[k]] <- datlist$Uplist[[id]][[im]]
      k <- k + 1
    }
  }
  types <- datlist$types
  
  tic <- proc.time()
  rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat, init='LFM', add_IC_iter=F)
  toc <- proc.time()
  time_smgfm <- toc[3] - tic[3]
  Lambda0 <- rep(datlist$sigma_eps, sum(sapply(XList, ncol)))
  hLambda <- 1/unlist(rlist$invLambdaf)
  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
  metricList$beta_norm[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- Reduce(cbind,datlist$beta0List)))
  metricList$Lambda_norm[i,1] <- normvec(Lambda0 - hLambda)
  metricList$Upsilon_tr[i,1] <- meanTr(hUplist, Uplist)
  metricList$F_tr[i,1] <- measurefun(rlist$M, datlist$F0)
  try({ metricList$U_tr[i,1] <-measurefun(Reduce(cbind,rlist$Xif), Reduce(cbind, datlist$U0))}, silent=T)
  metricList$timeMat[i,1] <- time_smgfm
  
  ## GFM
  tic <- proc.time()
  res_gfm <- gfm(XList, types=types, q=q)
  toc <- proc.time()
  time_gfm <- toc[3] - tic[3]
  metricList$Upsilon_tr[i,2] <-meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), Uplist)
  metricList$F_tr[i,2] <- measurefun(res_gfm$hH, datlist$F0)
  metricList$timeMat[i,2] <- time_gfm
  
  
  ### MRRR
  family_use <- list( poisson(), binomial())
  familygroup <- lapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
  res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
                       familygroup = unlist(familygroup), maxIter=10) #100
  
  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- Reduce(cbind,datlist$beta0List)))
  metricList$timeMat[i,3] <- res_mrrr$time_use
  
  
  ## LFR
  XtmpList <- XList
  XtmpList[[1]] <- log(1+XList[[1]]) # log-transformation
  XcList <- list(Reduce(cbind, XtmpList))
  try({
    res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 500
    metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
    metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf-  Reduce(cbind,datlist$beta0List)))
    metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, Uplist)
    metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
    metricList$timeMat[i,4] <- res_rf$time.use
  }, silent=T)
  # ## generalizedPCA
  # library(generalizedPCA)
  # t1 <- proc.time()
  # res_gPCA <- generalizedPCA(datlist$X, k=q, family = "gaussian", quiet = F)
  # t2 <- proc.time()
  # time_gPCA <- t2[3]-t1[3]
  # metricList$Upsilon_tr[i,5] <-meanTr(mat2list(cbind(res_gPCA$mu,res_gPCA$U), pvec), datlist$Uplist)
  # metricList$F_tr[i,5] <- measurefun(res_gPCA$PCs, datlist$F0)
  # metricList$timeMat[i,5] <- time_gPCA
}

save(metricList, file=paste0('metricList_simu3_pb.rds'))
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)







# MultiModality Gaussian, Poisson and Binomial data ---------------------------------------------

library(GFM)
pveclist <- list('gaussian'=c(50, 150),'poisson'=c(50, 150),
                 'binomial'=c(100,60))
q <- 6
sigmavec <- rep(1,3)
N <- 100
pvec <- unlist(pveclist)
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
  
  datlist <- gendata(pveclist = pveclist, seed = i, n = 300,d = 3,
                           q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
                           sigmavec=sigmavec, sigma_eps=1)
  str(datlist)
  XList <- datlist$XList
  Z <- datlist$Z
  numvarmat <- datlist$numvarmat
  Uplist <- list()
  k <- 1
  for(id in 1:length(datlist$Uplist)){
    for(im in 1:length(datlist$Uplist[[id]])){
      Uplist[[k]] <- datlist$Uplist[[id]][[im]]
      k <- k + 1
    }
  }
  types <- datlist$types
  
  tic <- proc.time()
  rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat, init='LFM', add_IC_iter=F)
  toc <- proc.time()
  time_smgfm <- toc[3] - tic[3]
  Lambda0 <- rep(datlist$sigma_eps, sum(sapply(XList, ncol)))
  hLambda <- 1/unlist(rlist$invLambdaf)
  hUplist <- lapply(seq_along(rlist$Bf), function(m) cbind(rlist$muf[[m]], rlist$Bf[[m]]))
  metricList$beta_norm[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- Reduce(cbind,datlist$beta0List)))
  metricList$Lambda_norm[i,1] <- normvec(Lambda0 - hLambda)
  metricList$Upsilon_tr[i,1] <- meanTr(hUplist, Uplist)
  metricList$F_tr[i,1] <- measurefun(rlist$M, datlist$F0)
  try({ metricList$U_tr[i,1] <-measurefun(Reduce(cbind,rlist$Xif), Reduce(cbind, datlist$U0))}, silent=T)
  metricList$timeMat[i,1] <- time_smgfm
  
  ## GFM
  tic <- proc.time()
  res_gfm <- gfm(XList, types=types, q=q)
  toc <- proc.time()
  time_gfm <- toc[3] - tic[3]
  metricList$Upsilon_tr[i,2] <-meanTr(mat2list(cbind(res_gfm$hmu,res_gfm$hB), pvec), Uplist)
  metricList$F_tr[i,2] <- measurefun(res_gfm$hH, datlist$F0)
  metricList$timeMat[i,2] <- time_gfm
  
  
  ### MRRR
  family_use <- list(gaussian(), poisson(), binomial())
  familygroup <- sapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
  res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
                       familygroup = unlist(familygroup), maxIter=100) #
  
  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- Reduce(cbind,datlist$beta0List)))
  metricList$timeMat[i,3] <- res_mrrr$time_use
  
  
  ## LFR
  XtmpList <- XList
  XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
  XcList <- list(Reduce(cbind, XtmpList))
  try({
    res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 500
    metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
    metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf-  Reduce(cbind,datlist$beta0List)))
    metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, Uplist)
    metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
    metricList$timeMat[i,4] <- res_rf$time.use
  }, silent=T)
  # ## generalizedPCA
  # library(generalizedPCA)
  # t1 <- proc.time()
  # res_gPCA <- generalizedPCA(datlist$X, k=q, family = "gaussian", quiet = F)
  # t2 <- proc.time()
  # time_gPCA <- t2[3]-t1[3]
  # metricList$Upsilon_tr[i,5] <-meanTr(mat2list(cbind(res_gPCA$mu,res_gPCA$U), pvec), datlist$Uplist)
  # metricList$F_tr[i,5] <- measurefun(res_gPCA$PCs, datlist$F0)
  # metricList$timeMat[i,5] <- time_gPCA
  print(i)
}

save(metricList, file=paste0('metricList_simu3_npb.rds'))
sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)
