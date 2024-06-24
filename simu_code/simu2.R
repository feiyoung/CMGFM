setwd("D:\\WorkSpace\\ScientificResearch\\idea\\SpatialGFM\\Rcode")
source("selfdefinedFuncs.R")

rhovecList  <- list(c(1, 0.2), c(1.5, 0.2), c(1, 0.5))
ja <- commandArgs(TRUE)
ja <- as.integer(ja)

# ja <- 3

library(GFM)
q <- 6
rho_vec<- rhovecList[[ja]]
N <- 100
pvec <- c(150, 200)
methodNames <- c("CMGFM", "GFM", "MRRR", "LFR" , "GPCA", "COAP")
n_methods <- length(methodNames)
metricList <- list(F_tr = matrix(NA,N, n_methods), 
                   U_tr = matrix(NA,N, n_methods), 
                   Upsilon_tr = matrix(NA, N, n_methods),
                   beta_norm=matrix(NA, N, n_methods),
                   beta_normres=matrix(NA, N, n_methods),
                   Lambda_norm=matrix(NA, N, n_methods),
                   timeMat = matrix(NA, N, n_methods))
for(ii in seq_along(metricList)) colnames(metricList[[ii]]) <- methodNames

for(i in 1:N){
  #i <- 1
  
  datlist <- gendata_pois(seed = i, n = 300, pvec = pvec, d = 3,
                          q = q, rho = rep(rho_vec[1],length(pvec)), rho_z=rho_vec[2],
                          sigmavec=rep(1, length(pvec)), sigma_eps=1)
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
  metricList$beta_normres[i,1] <- normvec(as.vector(Reduce(cbind,rlist$betaf)- datlist$beta))/ normvec(as.vector(datlist$beta))
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
                       familygroup = unlist(familygroup), maxIter=100) #100

  metricList$Upsilon_tr[i,3] <-meanTr(mat2list(cbind(res_mrrr$hmu,res_mrrr$hB), pvec), datlist$Uplist)
  metricList$F_tr[i,3] <- measurefun(res_mrrr$hH, datlist$F0)
  metricList$beta_norm[i,3] <- normvec(as.vector(res_mrrr$beta- datlist$beta))
  metricList$beta_normres[i,3] <- normvec(as.vector(res_mrrr$beta- datlist$beta))/ normvec(as.vector(datlist$beta))
  metricList$timeMat[i,3] <- res_mrrr$time_use


  ## LFR
  XcList <- list(Reduce(cbind, XList)) # 500
  res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=TRUE)
  metricList$F_tr[i,4] <- measurefun(res_rf$F, datlist$F0)
  metricList$beta_norm[i,4] <- normvec(as.vector(res_rf$beta_rf)- datlist$beta)
  metricList$beta_normres[i,4] <- normvec(as.vector(res_rf$beta_rf)- datlist$beta)/ normvec(as.vector(datlist$beta))
  metricList$Upsilon_tr[i,4] <- meanTr(res_rf$hUplist, datlist$Uplist)
  metricList$Lambda_norm[i,4] <- normvec(res_rf$lambdavec - hLambda)
  metricList$timeMat[i,4] <- res_rf$time.use


  ## generalizedPCA
  # install.packages("GFM")
  library(generalizedPCA)
  t1 <- proc.time()
  res_gPCA <- generalizedPCA_here(datlist$X, k=q, family = "poisson", quiet = F)
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
  metricList$beta_normres[i,6] <- normvec(as.vector(reslist$hbeta) - datlist$beta)/ normvec(as.vector(datlist$beta))
  metricList$Lambda_norm[i,6] <- normvec(as.vector(1/reslist$invLambda) - hLambda)
}

save(metricList, file=paste0('metricList_simu1_pois_rhovec', ja ,'.rds'))

sapply(metricList, colMeans, na.rm=T)
sapply(metricList, colSD)

