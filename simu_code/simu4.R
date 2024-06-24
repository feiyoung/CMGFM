setwd("D:\\WorkSpace\\ScientificResearch\\idea\\SpatialGFM\\Rcode")
source("selfdefinedFuncs.R")
library(GFM)

q <- 6
sigmavec <- rep(1, 3)
N <- 10
methodNames <- c("CMGFM", "GFM", "MRRR", "LFR", "GPCA", "COAP")
n_methods <- length(methodNames)



# Fix p, while change n ---------------------------------------------------
nvec <- c(500, 1000, 3000, 5000, 8000, 10000)
pvec2 <- c(50, 150) #c(50, 150) # 
pveclist <- list("gaussian"=pvec2, 'poisson'=pvec2,
                 'binomial' = pvec2)
timeArray <- array(dim=c(N, n_methods,  length(nvec)))
colnames(timeArray) <- methodNames
for(jn in c(1, 5:length(nvec))){
  # jn <- 1
  n <- nvec[jn]
  message('jn  = ', jn, ', n=', n)
  for(i in 1:N){
    #i <- 1
    
    datlist <- gendata(pveclist = pveclist, seed = i, n = n,d = 3,
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
    
    library(CMGFM)
    tic <- proc.time()
    rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat, init='LFM', add_IC_iter=F)
    toc <- proc.time()
    time_smgfm <- toc[3] - tic[3]
    timeArray[i,1,jn] <- time_smgfm
    
    ## GFM
    tic <- proc.time()
    res_gfm <- gfm(XList, types=types, q=q, dc_eps = 1e-20)
    toc <- proc.time()
    time_gfm <- toc[3] - tic[3]
    timeArray[i,2,jn] <- time_gfm
    
    
    ### MRRR
    family_use <- list(gaussian(), poisson(), binomial())
    familygroup <- sapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
    res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
                         familygroup = unlist(familygroup), maxIter=100) #
    
    timeArray[i,3,jn] <- res_mrrr$time_use
    
    
    # ## LFR
    # XtmpList <- XList
    # XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
    # XcList <- list(Reduce(cbind, XtmpList))
    # try({
    #   res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 500
    #   timeArray[i,4,jn] <- res_rf$time.use
    # }, silent=T)
    ## generalizedPCA
    library(generalizedPCA)
    Xall <- cbind(datlist$XList[[2]], datlist$XList[[2]], datlist$XList[[2]])
    t1 <- proc.time()
    res_gPCA <- generalizedPCA_here(Xall, k=q, family = "poisson", quiet = F)
    t2 <- proc.time()
    time_gPCA <- t2[3]-t1[3]
    timeArray[i,5,jn] <- time_gPCA
    
    ### COAP
    tic <- proc.time()
    reslist <- COAP_run(Xall, Z = Z, matrix(colSums(numvarmat),1,2), q= q, rank = ncol(Z)+1, epsELBO = 1e-20)
    toc <- proc.time()
    time_apfactor <- toc[3] - tic[3]
    timeArray[i,6,jn] <- time_apfactor
    
  }
  save(timeArray, file=paste0('timeArray_simu4_exceptLFR_npb_changen.rds'))
}



apply(timeArray, c(2,3) , mean, na.rm=T)
apply(timeArray, c(2,3) , sd, na.rm=T)

### results combing
setwd("/share/analysisdata/liuw/Others/CMGFM/simu")
load("timeArray_simu4_exceptLFR_npb.rds")
load("timeArray_simu4_exceptLFR_npb_2start.rds")
apply(timeArray, c(2,3) , mean, na.rm=T)

N <- 10
### Run LFR
for(i in 1:N){
  #i <- 1
 for(jn in 2:length(nvec)){
  # jn <- 1
  n <- nvec[jn]
  message('jn  = ', jn, ', n=', n)
 
    
    datlist <- gendata(pveclist = pveclist, seed = i, n = n,d = 3,
                       q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
                       sigmavec=sigmavec, sigma_eps=1)
    str(datlist)
    XList <- datlist$XList
    Z <- datlist$Z
    numvarmat <- datlist$numvarmat
    
    
    ## LFR
    XtmpList <- XList
    XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
    XcList <- list(Reduce(cbind, XtmpList))
    try({
      res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 
      timeArray[i,4,jn] <- res_rf$time.use
    }, silent=T)

    
    
  }
  save(timeArray, file=paste0('timeArray_simu4_LFR_npb_changen.rds'))
}



# Fix n while change pvec -------------------------------------------------

cvec <- c(1, 3, 4, 5, 6, 8)
n <- 200
N <- 10
timeArray <- array(dim=c(N, n_methods,  length(cvec)))
colnames(timeArray) <- methodNames
for(jn in c(6:length(cvec))){
  # jn <- 1
  
  pvec2 <- c(50, 50)* cvec[jn] #c(50, 150) # 
  pveclist <- list("gaussian"=pvec2, 'poisson'=pvec2,
                   'binomial' = pvec2)
  message('jn  = ', jn, ', n=', n, ", p=", sum(pvec2)*3)
  for(i in 1:N){
    #i <- 1
    
    datlist <- gendata(pveclist = pveclist, seed = i, n = n,d = 3,
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
    
    library(CMGFM)
    tic <- proc.time()
    rlist <- CMGFM(XList, Z, types=types, q=q, numvarmat, init='LFM', add_IC_iter=F)
    toc <- proc.time()
    time_smgfm <- toc[3] - tic[3]
    timeArray[i,1,jn] <- time_smgfm
    
    ## GFM
    tic <- proc.time()
    res_gfm <- gfm(XList, types=types, q=q, dc_eps = 1e-20)
    toc <- proc.time()
    time_gfm <- toc[3] - tic[3]
    timeArray[i,2,jn] <- time_gfm
    
    
    ### MRRR
    family_use <- list(gaussian(), poisson(), binomial())
    familygroup <- sapply(seq_along(datlist$XList), function(j) rep(j, ncol(datlist$XList[[j]])))
    res_mrrr <- mrrr_run(Reduce(cbind, datlist$XList),  Z = Z, numvarmat, rank0=q, family=family_use,
                         familygroup = unlist(familygroup), maxIter=100) #
    
    timeArray[i,3,jn] <- res_mrrr$time_use
    
    
    # ## LFR
    # XtmpList <- XList
    # XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
    # XcList <- list(Reduce(cbind, XtmpList))
    # try({
    #   res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=500, log.transform=FALSE) # 500
    #   timeArray[i,4,jn] <- res_rf$time.use
    # }, silent=T)
    ## generalizedPCA
    library(generalizedPCA)
    Xall <- cbind(datlist$XList[[2]], datlist$XList[[2]], datlist$XList[[2]])
    t1 <- proc.time()
    res_gPCA <- generalizedPCA_here(Xall, k=q, family = "poisson", quiet = F)
    t2 <- proc.time()
    time_gPCA <- t2[3]-t1[3]
    timeArray[i,5,jn] <- time_gPCA
    
    ### COAP
    tic <- proc.time()
    reslist <- COAP_run(Xall, Z = Z, matrix(colSums(numvarmat),1,2), q= q, rank = ncol(Z)+1, epsELBO = 1e-20)
    toc <- proc.time()
    time_apfactor <- toc[3] - tic[3]
    timeArray[i,6,jn] <- time_apfactor
    
  }
  save(timeArray, file=paste0('timeArray_simu4_changep_exceptLFR_npb.rds'))
}



apply(timeArray, c(2,3) , mean, na.rm=T)
apply(timeArray, c(2,3) , sd, na.rm=T)


### Run LFR
timeArray <- array(dim=c(N, n_methods,  length(cvec)))
colnames(timeArray) <- methodNames
for(jn in 1:length(cvec)){
  # jn <- 3
  
  pvec2 <- c(50, 50)* cvec[jn] #c(50, 150) # 
  pveclist <- list("gaussian"=pvec2, 'poisson'=pvec2,
                   'binomial' = pvec2)
  message('jn  = ', jn, ', n=', n, ", p=", sum(pvec2)*3)
  for(i in 1:N){
    #i <- 2
    
    datlist <- gendata(pveclist = pveclist, seed = i, n = n,d = 3,
                       q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
                       sigmavec=sigmavec, sigma_eps=1)
    str(datlist)
    XList <- datlist$XList
    Z <- datlist$Z
    numvarmat <- datlist$numvarmat
    
    
    ## LFR
    XtmpList <- XList
    XtmpList[[2]] <- log(1+XList[[2]]) # log-transformation
    XcList <- list(Reduce(cbind, XtmpList))
    try({
      res_rf <- MSFR_run(XcList, Z, numvarmat,q,  maxIter=2, log.transform=FALSE,constraint = "block_lower2") # 
      timeArray[i,4,jn] <- res_rf$time.use*250
    }, silent=T)
    
    
    
  }
  save(timeArray, file=paste0('timeArray_simu4_changep_LFR_npb.rds'))
}


