
# Data generation functions -----------------------------------------------
gendata_gauss <- function (seed = 1, n = 300, pvec = c(50, 150), d = 3,
                           q = 6, rho = rep(1,length(pvec)), rho_z=1,
                           sigmavec=rep(0.5,length(pvec)), sigma_eps=1){
  
  # seed = 1; width=20; height=15; pvec = c(50, 150); 
  # q = 6; rho = rep(1, length(pvec)); sigma_eps=1
  ## rho corresponds to the parameters of signal strength for each modality
  if(length(rho) != length(pvec)) stop("gendata_sp_gauss: the length of rho must be equal to that of pvec")
  if(length(sigmavec) != length(pvec)) stop("gendata_sp_gauss: the length of sigmavec must be equal to that of pvec")
  
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  M <- length(pvec)
  p <- sum(pvec)
  set.seed(1)
  betamat <- matrix(rnorm(d*M, sd=2), d, M)
  mu0 <- 0.4 * rnorm(p)
  mu0list <- vec2list(mu0, pvec)
  glist <- list()
  pcum <- c(0, cumsum(pvec))
  
  U0 <- matrix(0, n, M)
  
  B0list <- list()
  Uplist <- list()
  for(m in seq_along(pvec)){
    glist[[m]] <- (pcum[m]+1):pcum[m+1]
    
    Z1 <- matrix(rnorm(pvec[m] * q), pvec[m], q)
    svdZ <- svd(Z1) 
    B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:q])) * rho[m]
    B <- B %*% Diag(sign(B[1,]))
    B0list[[m]] <- B
    Uplist[[m]] <- cbind(mu0list[[m]], B)
    U0[,m] <- rnorm(n, sd=sqrt(sigmavec[m]))
  }
  Umat <- NULL
  betaMat <- NULL
  for(m in 1:M){
    Umat <- cbind(Umat, matrix(U0[,m], n, pvec[m]))
    betaMat <- cbind(betaMat, matrix(betamat[,m], d, pvec[m]))
  }
  
  set.seed(seed)
  F0 <- matrix(rnorm(n*q), n, q)
  Z <- matrix(runif(n*d, -3, 3), n, d)*rho_z
  
  B0 <- Reduce(rbind, B0list)
  X <- Z%*% betaMat + F0 %*% t(B0) + Umat + matrix(mu0, n, p, byrow = T) + mvrnorm(n,rep(0, p), sigma_eps*diag(p))
  XList <- list()
  for(m in seq_along(pvec)){
    XList[[m]] <- X[,glist[[m]]]
  }
  Alist <- list(matrix(0, n, length(pvec)))
  
  return(list(X=X, XList = XList, Z=Z, Alist=Alist, B0list = B0list,
              mu0 = mu0, U0=U0, F0=F0, Uplist=Uplist, beta= betamat,
              sigma_eps=sigma_eps,numvarmat= matrix(pvec, nrow=1)))
}

gendata_pois <- function (seed = 1, n = 300, pvec = c(50, 150), d=3, 
                          q = 6, rho = rep(1,length(pvec)), rho_z=1, 
                          sigmavec=rep(0.5,length(pvec)), sigma_eps=1){
  
  # seed = 1; width=20; height=15; pvec = c(50, 150); alpha0=0.5 
  # q = 6; rho = rep(1, length(pvec)); sigma_eps=1; sigmavec=rep(0.5,length(pvec))
  ## rho corresponds to the parameters of signal strength for each modality
  
  if(length(rho) != length(pvec)) stop("gendata_sp_pois: the length of rho must be equal to that of pvec")
  if(length(sigmavec) != length(pvec)) stop("gendata_sp_pois: the length of sigmavec must be equal to that of pvec")
  
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  M <- length(pvec)
  p <- sum(pvec)
  set.seed(1)
  betamat <- matrix(rnorm(d*M, sd=2), d, M)
  mu0 <- 0.4 * rnorm(p)
  mu0list <- vec2list(mu0, pvec)
  glist <- list()
  pcum <- c(0, cumsum(pvec))
  
  U0 <- matrix(0, n, M)
  
  B0list <- list()
  Uplist <- list()
  for(m in seq_along(pvec)){
    glist[[m]] <- (pcum[m]+1):pcum[m+1]
    
    Z1 <- matrix(rnorm(pvec[m] * q), pvec[m], q)
    Z1 <-  Z1
    svdZ <- svd(Z1) 
    B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:q])) * rho[m]
    B <- B %*% Diag(sign(B[1,]))
    B0list[[m]] <- B
    Uplist[[m]] <- cbind(mu0list[[m]], B)
    U0[,m] <- rnorm(n, sd=sqrt(sigmavec[m]))
  }
  Umat <- NULL
  betaMat <- NULL
  for(m in 1:M){
    Umat <- cbind(Umat, matrix(U0[,m], n, pvec[m]))
    betaMat <- cbind(betaMat, matrix(betamat[,m], d, pvec[m]))
  }
  
  set.seed(seed)
  F0 <- matrix(rnorm(n*q), n, q)
  Z <- matrix(runif(n*d, -3, 3), n, d)*rho_z
  
  B0 <- Reduce(rbind, B0list)
  Mu <- Z%*% betaMat + F0 %*% t(B0) + Umat + matrix(mu0, n, p, byrow = T) + mvrnorm(n,rep(0, p), sigma_eps*diag(p))
  X <- matrix(rpois(n * p, lambda = exp(Mu)), n, p)
  XList <- list()
  for(m in seq_along(pvec)){
    XList[[m]] <- X[,glist[[m]]]
  }
  Alist <- list(matrix(0, n, length(pvec)))
  
  return(list(X=X, XList = XList, Z=Z, Alist=Alist, B0list = B0list,
              mu0 = mu0, U0=U0, F0=F0, Uplist=Uplist, beta= betamat,
              sigma_eps=sigma_eps,numvarmat= matrix(pvec, nrow=1)))
}

gendata_bino <- function (seed = 1, n = 300, pvec = c(50, 150), d=3, 
                          q = 6, rho = rep(1,length(pvec)), rho_z=1, 
                          sigmavec=rep(0.5,length(pvec)), sigma_eps=1, n_bin=1){
  
  # seed = 1; width=20; height=15; pvec = c(50, 150); 
  # q = 6; rho = rep(1, length(pvec)); sigma_eps=1
  ## rho corresponds to the parameters of signal strength for each modality
  
  if(length(rho) != length(pvec)) stop("gendata_sp_pois: the length of rho must be equal to that of pvec")
  if(length(sigmavec) != length(pvec)) stop("gendata_sp_pois: the length of sigmavec must be equal to that of pvec")
  
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  M <- length(pvec)
  p <- sum(pvec)
  set.seed(1)
  betamat <- matrix(rnorm(d*M, sd=2), d, M)
  mu0 <- 0.4 * rnorm(p)
  mu0list <- vec2list(mu0, pvec)
  glist <- list()
  pcum <- c(0, cumsum(pvec))
  
  U0 <- matrix(0, n, M)
  
  B0list <- list()
  Uplist <- list()
  for(m in seq_along(pvec)){
    glist[[m]] <- (pcum[m]+1):pcum[m+1]
    
    Z1 <- matrix(rnorm(pvec[m] * q), pvec[m], q)
    svdZ <- svd(Z1) 
    B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:q])) * rho[m]
    B <- B %*% Diag(sign(B[1,]))
    B0list[[m]] <- B
    Uplist[[m]] <- cbind(mu0list[[m]], B)
    U0[,m] <- rnorm(n, sd=sqrt(sigmavec[m]))
  }
  Umat <- NULL
  betaMat <- NULL
  for(m in 1:M){
    Umat <- cbind(Umat, matrix(U0[,m], n, pvec[m]))
    betaMat <- cbind(betaMat, matrix(betamat[,m], d, pvec[m]))
  }
  
  set.seed(seed)
  F0 <- matrix(rnorm(n*q), n, q)
  Z <- matrix(runif(n*d, -3, 3), n, d)*rho_z
  
  B0 <- Reduce(rbind, B0list)
  Mu <- Z%*% betaMat + F0 %*% t(B0) + Umat + matrix(mu0, n, p, byrow = T) + mvrnorm(n,rep(0, p), sigma_eps*diag(p))
  
  X <- matrix(rbinom(n * p, size=n_bin,prob  = 1/(1+exp(-Mu))), n, p)
  XList <- list()
  for(m in seq_along(pvec)){
    XList[[m]] <- X[,glist[[m]]]
  }
  Alist <- list(matrix(0, n, length(pvec)))
  
  return(list(X=X, XList = XList, Z=Z, Alist=Alist, B0list = B0list,
              mu0 = mu0, U0=U0, F0=F0, Uplist=Uplist, beta= betamat,
              sigma_eps=sigma_eps,numvarmat= matrix(pvec, nrow=1)))
}



gendata <- function (seed = 1, n=300, pveclist = list('gaussian'=c(50, 150),'poisson'=c(50),
                                                      'binomial'=c(100,60)),
                     q = 6,  d= 3, rho = rep(1,length(pveclist)), rho_z=1,
                     sigmavec=rep(0.5, length(pveclist)), n_bin=1, sigma_eps=1){
  
  # seed = 1; width=20; height=15; pvec = c(50, 150); 
  # q = 6; rho = rep(1, length(pvec)); sigma_eps=1
  ## rho corresponds to the parameters of signal strength for each modality
  
  if(length(rho) != length(pveclist)) stop("gendata_sp: the length of rho must be equal to that of pveclist")
  if(length(sigmavec) != length(pveclist)) stop("gendata_sp: the length of sigmavec must be equal to that of pveclist")
  
  library(MASS)
  library(LaplacesDemon)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  names_pvec <- names(pveclist)
  if(!all(names_pvec %in% c("gaussian", "poisson", "binomial"))) stop("gendata_sp: the names of pveclist must be in ('n', 'p', 'b')!")
  
  n_vars <- sapply(pveclist, length)
  numvarmat <- matrix(0, nrow=length(n_vars), ncol=max(n_vars))
  for(id in seq_along(n_vars)){
    numvarmat[id, 1:n_vars[id]] <- pveclist[[id]]
  }
  row.names(numvarmat) <- names_pvec
  
  
  set.seed(seed)
  F0 <- matrix(rnorm(n*q), n, q)
  Z <- matrix(runif(n*d, -3, 3), n, d)*rho_z
  
  ## generate mu0 and B0
  p_all <- sum(unlist(pveclist))
  mu0List <- list()
  XList <- list()
  U0List <- list()
  B0List <- list()
  beta0List <- list()
  UpList <- list()
  set.seed(1)
  for(di in 1:length(pveclist)){
    
    pvec <- pveclist[[di]]
    p <- sum(pvec)
    M <- length(pvec)
    betamat <- matrix(rnorm(d*M, sd=2), d, M)
    beta0List[[di]] <- betamat
    mu0 <- 0.4 * rnorm(p)
    mu0list <- vec2list(mu0, pvec)
    mu0List[[di]] <- mu0
    glist <- list()
    pcum <- c(0, cumsum(pvec))
    
    B0list <- list()
    Uplist <- list()
    for(m in seq_along(pvec)){
      glist[[m]] <- (pcum[m]+1):pcum[m+1]
      
      Z1 <- matrix(rnorm(pvec[m] * q), pvec[m], q)
      svdZ <- svd(Z1) 
      B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:q])) * rho[di] # control the signal for each type of variables 
      B <- B %*% Diag(sign(B[1,]))
      B0list[[m]] <- B
      Uplist[[m]] <- cbind(mu0list[[m]], B)
      
    }
    B0List[[di]] <- B0list
    UpList[[di]] <- Uplist
    set.seed(seed)
    U0 <- matrix(0, n, M)
    for(m in seq_along(pvec)){
      U0[,m] <- rnorm(n, sd=sqrt(sigmavec[di]))
    }
    U0List[[di]] <- U0
    Umat <- NULL
    betaMat <- NULL
    for(m in 1:M){
      Umat <- cbind(Umat, matrix(U0[,m], n, pvec[m]))
      betaMat <- cbind(betaMat, matrix(betamat[,m], d, pvec[m]))
    }
    
    B0 <- Reduce(rbind, B0list)
    set.seed(seed)
    X <- Z %*% betaMat + F0 %*% t(B0) + Umat + matrix(mu0, n, p, byrow = T) + mvrnorm(n,rep(0, p), sigma_eps*diag(p))
    if(names_pvec[di] == 'poisson'){
      X <- matrix(rpois(n * p, lambda = exp(X)), n, p)
    }
    if(names_pvec[di] == 'binomial'){
      
      X <- matrix(rbinom(n * p, size=n_bin,prob  = 1/(1+exp(-X))), n, p)
    }
    
    XList[[di]] <- X
  }
  
  Alist <- list(matrix(0, n, sum(numvarmat>0)))
  
  return(list(XList = XList, Z=Z, types=names_pvec, Alist=Alist, numvarmat=numvarmat, 
              B0List = B0List, mu0List = mu0List, beta0List=beta0List,Uplist = UpList,  U0List=U0List,
              F0=F0, sigma_eps=sigma_eps))
}



cor.mat <- function(p, rho, type='toeplitz'){
  
  mat <- diag(p)
  if(type=='toeplitz'){
    for(i in 2:p){
      for(j in 1:i){
        mat[i,j] <- mat[j,i] <- rho^(abs(i-j))
      }
    }
  }
  if(type=='identity'){
    mat[mat==0] <- rho
  }
  return(mat)
}
cov.mat <- function(sdvec, rho, type='toeplitz'){
  p <- length(sdvec)
  cormat <- cor.mat(p, rho, type)
  covmat <- matrix(0, p, p)
  for(i in 1:p){
    for(j in 1:p){
      covmat[i,j] <- sdvec[i]*sdvec[j]*cormat[i,j]
    }
  }
  return(covmat)
}




# Unified function -------------------------------------------------------
library(GFM)
Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}
colSD <- function(mat, na.rm=TRUE){
  apply(mat, 2, sd, na.rm=na.rm)
}
normvec <- function(vec, norm=c('L1', 'L2')){
  
  norm <- match.arg(norm)
  if(norm=='L1'){
    val = mean(abs(vec))
  }
  if(norm=='L2'){
    val = mean(sum(vec^2))
  }
  return(val)
}

meanTr <- function(hBlist, Blist,  type='trace_statistic'){
  ###It is noted that the trace statistics is not symmetric, the true value must be in last
  trvec <- sapply(1:length(Blist), function(j) measurefun(hBlist[[j]], Blist[[j]], type = type))
  return(mean(trvec))
}
allTr <- function(hBlist, Blist,  type='trace_statistic'){
  trvec <- sapply(1:length(Blist), function(j) measurefun(Blist[[j]], hBlist[[j]], type = type))
  return(trvec)
}
mat2list <- function(B, pvec, by_row=TRUE){
  Blist <- list()
  pcum = c(0, cumsum(pvec))
  for(i in 1:length(pvec)){
    if(by_row){
      Blist[[i]] <- B[(pcum[i]+1):pcum[i+1],]
    }else{
      Blist[[i]] <- B[, (pcum[i]+1):pcum[i+1]]
    }
  }
  return(Blist)
}
vec2list <- function(y_int, nvec){
  
  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")
  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){
    
    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}

# Compared methods --------------------------------------------------------

## Compare with MRRR
mrrr_run <- function(Y, Z, numvarmat, rank0,family=list(poisson()),
                     familygroup,  epsilon = 1e-4, sv.tol = 1e-2,
                     lambdaSVD=0.1, maxIter = 2000, trace=TRUE, trunc=500){
  # epsilon = 1e-4; sv.tol = 1e-2; maxIter = 30; trace=TRUE,lambdaSVD=0.1
  
  require(rrpack)
  q <- rank0
  n <- nrow(Y); p <- ncol(Y)
  X <- cbind(cbind(1, Z),  diag(n))
  d <- ncol(Z)
  ## Trunction
  Y[Y>trunc] <- trunc
  Y[Y< -trunc] <- -trunc
  
  tic <- proc.time()
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  pcums <- cumsum(pvec)
  idxlist <- list()
  idxlist[[1]] <- 1:pcums[1]
  if(length(pvec)>1){
    for(i in 2:length(pvec)){
      idxlist[[i]] <- (pcums[i-1]+1):pcums[i] 
    }
  }
  
  svdX0d1 <- svd(X)$d[1]
  init1 = list(kappaC0 = svdX0d1 * 5)
  offset = NULL
  control = list(epsilon = epsilon, sv.tol = sv.tol, maxit = maxIter,
                 trace = trace, gammaC0 = 1.1, plot.cv = TRUE,
                 conv.obj = TRUE)
  res_mrrr <- mrrr(Y=Y, X=X[,-1], family = family, familygroup = familygroup,
                   penstr = list(penaltySVD = "rankCon", lambdaSVD = lambdaSVD),
                   control = control, init = init1, maxrank = rank0+d) #
  
  hmu <- res_mrrr$coef[1,]
  hbeta <- t(res_mrrr$coef[2:(d+1),])
  hTheta <- res_mrrr$coef[-c(1:(d+1)),]
  
  hbeta_rf <- NULL
  for(i in seq_along(pvec)){
    hbeta_rf <- cbind(hbeta_rf, colMeans(hbeta[idxlist[[i]],]))
  }
  
  # Matrix::rankMatrix(hTheta)
  svd_Theta <- svd(hTheta, nu=q, nv=q)
  hH <- svd_Theta$u
  hB <- svd_Theta$v %*% Diag(svd_Theta$d[1:q])
  toc <- proc.time()
  time_mrrr <- toc[3] - tic[3]
  
  return(list(hH=hH, hB=hB, hmu= hmu, beta=hbeta_rf, time_use=time_mrrr))
}


## Factor regression
MSFR_run <- function(XcList, Z, numvarmat, q, qs=1, maxIter=1e4, log.transform=FALSE,constraint = "block_lower1", dir_source=NULL){
  
  require(MSFA)
  require(psych)
  if(is.null(dir_source)){
    source("MSFR_main_R_MSFR_V1.R")
  }else{
    source(paste0(dir_source, "/MSFR_main_R_MSFR_V1.R"))
  }
  
  #fa <- psych::fa
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  pcums <- cumsum(pvec)
  idxlist <- list()
  idxlist[[1]] <- 1:pcums[1]
  if(length(pvec)>1){
    for(i in 2:length(pvec)){
        idxlist[[i]] <- (pcums[i-1]+1):pcums[i] 
    }
  }
  ZList <- list(cbind(1,Z))
  
  B_s <- ZList
  if(log.transform){
    X_s <- lapply(XcList, function(x) log(1+x))
  }else{
    X_s <- XcList #
  }
  t1<- proc.time()
  test <- start_msfa(X_s, B_s, 5, k=q, j_s=qs, constraint = constraint, method = "adhoc")
  obj <- try({
    EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE, constraint = constraint, nIt=maxIter)
  }, silent=TRUE)
  t2<- proc.time()
  if(class(obj) == "try-error"){
    EM_beta <- test
  }else{
    EM_beta$F <- t(EM_beta$E_f[[1]])
    EM_beta$H <- t(EM_beta$E_l[[1]])
    
    hbeta_rf <- NULL
    for(i in seq_along(pvec)){
      hbeta_rf <- cbind(hbeta_rf, colMeans(EM_beta$beta[idxlist[[i]],-1]))
    }
    EM_beta$Upsilon <- cbind(EM_beta$beta[,1], EM_beta$Phi)
    hUplist <- list()
    for(i in seq_along(pvec)){
      hUplist[[i]] <- EM_beta$Upsilon[idxlist[[i]],]
    }
    EM_beta$hUplist <- hUplist
    EM_beta$beta_rf <- hbeta_rf
    EM_beta$lambdavec <- EM_beta$psi_s[[1]]
  }
  
  
  EM_beta$time.use <- t2[3] - t1[3]
  return(EM_beta)
}

COAP_run <- function(X, Z, numvarmat, q, rank,  epsELBO = 1e-7, ...){
  
  
  library(COAP)
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  pcums <- cumsum(pvec)
  idxlist <- list()
  idxlist[[1]] <- 1:pcums[1]
  if(length(pvec)>1){
    for(i in 2:length(pvec)){
      idxlist[[i]] <- (pcums[i-1]+1):pcums[i] 
    }
  }
  
  tic <- proc.time()
  reslist <- COAP::RR_COAP(X, Z = cbind(1,Z), q= q, rank_use = rank, ...)
  toc <- proc.time()
  time_apfactor <- toc[3] - tic[3]
  reslist$time.use <- time_apfactor
  reslist$hmu <- reslist$bbeta[,1]
  
  hbeta_rf <- NULL
  for(i in seq_along(pvec)){
    hbeta_rf <- cbind(hbeta_rf, colMeans(reslist$bbeta[idxlist[[i]],-1]))
  }
  reslist$hbeta <- hbeta_rf
  
  return(reslist)
  
}



generalizedPCA_here <- function (x, k = 2, M = 4, family = c("gaussian", "binomial", 
                                                             "poisson", "multinomial"), weights, quiet = TRUE, 
                                 majorizer = c("row", "all"), partial_decomp = FALSE, 
                                 max_iters = 1000, conv_criteria = 1e-05, random_start = FALSE, 
                                 start_U, start_mu, main_effects = TRUE, normalize = FALSE, 
                                 validation, val_weights){
  # x <- tX
  # k = 2; M = 4; family =  "poisson";  quiet = TRUE; 
  #                          majorizer = "row"; partial_decomp = FALSE; 
  #                          max_iters = 1000; conv_criteria = 1e-05; random_start = FALSE; 
  #                          main_effects = TRUE; normalize = FALSE
  #                          start_U; start_mu; 
  
  family = match.arg(family)
  # generalizedPCA:::check_family(x, family)
  majorizer = match.arg(majorizer)
  x = as.matrix(x)
  missing_mat = is.na(x)
  n = nrow(x)
  d = ncol(x)
  ones = rep(1, n)
  if (missing(weights)) {
    weights = 1
    sum_weights = sum(!is.na(x))
  }
  # else {
  #   weights[is.na(x)] <- 0
  #   if (any(is.na(weights))) {
  #     stop("Can't have NA in weights")
  #   }
  #   if (any(weights < 0)) {
  #     stop("weights must be non-negative")
  #   }
  #   if (!all(dim(weights) == dim(x))) {
  #     stop("x and weights are not the same dimension")
  #   }
  #   sum_weights = sum(weights)
  # }
  saturated_natural_parameters <- generalizedPCA:::saturated_natural_parameters
  if (main_effects) {
    if (length(weights) == 1) {
      weighted_col_means = colMeans(x, na.rm = TRUE)
    }
    else {
      weighted_col_means = colSums(x * weights, na.rm = TRUE)/colSums(weights)
    }
    null_theta = as.numeric(saturated_natural_parameters(matrix(weighted_col_means, 
                                                                1), family, M))
  }else {
    null_theta = rep(0, d)
  }
  if (normalize) {
    if (any(apply(x, 2, stats::var, na.rm = TRUE) == 0)) {
      stop("At least one variable has variance of 0. Cannot normalize")
    }
    eta_sat_nat = saturated_natural_parameters(x, family, 
                                               M = Inf)
    norms = sapply(1:d, function(j) {
      2 * (exp_fam_log_like(x[, j], eta_sat_nat[, j], family, 
                            weights) - exp_fam_log_like(x[, j], rep(null_theta[j], 
                                                                    n), family, weights))
    })/n
    if (any(norms <= 0)) {
      stop("Normalization caused weights to be <= 0")
    }
    if (length(weights) == 1) {
      weights = outer(ones, 1/norms)
    }
    else {
      weights = sweep(weights, 2, 1/norms, "*")
    }
  }
  if (M == 0) {
    if (any(is.na(x))) {
      stop("Cannot solve for M with missing weights")
    }
    M = 4
    solve_M = TRUE
    if (!missing(validation)) {
      if (ncol(validation) != ncol(x)) {
        stop("validation does not have the same variables as x")
      }
      if (missing(val_weights)) {
        val_weights = 1
      }
      else {
        if (!all(dim(val_weights) == dim(validation))) {
          stop("validation and val_weights are not the same dimension")
        }
      }
      validation = as.matrix(validation)
      M_mat = exp_fam_sat_ind_mat(validation, family)
    }
    else {
      M_mat = exp_fam_sat_ind_mat(x, family)
    }
  }else {
    solve_M = FALSE
  }
  if (family == "gaussian" & all(weights == 1) & sum(missing_mat) == 
      0) {
    max_iters = 0
  }
  
  if (main_effects) {
    # if (!missing(start_mu)) {
    #   mu = start_mu
    # }else {
    eta = generalizedPCA:::saturated_natural_parameters(x, family, M = M)
    is.na(eta[is.na(x)]) <- TRUE
    mu = colMeans(eta, na.rm = TRUE)
    # }
  }else {
    mu = rep(0, d)
  }
  eta = saturated_natural_parameters(x, family, M = M) + missing_mat * 
    outer(ones, mu)
  eta_centered = scale(eta, center = mu, scale = FALSE)
  if (!missing(start_U)) {
    U = sweep(start_U, 2, sqrt(colSums(start_U^2)), "/")
  }else if (random_start) {
    U = matrix(rnorm(d * k), d, k)
    U = qr.Q(qr(U))
  }else {
    if (partial_decomp) {
      udv = RSpectra::svds(scale(eta, center = mu, scale = normalize),
                           k, nu = k, nv = k)
    }
    else {
      udv = svd(scale(eta, center = mu, scale = normalize))
    }
    U = matrix(udv$v[, 1:k], d, k)
  }
  exp_fam_log_like <- generalizedPCA:::exp_fam_log_like
  eta_sat_nat = saturated_natural_parameters(x, family, M = Inf)
  sat_loglike = exp_fam_log_like(x, eta_sat_nat, family, weights)
  loss_trace = numeric(max_iters + 1)
  theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
  loglike <- exp_fam_log_like(x, theta, family, weights)
  loss_trace[1] = 2 * (sat_loglike - loglike)/sum_weights
  ptm <- proc.time()
  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }
  for (m in seq_len(max_iters)) {
    # m <- 1
    last_U = U
    last_M = M
    last_mu = mu
    if (solve_M) {
      gpca_obj = structure(list(mu = mu, U = U, M = M, 
                                family = family), class = "gpca")
      if (missing(validation)) {
        fitted_theta = predict(gpca_obj, newdata = x, 
                               type = "link")
      }
      else {
        fitted_theta = predict(gpca_obj, newdata = validation, 
                               type = "link")
      }
      fitted_mean = exp_fam_mean(fitted_theta, family)
      if (missing(validation)) {
        M_slope = sum(((fitted_mean - x) * weights * 
                         (M_mat %*% tcrossprod(U)))[!is.na(M_mat)])
        fitted_var = exp_fam_variance(fitted_theta, family, 
                                      weights)
      }
      else {
        M_slope = sum(((fitted_mean - validation) * val_weights * 
                         (M_mat %*% tcrossprod(U)))[!is.na(M_mat)])
        fitted_var = exp_fam_variance(fitted_theta, family, 
                                      val_weights)
      }
      M_curve = sum((fitted_var * (M_mat %*% tcrossprod(U))^2)[!is.na(M_mat)])
      M = max(M - M_slope/M_curve, 0)
      eta = saturated_natural_parameters(x, family, M = M) + 
        missing_mat * outer(ones, mu)
      eta_centered = scale(eta, center = mu, scale = FALSE)
      theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
    }
    exp_fam_mean <- generalizedPCA:::exp_fam_mean
    exp_fam_variance <- generalizedPCA:::exp_fam_variance
    first_dir = exp_fam_mean(theta, family)
    second_dir = exp_fam_variance(theta, family, weights)
    if (majorizer == "row") {
      W = apply(second_dir, 1, max)
    } else if (majorizer == "all") {
      W = rep(max(second_dir), n)
    }
    Z = as.matrix(theta + weights * (x - first_dir)/outer(W, 
                                                          rep(1, d)))
    Z[is.na(x)] <- theta[is.na(x)]
    if (main_effects) {
      mu = as.numeric(colSums((Z - eta %*% tcrossprod(U)) * 
                                W)/sum(W))
    }
    eta = saturated_natural_parameters(x, family, M = M) + 
      missing_mat * outer(ones, mu)
    eta_centered = scale(eta, center = mu, scale = FALSE)
    mat_temp = t(eta_centered * W) %*% scale(Z, center = mu, 
                                             scale = FALSE)
    mat_temp = mat_temp + t(mat_temp) - t(eta_centered * 
                                            W) %*% eta_centered
    repeat {
      if (partial_decomp) {
        eig = RSpectra::eigs_sym(mat_temp, k = min(k + 
                                                     2, d))
      }
      else {
        eig = eigen(mat_temp, symmetric = TRUE)
      }
      U = matrix(eig$vectors[, 1:k], d, k)
      theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
      this_loglike <- exp_fam_log_like(x, theta, family, 
                                       weights)
      if (!partial_decomp | this_loglike >= loglike) {
        loglike = this_loglike
        break
      }
      else {
        partial_decomp = FALSE
        if (!quiet) {
          cat("RSpectra::eigs_sym was too inaccurate in iteration ", 
              m, ". Switched to base::eigen")
        }
      }
    }
    loss_trace[m + 1] = 2 * (sat_loglike - loglike)/sum_weights
    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters/m * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(m, "  ", loss_trace[m + 1], "")
      cat(round(time_elapsed/3600, 1), "hours elapsed. Max", 
          round(time_remain/3600, 1), "hours remain.\n")
    }
    if (m > 4) {
      if (abs(loss_trace[m] - loss_trace[m + 1]) < 2 * 
          conv_criteria) {
        break
      }
    }
  }
  if (max_iters > 0 && (loss_trace[m + 1] - loss_trace[m]) > 
      (1e-10)) {
    U = last_U
    mu = last_mu
    M = last_M
    m = m - 1
    eta = saturated_natural_parameters(x, family, M = M) + 
      missing_mat * outer(ones, mu)
    eta_centered = scale(eta, center = mu, scale = FALSE)
    if (family != "poisson") {
      warning("Deviance increased in last iteration.\nThis should not happen!")
    }
    else {
      message("Deviance increased in last iteration.")
    }
  } else if (max_iters == 0) {
    m = 0
  }
  null_loglike = exp_fam_log_like(x, outer(ones, null_theta), 
                                  family, weights)
  object <- list(mu = mu, U = U, PCs = eta_centered %*% U, M = M, family = family,
                 iters = m, loss_trace = loss_trace[1:(m + 1)], 
                 prop_deviance_expl = 1 - (loglike - sat_loglike)/(null_loglike - sat_loglike))
  
  class(object) <- "gpca"
  object
}


