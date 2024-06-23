# generate man files
# devtools::document()
# R CMD check --as-cran CMGFM_1.1.tar.gz
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("simu")
# pkgdown::build_article("mouseSpleen")

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('html_document'))

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('pdf_document'), clean = F)


get_initials <- function(X, q){
  library(irlba)
  n <- nrow(X); p <- ncol(X)
  mu <- colMeans(X)
  X <- X - matrix(mu, nrow=n, ncol=p, byrow=TRUE)
  svdX  <- irlba(A =X, nv = q)
  PCs <- sqrt(n) * svdX$u
  loadings <- svdX$v %*% diag(svdX$d[1:q]) / sqrt(n)
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(hH = PCs, hB = loadings, hmu=mu,sigma2vec = Lam_vec))
  
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






#' Fit the CMGFM model
#' @description Fit the covariate-augumented generalized factor model
#' @param XList a list consisting of multiple matrices in which each matrix has the same type of values, i.e., continuous, or count, or binomial/binary values.
#' @param Z a matrix, the fixed-dimensional covariate matrix with control variables.
#' @param types a string vector, specify the variable type in each matrix in \code{XList};
#' @param numvarmat a \code{length(types)}-by-d matrix, specify the number of variables in modalities that belong to the same type.
#' @param Alist an optional vector, the offset for each unit; default as full-zero vector.
#' @param init an optional character, specify the method in initialization.
#' @param q an optional string, specify the number of factors; default as 15.
#' @param epsELBO  an optional positive value, tolerance of relative variation rate of the evidence lower bound value, default as '1e-8'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param add_IC_iter a logical value, add the identifiability condition in iterative algorithm or add it after algorithm converges; default as FALSE.
#' @param seed an integer, set the random seed in initialization, default as 1;
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{betaf} - the estimated regression coefficient vector for each modality; 
#'   \item \code{Bf} - the estimated loading matrix for each modality; 
#'   \item \code{M} - the estimated modality-shared factor matrix; 
#'   \item \code{Xif} - the estimated modality-specified factor vector;
#'   \item \code{S} - the estimated covariance matrix of modality-shared latent factors;
#'   \item \code{Om} - the posterior variance of modality-specified latent factors;
#'   \item \code{muf} - the estimated intercept vector for each modality;
#'   \item \code{Sigmam} - the variance of modality-specified factors;
#'   \item \code{invLambdaf} - the inverse of the estimated variances of error for each modality.
#'   \item \code{ELBO} -  the ELBO value when algorithm stops;
#'   \item \code{ELBO_seq} - the sequence of ELBO values.
#'   \item \code{time_use} - the running time in model fitting;
#' }
#' @details None
#' @seealso None
#' @references None
#' @export
#' @useDynLib CMGFM, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#' @importFrom  GFM gfm
#'
#'
#' @examples
#' pveclist <- list('gaussian'=c(50, 150),'poisson'=c(50, 150),
#'    'binomial'=c(100,60))
#' q <- 6
#' sigmavec <- rep(1,3)
#' pvec <- unlist(pveclist)
#' datlist <- gendata_cmgfm(pveclist = pveclist, seed = 1, n = 300,d = 3,
#'                          q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
#'                          sigmavec=sigmavec, sigma_eps=1)
#' XList <- datlist$XList
#' Z <- datlist$Z
#' numvarmat <- datlist$numvarmat
#' types <- datlist$types
#' rlist <- CMGFM(XList, Z, types=types, numvarmat, q=q)
#' str(rlist)
#' 
CMGFM <- function(XList, Z, types, numvarmat, q=15, Alist=NULL, init=c("LFM", "GFM","random"),
                  maxIter=30, epsELBO=1e-8,verbose=TRUE, add_IC_iter=FALSE, seed=1){
  
  
  ### check the arguments:
  if(!inherits(XList, "list")) stop("CMGFM: XList must be a list!")
  if(!inherits(XList[[1]], "matrix")) stop("CMGFM: Each component of XList must be a matrix!")
  if(q<1 ) stop("CMGFM: q must be other positive integer!")
  if(!inherits(numvarmat, "matrix"))  stop("CMGFM: numvarmat must be a matrix!")
  if((!is.null(Alist)) && (!inherits(Alist, "list"))) stop("CMGFM: Alist must be NULL or a list!")
  if(any(sapply(XList, function(x) sum(is.na(x))))) stop("CMGFM:  Each component of XList can not include missing values!")
  ####
  
  init <- match.arg(init)
  n <- nrow(XList[[1]]);p_min <- min(sapply(XList, ncol))
  if(p_min<2) stop("CMGFM: Number of variables in each modality must be at least no less than 2!")
  if(q>= p_min) stop("CMGFM: number of variables in each modality must be greater than q")
  if(n <2) stop("CMGFM: number of observations must be at least no less than 2!")
  ## Initialize for the algorithm
  if(length(XList) != nrow(numvarmat)) stop("CMGFM: the length of XList must be equal to the nrow of numvarmat!")
  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types]) ## match IDs
  n <- nrow(XList[[1]])
  
  ## Initialize for the algorithm
  Slist_y_int <- invLambdalist_int <- Blist_int<-mulist_int <- list()
  Sigmamlist_int <- Mulist_y_int <- betalist_int<-list()
  
  
  # list(matrix(0, nrow=n, ncol=2), matrix(0, n,2), matrix(0, n,1))
  d <- ncol(Z)
  pvec <- rowSums(numvarmat)
  for(id in seq_along(XList)){
    pm <- pvec[id]
    if(typeID[id]==2){
      Mulist_y_int[[id]] <-log(XList[[id]]+1)
    }else{
      Mulist_y_int[[id]] <- XList[[id]]
      
    }
    Slist_y_int[[id]] <- matrix(0, nrow=n, ncol=pm)
    
    invLambdalist_int[[id]] <- rep(1, pm)
    m1 <- sum(numvarmat[id,]>0)
    Sigmamlist_int[[id]] <- matrix(1, 1, m1)
    if(is.null(Alist) && id ==1){
      Alist <- list()
    }
    Alist[[id]] <- matrix(0, nrow=n, ncol=m1)
    betalist_int[[id]] <- matrix(0, nrow=d, ncol=m1)
  }
  
  
  if(init == 'random'){
    set.seed(seed)
    for(id in seq_along(XList)){
      pm <- pvec[id]
      Blist_int[[id]] <- matrix(1, pm, q)
      mulist_int[[id]] <- rep(0, pm)
      M_int <- matrix(0, n, q)
    }
  }else if(init == 'LFM'){
    Fac <- get_initials(Reduce(cbind, Mulist_y_int), q= q)
    M_int <- Fac$hH 
    Blist_int <- mat2list(Fac$hB, pvec)
    mulist_int <- vec2list(Fac$hmu, pvec)
  }else if(init == "GFM"){
    Fac <- gfm(XList, types= types, q=q, verbose = FALSE, maxIter = 5)
    M_int <- Fac$hH
    Blist_int <- mat2list(Fac$hB, pvec)
    mulist_int <- vec2list(Fac$hmu, pvec)
  }
  
  M_int <- matrix(0, n, q)
  S_int <- diag(rep(1,q))
  Xi_int <- 0.5*M_int[,1]
  O_int <- 1
  
  tic <- proc.time()
  rlist <- vb_cmgfmcpp(XList, typeID, numvarmat, Alist, Z, Mulist_y_int, Slist_y_int, 
                       Sigmamlist_int, invLambdalist_int, Blist_int, mulist_int, 
                       betalist_int, M_int, S_int, Xi_int, O_int, epsELBO, maxIter, 
                       verbose, add_IC_inter = add_IC_iter, update_sigma = TRUE) 
  toc <- proc.time()
  
  ## Reorganize the results
  ind <- sapply(rlist$Bf, function(x) dim(x)[1]>1)
  rlist$numvarmat <- numvarmat
  rlist$Bf <- rlist$Bf[ind]
  rlist$betaf <-  rlist$betaf[ind]
  rlist$Xif <- rlist$Xif[ind]
  rlist$muf <- rlist$muf[ind]
  rlist$invLambdaf <- rlist$invLambdaf[ind]
  rlist$ELBO_seq <- rlist$ELBO_seq[-1]
  rlist$time_use <- toc[3] - tic[3]
  
  return(rlist)
}


### Select the number of factors
#' Select the number of factors
#' @description Select the number of factors using maximum singular value ratio based method
#' @param XList a list consisting of multiple matrices in which each matrix has the same type of values, i.e., continuous, or count, or binomial/binary values.
#' @param Z a matrix, the fixed-dimensional covariate matrix with control variables.
#' @param types a string vector, specify the variable type in each matrix in \code{XList};
#' @param numvarmat a \code{length(types)}-by-d matrix, specify the number of variables in modalities that belong to the same type.
#' @param Alist an optional vector, the offset for each unit; default as full-zero vector.
#' @param q_max an optional string, specify the maximum number of factors; default as 20.
#' @param threshold  an optional positive value, a cutoff to filter the singular values that are smaller than it.
#' @param ... other arguments passed to CMGFM
#' @return return the estimated number of factors.
#' @details None
#' @seealso None
#' @references None
#' @export
#' @useDynLib CMGFM, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#' @importFrom  GFM gfm
#'
#'
#' @examples
#' pveclist <- list('gaussian'=c(50, 150),'poisson'=c(50, 150),
#'    'binomial'=c(100,60))
#' q <- 6
#' sigmavec <- rep(1,3)
#' pvec <- unlist(pveclist)
#' datlist <- gendata_cmgfm(pveclist = pveclist, seed = 1, n = 300,d = 3,
#'                          q = q, rho = rep(1,length(pveclist)), rho_z=0.2,
#'                          sigmavec=sigmavec, sigma_eps=1)
#' XList <- datlist$XList
#' Z <- datlist$Z
#' numvarmat <- datlist$numvarmat
#' types <- datlist$types
#' hq <- MSVR(XList, Z, types=types, numvarmat, q_max=20)
#' 
#' print(c(q_true=q, q_est=hq))

MSVR <- function(XList, Z, types, numvarmat, Alist=NULL, q_max=20, threshold=1e-5, ...){
  rlist <- CMGFM(XList, Z, types, q=q_max, numvarmat, Alist=Alist,...)
  getratio <- function(B, threshold){
    B_svalues <- svd(B)$d
    B_svalues <- B_svalues[B_svalues>threshold]
    ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
    nn <- length(ratio_fac)
    return(which.max(ratio_fac[-nn]))
  }
  qvec <- sapply(rlist$Bf, getratio, threshold=threshold)
  return(max(qvec))
}
