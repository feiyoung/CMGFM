Diag <- function (vec){
  q <- length(vec)
  if (q > 1) {
    y <- diag(vec)
  }
  else {
    y <- matrix(vec, 1, 1)
  }
  return(y)
}

cor.mat <- function (p, rho, type = "toeplitz"){
  if (p == 1) 
    return(matrix(1, 1, 1))
  mat <- diag(p)
  if (type == "toeplitz") {
    for (i in 2:p) {
      for (j in 1:i) {
        mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
      }
    }
  }
  if (type == "identity") {
    mat[mat == 0] <- rho
  }
  return(mat)
}


#' Generate simulated data
#' @description Generate simulated data from covariate-augumented generalized factor model
#' @param seed a positive integer, the random seed for reproducibility of data generation process.
#' @param n a positive integer, specify the sample size. 
#' @param pveclist a named list, specify the number of modalities for each variable type and dimension of variables in each modality.
#' @param q a positive integer,  specify the number of modality-shared factors.
#' @param d a positive integer,  specify the dimension of covariate matrix.
#' @param rho a numeric vector with length \code{length(pveclist)} and positive elements, specify the signal strength of loading matrix for each modality with the same variable type. 
#' @param rho_z a positive real, specify the signal strength of covariates.
#' @param sigmavec a positive vector with length \code{length(pveclist)}, the variance of modality-specified latent factors.
#' @param n_bin a positive integer, specify the number of trails in Binomial distribution.
#' @param sigma_eps a positive real, the variance of overdispersion error.
#' @param seed.para a positive integer, the random seed for reproducibility of data generation process by fixing the regression coefficient vector and loading matrices.
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{XList} - a list consisting of multiple matrices in which each matrix has the same type of values, i.e., continuous, or count, or binomial/binary values.
#'   \item \code{Z} - a matrix, the fixed-dimensional covariate matrix with control variables;
#'   \item \code{Alist} - the the offset vector for each modality; 
#'   \item \code{B0list} - the true loading matrix for each modality; 
#'   \item \code{mu0} - the true intercept vector for each modality;
#'   \item \code{U0} - the  modality-specified factor vector;
#'   \item \code{F0} -  the  modality-shared factor matrix; 
#'   \item \code{Uplist} - the true intercept-loading matrix for each modality;
#'   \item \code{beta} - the true regression coefficient vector for each modality; 
#'   \item \code{sigma_eps} -  the standard deviation of error term;
#'   \item \code{numvarmat} - a length(types)-by-d matrix, the number of variables in modalities that belong to the same type.
#' }
#' @details None
#' @seealso \code{\link{CMGFM}}
#' @references None
#' @export
#' @importFrom  MASS mvrnorm
#' @importFrom  stats rnorm rpois rbinom runif
#'
#' @examples
#' n <- 300; 
#' pveclist = list('gaussian'=c(50, 150),'poisson'=c(50),'binomial'=c(100,60))
#' d <- 20; q <- 6;
#' datlist <- gendata_cmgfm(n=n, pveclist=pveclist, q=q, d=d)
#' str(datlist)

gendata_cmgfm <- function (seed = 1, n=300, pveclist = list('gaussian'=c(50, 150),'poisson'=c(50),
                                                                       'binomial'=c(100,60)),
                                      q = 6,  d= 3, rho = rep(1,length(pveclist)), rho_z=1,
                                      sigmavec=rep(0.5, length(pveclist)), n_bin=1, sigma_eps=1, seed.para=1){
  
  # seed = 1; width=20; height=15; pvec = c(50, 150); 
  # q = 6; rho = rep(1, length(pvec)); sigma_eps=1
  ## rho corresponds to the parameters of signal strength for each modality
  
  if(length(rho) != length(pveclist)) stop("gendata_cmgfm: the length of rho must be equal to that of pveclist")
  if(length(sigmavec) != length(pveclist)) stop("gendata_cmgfm: the length of sigmavec must be equal to that of pveclist")
  

  names_pvec <- names(pveclist)
  if(!all(names_pvec %in% c("gaussian", "poisson", "binomial"))) stop("gendata_cmgfm: the names of pveclist must be in ('gaussian', 'poisson', 'binomial')!")
  
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
  set.seed(seed.para)
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

