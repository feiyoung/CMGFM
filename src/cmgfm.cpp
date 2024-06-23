// This script implement spatial multimodality linear factor model based on variational inference.
// Date: 2023-02-11

// Revised log:


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

// #define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//'
// diag(W0* Cki * W0^t)
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}

void update_muf(const field<mat>& Xf,const field<vec>& Af, const mat& Z, const mat& M,
                const field<vec>& Xif, const field<mat>& Bf, const field<vec>& betaf, field<vec>& muf){
  int id, im, pm,  d = Xf.n_cols, m=Xf.n_rows;
  mat mat_tmp;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        pm = Xf(im,id).n_cols;
        mat_tmp = (Xf(im,id) - repmat(Af(im,id)+Xif(im,id)+Z*betaf(im,id),1, pm) - M*Bf(im,id).t());
        mat vec_tmp = mean(mat_tmp,0);
        muf(im,id) = conv_to< colvec >::from(vec_tmp);
      }
    }
  }
}

void update_betaf(const field<mat>& Xf,const field<vec>& Af, const mat& Z, const mat& M,
                const field<vec>& Xif, const field<vec>& muf, const field<mat>& Bf,
                const field<vec>& invLambdaf, field<vec>& betaf){
  
  int id, im, pm,  d = Xf.n_cols, m=Xf.n_rows, n=Xf(0).n_rows;
  mat mat_tmp;
  vec vec_tmp;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        pm = Xf(im,id).n_cols;
        mat_tmp = (Xf(im,id) - repmat(Af(im,id)+Xif(im,id),1, pm) - repmat(muf(im,id).t(), n, 1) - M*Bf(im,id).t());
        vec_tmp = (Z.t()*Z).i() * Z.t() * sum(mat_tmp % repmat(invLambdaf(im, id).t(),n,1), 1) / accu(invLambdaf(im, id));
        betaf(im,id) = vec_tmp;
      }
    }
  }
}




void update_Bf(const field<mat>& Xf,const field<vec>& Af, const mat& Z, const mat& M, const field<vec>& Xif,
               const field<vec>& betaf, const field<vec>& muf, const mat& S,  field<mat>& Bf){
  int id, im, pm,  d = Xf.n_cols, m=Xf.n_rows, n=Xf(0).n_rows;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        pm = Xf(im,id).n_cols;
        mat mat_tmp = (Xf(im,id) - repmat(Af(im,id)+Xif(im,id)+Z*betaf(im,id),1, pm) - repmat(muf(im,id).t(), n, 1));
        Bf(im,id) = mat_tmp.t()* M * inv(n*S+M.t() * M);
      }
    }
  }
  
}
void update_invLambdaf(const field<mat>& Xf,const field<vec>& Af, const mat& Z, const mat& M, const field<vec>& Xif,
                       const field<vec>& muf, const field<vec>& betaf,const field<mat>& Bf, const mat& S, const  mat& Om,
                       const field<mat>& Sf_y, field<vec>& invLambdaf){
  
  int id, im, pm,  d = Xf.n_cols, m=Xf.n_rows, n=Xf(0).n_rows;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        pm = Xf(im,id).n_cols;
        mat mat_tmp = (Xf(im,id) - repmat(Af(im,id)+Xif(im,id)+Z*betaf(im,id),1, pm) - repmat(muf(im,id).t(), n, 1)- M*Bf(im,id).t());
        vec vec_tmp = decomp(S, Bf(im,id)) + trans(mean(Sf_y(im,id))) + Om(im,id);
        vec_tmp = vec_tmp + trans(mean(mat_tmp % mat_tmp));
        invLambdaf(im,id) = 1.0 / vec_tmp;
      }
    }
  }
}




void update_Sigmam(const field<vec>& Xif, const  mat& Om, const field<mat>&Bf,  mat& Sigmam){
  int im, id, m = Sigmam.n_rows, d = Sigmam.n_cols;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        Sigmam(im, id) = mean(Xif(im, id)% Xif(im, id)) + Om(im, id);
        
      }
    }
  }
}

void add_IC_Orth(field<mat>& Bf){
  // Try the orthogonal matrix method
  int id, im, qs1,  d = Bf.n_cols, m=Bf.n_rows;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        
        qs1 = Bf(im,id).n_cols;
        mat U1, V1;
        vec s1;
        svd(U1, s1, V1, Bf(im,id));
        vec signU1 = sign(U1.row(0).t());
        Bf(im, id) = U1.cols(0,qs1-1) * diagmat(s1 % signU1.subvec(0,qs1-1));
      }
    }
  }
}

void add_IC_beta(const mat&Z, const field<mat>& Bf, field<vec>& betaf, field<vec>& Xif){
  
  int id, im, qs1,  d = Bf.n_cols, m=Bf.n_rows;
  vec vec_tmp;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        vec_tmp = (Z.t()*Z).i() * Z.t() * Xif(im, id);
        betaf(im,id) += vec_tmp;
        Xif(im, id) -= Z*vec_tmp;
      }
    }
  }
  
}
double ELBOfun(const field<mat>& Xf, const field<vec>& Af, const vec& typeID, const mat& Z,
               const field<mat>& Muf_y,  const field<mat>& Sf_y, const mat& M, const mat& S, 
               const field<vec>& Xif,  const mat& Om,
               const field<mat>& Bf, const field<vec>& betaf, const mat& Sigmam, const field<vec>& muf,
               const field<vec>& invLambdaf, const bool& update_sigma){
  // Basic information:
  int id, im, pm,  d = Xf.n_cols, m=Xf.n_rows, n=Xf(0).n_rows;
  double tmp_v;
  
  // log P(X|Y) 
  double val_tmp0 = 0.0;
  for(id=0; id <d; ++id){ // Loop for all modality group
    if(typeID(id)==2){
      for(im=0; im<m; ++im){ // loop for all modalities in Count-type group.
        if(Bf(im,id).n_cols>1){ 
          val_tmp0 += accu(Xf(im,id) % Muf_y(im,id)-exp(Muf_y(im,id) + Sf_y(im,id)/2));
        }
      }
    }
    if(typeID(id)==3){
      for(im=0; im<m; ++im){ // loop for all modalities in Binomial-type group.
        if(Bf(im,id).n_cols>1){ 
          rowvec nvec = max(Xf(im,id));
          val_tmp0 += accu((Xf(im,id)-repmat(nvec, n, 1)) % Muf_y(im,id));  //
          val_tmp0 += -accu(repmat(nvec, n, 1) %  exp(-Muf_y(im,id) + Sf_y(im,id)/2)); // term2 is not written now!
        }
      }
    }
  }
  
  
  // log P(Y|M,Xi, Z)
  //Rprintf("Good EBLO1!\n");
  double val_tmp1 = 0.0; //
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){ 
        pm = Xf(im,id).n_cols;
        mat mat_tmp = (Muf_y(im,id) - repmat(Af(im,id)+Xif(im,id)+Z*betaf(im,id),1, pm) - repmat(muf(im,id).t(), n, 1)- M*Bf(im,id).t());
        vec vec_tmp = n*decomp(S, Bf(im,id)) + trans(sum(Sf_y(im,id))) + n* Om(im,id);
        val_tmp1 += accu(vec_tmp % invLambdaf(im,id)) + accu(sum(mat_tmp % mat_tmp) % invLambdaf(im,id).t()) - n*accu(log(invLambdaf(im,id)));
      }
    }
  }
  val_tmp1 = -0.5* val_tmp1;
  //Rprintf("Good EBLO2!\n");
  // log P(Xi)
  double val_tmp2 = 0.0;
  if(update_sigma){
    for(id=0; id <d; ++id){ // Loop for all modality group
      for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id).n_cols>1){ 
          tmp_v = log(Sigmam(im,id)+1e-10);
          val_tmp2 += n * tmp_v + accu(Xif(im,id) % Xif(im,id) + Om(im,id))/(Sigmam(im,id));
        }
      }
    }
  }
  val_tmp2 = -0.5* val_tmp2;
  double val_tmp3 = 0.0;
  val_tmp3 =  accu(M%M) + n*trace(S);
  val_tmp3 = -0.5* val_tmp3;
  
  
  // Entropy
  tmp_v = log(det(S));
  double entropy= n* tmp_v;
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){ 
        tmp_v = log(Om(im,id));
        entropy += n*tmp_v + accu(log(Sf_y(im,id)+1e-10));
        
      }
    }
  }
  entropy = 0.5*entropy;
  //Rprintf("Good EBLO4!\n");
  
  return val_tmp0+val_tmp1+val_tmp2+val_tmp3+ entropy;
}



void VB_Estep(const field<mat>& Xf, const field<vec>& Af,const vec& typeID, const mat& Z,
              field<mat>& Muf_y,  field<mat>& Sf_y, mat& M, mat& S, 
              field<vec>& Xif,  mat& Om,
              const field<mat>& Bf, const field<vec> betaf, const mat& Sigmam, const field<vec>& muf,
              const field<vec>& invLambdaf){
  
  // Basic information:
  int id, im, pm,  d = Xf.n_cols, m=Xf.n_rows, n=Xf(0,0).n_rows, q= Bf(0,0).n_cols;
  // vec p_vec(d+1, fill::zeros);
  // Mu_y is a n*p matrix
  //mat Zm = (M+Xi(m)) * Bf(m).t() + repmat(muf(i,m).t(), n, 1) + repmat(Af(m),1, pm);
  //mat invLamMat = repmat(invLambda.t(), n, 1);
  // Rprintf("Good Estep1!\n");
  // double elbo0 =  ELBOfun( Xf, Af, typeID, Muf_y, Sf_y,  M, S, Muf, Xif, Of, Bf, Sigmamf, muf, invLambdaf, tau, alpha);
  // Update Muf_y and Sf_y
  for(id=0; id <d; ++id){ // Loop for all modality group
    // Initialize in the external if typeID(id)==1
    // Rprintf("Good Estep1.1!\n");
    if(typeID(id)==2){
      for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id).n_cols>1){
          pm = Xf(im,id).n_cols;
          mat Mu_y1 = Muf_y(im,id);
          mat tZ = repmat(Af(im,id) + Xif(im,id)+Z*betaf(im, id),1, pm) + repmat(muf(im,id).t(), n, 1)+ M*Bf(im,id).t();
          Muf_y(im,id) = (Xf(im,id) -  exp(Mu_y1) % (1-Mu_y1) + repmat(invLambdaf(im,id).t(), n, 1)% tZ) /
            (exp(Mu_y1) +  repmat(invLambdaf(im,id).t(), n, 1) );
          Sf_y(im,id)= 1.0 / (exp(Muf_y(im,id)) + repmat(invLambdaf(im,id).t(), n, 1) );
        }
      }
    }
    if(typeID(id)==3){
      for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
        if(Bf(im,id).n_cols>1){
          pm = Xf(im,id).n_cols;
          rowvec nvec = max(Xf(im,id));
          mat tZ = repmat(Af(im,id)+ Xif(im,id)+Z*betaf(im, id),1, pm) + repmat(muf(im,id).t(), n, 1)+ M*Bf(im,id).t();
          mat Mu_y2 = Muf_y(im,id); // take out this submatrix
          Mu_y2 = (Xf(im,id)- (1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + repmat(invLambdaf(im,id).t(), n, 1)% tZ) /
            ((1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + repmat(invLambdaf(im,id).t(), n, 1) );
          Muf_y(im,id) = Mu_y2;
          Sf_y(im,id) = 1.0 / ( (1/(1+exp(-Mu_y2))) %(1- 1/(1+exp(-Mu_y2))) % repmat(nvec,n, 1) + repmat(invLambdaf(im,id).t(), n, 1));
        }
      }
    }
  }
  // double elbo1 =  ELBOfun( Xf, Af, typeID, Muf_y, Sf_y,  M, S, Muf, Xif, Of, Bf, Sigmamf, muf, invLambdaf, tau, alpha);
  // Rprintf("dMu_y= %4f\n", elbo1 - elbo0);
  // Update Si and mi
  mat Si_inv(q,q, fill::zeros);
  mat M_tmp(n,q, fill::zeros);
  mat mat_tmp2;
  double val_tmp1;
  //mat mat_tmp3;
  //field<mat> Om_inv(m,d);
  for(id=0; id <d; ++id){ // Loop for all modality group
    for(im=0; im<m; ++im){ // loop for all modalities in each modality group.
      if(Bf(im,id).n_cols>1){
        pm = Xf(im,id).n_cols;
        mat BfL = (Bf(im,id) % repmat(invLambdaf(im,id), 1, q));
        val_tmp1 = accu(invLambdaf(im,id));
        Si_inv += Bf(im,id).t() *BfL;
        mat_tmp2 = (Muf_y(im,id) - repmat(Af(im,id)+Z*betaf(im, id),1, pm) - repmat(muf(im,id).t(), n, 1))*BfL;
        // M_tmp += (Xf(im,id) - repmat(Af(im,id),1, pm) - repmat(muf(im,id).t(), n, 1)-Xif(im,id)*Bf(im,id).t())*BfL;
        mat mat_tmp1 = Muf_y(im,id) - repmat(Af(im,id)+Z*betaf(im, id),1, pm) - repmat(muf(im,id).t(), n, 1) - M*Bf(im,id).t();
        //Om_inv(im,id) = mat_tmp1 + diagmat(1.0/Sigmamf(im,id));
        //Of(im,id) = Om_inv(im,id).i();
        Om(im,id) = 1/(val_tmp1 + 1.0/Sigmam(im,id));
        // mat_tmp3 = Om_inv(im,id).i(); // Only update Xif
        //Xif(im,id) = (mat_tmp2 - (M*Bf(im,id).t())*BfL)* mat_tmp3;
        Xif(im,id) =Om(im,id)* sum(mat_tmp1 % repmat(invLambdaf(im,id).t(),n , 1), 1); // sum by row
        M_tmp += mat_tmp2 -  repmat(Xif(im,id),1, pm)* BfL; // Must first update Xif, then update M_tmp
      }
    }
  }
  
  // double elbo2 =  ELBOfun( Xf, Af, typeID, Muf_y, Sf_y,  M, S, Muf, Xif, Of, Bf, Sigmamf, muf, invLambdaf, tau, alpha);
  // Rprintf("dXi= %4f\n", elbo2 - elbo1);
  // Rprintf("Good Estep2!\n");
  Si_inv +=  eye(q,q);
  S = Si_inv.i(); //
  // double elbo3 =  ELBOfun( Xf, Af, typeID, Muf_y, Sf_y,  M, S, Muf, Xif, Of, Bf, Sigmamf, muf, invLambdaf, tau, alpha);
  // Rprintf("dS= %4f\n", elbo3 - elbo2);
  M = M_tmp* S; // There is some problem in updating M!!!!! Why????????
  // double elbo4 =  ELBOfun(Xf, Af, typeID, Muf_y, Sf_y,  M, S, Muf, Xif, Of, Bf, Sigmamf, muf, invLambdaf, tau, alpha);
  // Rprintf("dM= %4f\n", elbo4 - elbo3);
}


// [[Rcpp::export]]
Rcpp::List vb_cmgfmcpp(const Rcpp::List& XList, const arma::vec& typeID,
                      const arma::mat& numvarmat, const Rcpp::List&  Alist,
                      const arma::mat& Z, const Rcpp::List& Mulist_y_int,
                      const Rcpp::List& Slist_y_int, const Rcpp::List& Sigmamlist_int,
                      const Rcpp::List& invLambdalist_int, const Rcpp::List& Blist_int,
                      const Rcpp::List& mulist_int, const Rcpp::List& betalist_int, 
                      const arma::mat& M_int, const arma::mat& S_int, 
                      const arma::vec& Xi_int, const double& O_int, 
                      const double& epsELBO, const int& maxIter, const bool& verbose,
                      const bool& add_IC_inter=true,
                      const bool& update_sigma=true){
  // typeID: 1 means Gaussian; 2 means Poisson; 3 means Binomial;
  // Xf field represents the modality group with same variable type.
  // numvarmat: the number ofvariables in each modality within modality group.
  // betalist_int[[t]] \in R^(d*M_t)
  
  int i, d = XList.length(), m= numvarmat.n_cols; // d is the number of variable types; m is the max number of modalities with same type.
  int d2 = typeID.n_elem;
  if(d != d2){
    stop("The length of XList must be equal to the length of typeID!");
  }
  // Rprintf("Good entry1!\n");
  // Rprintf("d=%d, m=%d", d, m);
  
  //alpha_grid.print();
  field<mat> Xf(m,d),  Muf_y(m,d),  Sf_y(m,d), Bf(m,d);
  field<vec> muf(m,d), invLambdaf(m,d),  Xif(m,d),  Af(m,d), betaf(m,d);
  mat Om(m,d, fill::zeros), Sigmam(m,d);
  for(i=0; i<d; ++i){ // put the modality group matrix of list into a field.
    mat Xtmp = XList[i];
    mat Mutmp = Mulist_y_int[i];
    mat Stmp = Slist_y_int[i];
    mat Btmp = Blist_int[i];
    mat Atmp = Alist[i];
    mat bbtmp = betalist_int[i];
    vec Sigmamtmp = Sigmamlist_int[i];
    vec mutmp = mulist_int[i];
    vec invLamtmp = invLambdalist_int[i];
    vec p_vec(m+1, fill::zeros); // initilize
    for(int im=0; im<m; ++im){
      if(numvarmat(i,im)>0){
        // Rprintf("i = %d, im=%d \n", i, im);
        Xif(im,i) = Xi_int;
        Om(im,i) = O_int;
        Sigmam(im,i) = Sigmamtmp(im);
        betaf(im,i) = bbtmp.col(im);
        p_vec(im+1) = p_vec(im) + numvarmat(i,im);
        Xf(im,i) = Xtmp.cols(p_vec(im), p_vec(im+1)-1);
        Muf_y(im,i) = Mutmp.cols(p_vec(im), p_vec(im+1)-1);
        // Rprintf("Good entry1!\n");
        Sf_y(im,i) = Stmp.cols(p_vec(im), p_vec(im+1)-1);
        Bf(im,i) = Btmp.rows(p_vec(im), p_vec(im+1)-1);
        muf(im,i) = mutmp.subvec(p_vec(im), p_vec(im+1)-1);
        invLambdaf(im,i) = invLamtmp.subvec(p_vec(im), p_vec(im+1)-1);
        // Rprintf("Good entry1.1!\n");
        Af(im,i) = Atmp.col(im);
      }else{
        Xf(im,i) = zeros(1,1); // fill 1-by-1 zero matrix for other empty position.
        Muf_y(im,i) = Xf(im,i);
        Sf_y(im,i) = Xf(im,i);
        Bf(im,i) = Xf(im,i);
        muf(im,i) =zeros(1,1);
        betaf(im,i) = zeros(1,1);
        invLambdaf(im,i) = muf(im,i);
        Sigmam(im,i) =  0;
        Af(im,i) = muf(im,i);
      }
      
    }
  }
  
  if(!update_sigma){
    Sigmam = zeros(m,d);
  }
  // Rprintf("Good entry2!\n");
  
  
  // Initialize
  mat  M(M_int), S(S_int);
  
  vec ELBO_vec(maxIter);
  ELBO_vec(0) = -1e20;
  mat dX, BL;
  int iter;
  
  for(iter = 1; iter < maxIter; ++iter){
    
    // Rprintf("E step starting!\n");
    // VB E-step
    VB_Estep(Xf, Af, typeID, Z, Muf_y, Sf_y, M, S, Xif,  Om, Bf, betaf, Sigmam, muf, invLambdaf);
    // Rprintf("Finish E step!\n");
    //VB M-step
    
    // double elbo1 = ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
    
    // update mu
    // Rprintf("update mu\n");
    update_muf(Muf_y, Af, Z, M, Xif, Bf, betaf, muf);
    // double elbo2 = ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
    // Rprintf("dmu= %4f \n", elbo2 - elbo1);
    
    // update beta
    // Rprintf("update beta\n");
    update_betaf(Muf_y, Af, Z, M, Xif, muf, Bf,invLambdaf, betaf);
    // double elbo21 = ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
    // Rprintf("dbeta= %4f \n", elbo21 - elbo2);
    
    //update B
    // Rprintf("update B\n");//there is some problem in updating lambda
    update_Bf(Muf_y, Af, Z, M, Xif, betaf, muf, S, Bf);
    if(add_IC_inter){
      add_IC_Orth(Bf); // Add identifiability condition
    }
    
    // double elbo3 = ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
    // Rprintf("dB= %4f \n", elbo3 - elbo21);
    
    // update Lambda
    // Rprintf("update Lambda\n"); // 
    update_invLambdaf(Muf_y, Af, Z, M, Xif, muf,betaf, Bf, S, Om, Sf_y, invLambdaf);
    // double elbo4 =   ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
    // Rprintf("dLambda= %4f\n", elbo4 - elbo3);
    
    // update Sigmam
    if(update_sigma){
      // Rprintf("update Sigmam\n");
      update_Sigmam(Xif, Om, Bf,Sigmam);
      // double elbo5 =  ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
      // Rprintf("dSigmam= %4f\n", elbo5 - elbo4);
    }
    
    
  
    ELBO_vec(iter) = ELBOfun(Xf, Af, typeID,Z, Muf_y, Sf_y,  M, S, Xif, Om, Bf, betaf, Sigmam, muf, invLambdaf,update_sigma);
    
    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n",
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }
  if(!add_IC_inter){
    add_IC_Orth(Bf); // Add identifiability condition
    // VB_Estep(Xf, Af, typeID, Z, Muf_y, Sf_y, M, S, Xif,  Om, Bf, betaf, Sigmam, muf, invLambdaf);
  }
  add_IC_beta(Z, Bf, betaf, Xif);
  
  // output return value
  List resList = List::create(
    Rcpp::Named("betaf") = betaf,
    Rcpp::Named("Bf") = Bf,
    Rcpp::Named("M") = M,
    Rcpp::Named("Xif") = Xif,
    Rcpp::Named("S") = S,
    Rcpp::Named("Om") = Om,
    Rcpp::Named("muf") = muf,
    Rcpp::Named("Sigmam") = Sigmam,
    Rcpp::Named("invLambdaf") = invLambdaf,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);
  
}

