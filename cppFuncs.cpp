#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]




LogicalVector isNA(NumericVector x) {
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  return out;
}


// [[Rcpp::export]]
uvec isNotNArowvec(const rowvec x) {
  
  int n = x.size();
  vec out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  uvec ids = find(out==0);
  //  return x.elem(ids);
  return(ids);
}

// [[Rcpp::export]]
uvec isNotNAcolvec(const colvec x) {
  
  int n = x.size();
  vec out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = NumericVector::is_na(x[i]);
  }
  uvec ids = find(out==0);
  //  return x.elem(ids);
  return(ids);
}



// [[Rcpp::export]]
mat fastLmResid(const mat& Y, const mat& X){
  
  int n = X.n_rows, k = X.n_cols;
  mat res;
  
  mat coef = Y*X*inv(trans(X)*X);    // fit model y ~ X
  res  = Y - coef*trans(X);           // residuals
  
  return res;
}


// [[Rcpp::export]]
mat fastLmResidWeighted(const mat& Y, const mat& X,  const rowvec& wa){
  
  int n = X.n_rows, k = X.n_cols;
  mat res;
  
  //rowvec ws=rowvec(wa.n_elem);
  //for (int j=0; j<ws.n_elem; j++){
  //  ws[j]=std::sqrt(wa[j]);
  //}
  mat W=diagmat(wa);
  
  // coeff=dat%*%W%*%modtmp %*% solve(t(modtmp) %*% W %*% modtmp)
  mat coef = Y*W*X*inv(trans(X)*W*X);    // fit model y ~ X
  res  = Y - coef*trans(X);           // residuals
  //res.each_row()%=ws;
  
  return res;
}  


// [[Rcpp::export]]
mat fastLmPredicted(const mat& Y, const mat& X){
  
  int n = X.n_rows, k = X.n_cols;
  mat res;
  
  mat coef = Y*X*inv(trans(X)*X);    // fit model y ~ X
  return coef*trans(X);           // residuals
  
  
}







// [[Rcpp::export]]
List fastLm(const mat& Y, const mat& X) {
  int n = X.n_rows, k = X.n_cols;
  
  //  coeff=data[i,]%*%mod %*% solve(t(mod) %*% mod)
  //    resid[i, ] = data[i,] -(coeff %*% t(mod))\
  
  
  mat coef = Y*X*inv(trans(X)*X);    // fit model y ~ X
  mat res  = Y - coef*trans(X);           // residuals
  
  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
  
  // colvec std_err = sqrt(s2 * diagvec(pinv(trans(X)*X)));
  
  return List::create(Named("coefficients") = coef,
                      //                   Named("stderr")       = std_err,
                      Named("df.residual")  = n - k);
}

// [[Rcpp::export]]
mat fastLmResidMat(const mat& Y, const mat& X) {
  uvec ids;
  mat rmat=mat(Y.n_rows, Y.n_cols);
  rmat.fill(datum::nan);
  uvec vec_i=uvec(1);
  int i,j;
  for (int i=0; i<Y.n_rows; i++){
    vec_i[0]=i;
    ids=isNotNArowvec(Y.row(i));
    if(ids.n_elem>X.n_cols){
      rmat.submat(vec_i, ids)=fastLmResid(Y.submat(vec_i, ids), X.rows(ids));
    }
  }
  return(rmat);
}


// [[Rcpp::export]]
mat fastLmPredictedMat(const mat& Y, const mat& X) {
  uvec ids;
  mat rmat=mat(Y.n_rows, Y.n_cols);
  rmat.fill(datum::nan);
  uvec vec_i=uvec(1);
  int i,j;
  for (int i=0; i<Y.n_rows; i++){
    vec_i[0]=i;
    ids=isNotNArowvec(Y.row(i));
    if(ids.n_elem>X.n_cols){
      rmat.submat(vec_i, ids)=fastLmPredicted(Y.submat(vec_i, ids), X.rows(ids));
    }
  }
  return(rmat);
}




// [[Rcpp::export]]
mat fastLmResidMatWeighted(const mat& Y, const mat& X, const mat& W) {
  uvec ids;
  mat rmat=mat(Y.n_rows, Y.n_cols);
  rmat.fill(datum::nan);
  uvec vec_i=uvec(1);
  int i,j;
  for (int i=0; i<Y.n_rows; i++){
    vec_i[0]=i;
    ids=isNotNArowvec(Y.row(i));
    if(ids.n_elem>X.n_cols){
      rmat.submat(vec_i, ids)=fastLmResidWeighted(Y.submat(vec_i, ids), X.rows(ids), W.submat(vec_i, ids));
    }
  }
  return(rmat);
}



// [[Rcpp::export]]
mat fastLmResidMatWeightedTrans(const mat& Y, const mat& X, const mat& W) {
  uvec ids;
  mat rmat=mat( Y.n_rows,Y.n_cols);   
  rmat.fill(datum::nan);
  uvec vec_i=uvec(1);
  mat Yt=trans(Y);

  int i,j;
  for (int i=0; i<Y.n_rows; i++){
    vec_i[0]=i;  
  //  cerr<<"Here"<<endl;
    ids=isNotNAcolvec(Yt.col(i));
//cerr<<"Here2"<<endl;
//cerr<<i<<endl;
    if(ids.n_elem>X.n_cols){
      rmat.submat( vec_i, ids)=(fastLmResidWeighted(Y.submat(  vec_i, ids), X.rows(ids), W.submat(vec_i, ids)));
    }
  }
  return(rmat);
}



// [[Rcpp::export]]
mat fastLmResidWeightedPredict(const mat& Y, const mat& X,  const rowvec& wa, const mat& newX){
  
  int n = X.n_rows, k = X.n_cols;
  mat res;
  

  mat W=diagmat(wa);
  

  mat coef = Y*W*X*inv(trans(X)*W*X);    // fit model y ~ X
  res  = coef*trans(newX);           // predictions
  
  
  return res;
}  


// [[Rcpp::export]]
mat fastLmResidMatWeightedPredict(const mat& Y, const mat& X, const mat& W, const mat& newX) {
  uvec ids;
  mat rmat=mat(Y.n_rows, newX.n_rows);
  
  uvec vec_i=uvec(1);
  int i,j;
  int nc=rmat.n_cols-1;
  int nrNew=newX.n_rows-1;
  
  for (int i=0; i<Y.n_rows; i++){
    vec_i[0]=i;
    ids=isNotNArowvec(Y.row(i));
    //    rmat.submat(vec_i, ids)=fastLmResidWeighted(Y.submat(vec_i, ids), X.rows(ids), W.submat(vec_i, ids));
    rmat.submat(span(i,i),span(0,nrNew))=fastLmResidWeightedPredict(Y.submat(vec_i,ids), X.rows(ids), W.submat(vec_i,ids), newX);
  }
  return(rmat);
}

/*** R
set.seed(123)
x=rnorm(10)
mod=model.matrix(~1+x)
y=rbind(rnorm(10))
w=rbind(rnorm(10)^2)
y[1,1]=NA
#make sure this is working
message("Running checks")
fastLmResidMatWeightedPredict(rbind(y), mod, rbind(w), rbind(mod, mod))

testpred=fastLmResidMatWeightedPredict(rbind(y), mod, rbind(w), rbind(mod))



#test normal lm
iinotna=which(!is.na(y))

#do normal lm no weights
lmres=lm(y[iinotna]~0+mod[iinotna,])
rvals=stats::resid(lmres)
myvals=naresidCPP(rbind(y[iinotna]), mod[iinotna,])
stopifnot(mean(abs(rvals-myvals))<1e-6)



#do normal lm with weights
lmres=lm(y[iinotna]~0+mod[iinotna,], weights = w[iinotna])
rvals=stats::resid(lmres)
myvals=naresidCPP(rbind(y[iinotna]), mod[iinotna,],weights = rbind(w[iinotna]))
#naresidCPP includes a stadartization step
stopifnot(mean(abs(rvals*sqrt(w[iinotna])-myvals))<1e-6)


#generate full predictions
lmpred= mod%*%coef(lmres)
#make sure they are the same
stopifnot(sum(abs(rbind(t(lmpred))-testpred))<1e-6)


myvalscpp=naresidCPP(rbind(y[iinotna]), mod[iinotna,],weights = rbind(w[iinotna]))
#naresidCPP also  includes a stadartization step
stopifnot(mean(abs(rvals*sqrt(w[iinotna])-myvals))<1e-6)


  */

