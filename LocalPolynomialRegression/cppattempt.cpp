// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// CHOLESKY DECOMPOSITION
// [[Rcpp::export]]
arma::vec choleskySolverArmadillo(arma::mat& X, 
                                  arma::vec& y) {
  // Cholesky decomposition of X^t X
  arma::mat L = arma::chol(X.t() * X);
  // Solve the system using Cholesky
  return solve(L, solve(L.t(), X.t() * y));
}

// QR DECOMPOSITION
// [[Rcpp::export]]
arma::vec qrSolverArmadillo(arma::mat& X, 
                            arma::vec& y){
  // QR decomposition. We need to initialize Q and R first
  arma::mat Q, R;
  qr_econ(Q, R, X); 
  // Solve the system
  return solve(R, Q.t() * y);
}

// MULTIVARIATE GAUSSIAN DENSITY
// [[Rcpp::export(name="rcppgaussian")]]
arma::vec dmvnInt(arma::mat& x, 
                  const arma::rowvec& x0, 
                  arma::mat& L)
{
  unsigned int d = x.n_cols;
  unsigned int m = x.n_rows;
  
  vec D = L.diag();
  vec out(m);
  vec z(d);
  
  double acc;
  unsigned int icol, irow, ii;
  for(icol = 0; icol < m; icol++)
  {
    for(irow = 0; irow < d; irow++)
    {
      acc = 0.0;
      for(ii = 0; ii < irow; ii++) acc += z.at(ii) * L.at(irow, ii);
      z.at(irow) = ( x.at(icol, irow) - x0.at(irow) - acc ) / D.at(irow);
    }
    out.at(icol) = sum(square(z));
  }
  
  out = exp( - 0.5 * out - ( (d / 2.0) * log(2.0 * M_PI) + sum(log(D)) ) );
  
  return out;
}


// LOCAL LEAST SQUARES FIT
// [[Rcpp::export(name="lsfit_rcpp")]]
int LocalLSFit(arma::vec& y, 
               arma::mat& x0, 
               arma::mat& X0, 
               arma::mat& x, 
               arma::mat& X, 
               arma::mat& H){
  // grab some values
  int nsub = x0.n_rows;
  // MULTIVARIATE NORMAL DENSITY EVALUATION
  arma::mat L = arma::chol(H, "lower");
  // CREATE THE MAIN LOOP
  for (int i=0; i < nsub; i++){
    // Grab x0[ii, ], X0[ii, ]
    arma::vec w = dmvnInt(x, x0.row(i), L);
    arma::mat W = diagmat(w);
  }
  return nsub;
}

//[[Rcpp::export]]
Rcpp::NumericVector something(arma::vec& y, 
              arma::mat& x0, 
              arma::mat& X0, 
              arma::mat& x, 
              arma::mat& X, 
              arma::mat& H){
  // Number of points we are doing local polynomial regression on
  int nsub = x0.n_rows;
  // Cholesky decomposition of the matrix L
  mat L = chol(H, "lower");
  // Store the final output here
  vec final_out(3, fill::zeros);  // should be nsub
  // Loop
  for (int i=0; i < 3; i++){
    mat W = diagmat(dmvnInt(x, x0.row(i), L));
    mat X_weighted = X.t()*W*X;
    vec y_weighted = X.t()*(W*y);
    vec beta = qrSolverArmadillo(X_weighted, y_weighted);
    rowvec X0_row = X0.row(i);
    vec result = X0_row * beta;
    final_out(i) = result(0);
  }
  return Rcpp::wrap(final_out);
}


//[[Rcpp::export]]
vec noloop(arma::vec& y, 
           arma::mat& x0, 
           arma::mat& X0, 
           arma::mat& x, 
           arma::mat& X, 
           arma::mat& H){
  // Number of points we are doing local polynomial regression on
  int nsub = x0.n_rows;
  // Cholesky decomposition of the matrix L
  mat L = chol(H, "lower");
  // Store the final output here
  vec out(2, fill::zeros);  // should be nsub
  // Loop
  for (int i=0; i < 2; i++){
    mat W = diagmat(dmvnInt(x, x0.row(i), L));
    mat X_w = X.t()*W*X;
    vec y_w = X.t()*(W*y);
    out(i) = dot(X0.row(i), 
                 qrSolverArmadillo(X_w, y_w)
              );
  }
  return out;
}


//[[Rcpp::export]]
vec noloop(arma::vec& y, 
           arma::mat& x0, 
           arma::mat& X0, 
           arma::mat& x, 
           arma::mat& X, 
           arma::mat& H){
  // Number of points we are doing local polynomial regression on
  int nsub = x0.n_rows;
  // Cholesky decomposition of the matrix L
  mat L = chol(H, "lower");
  // Store the final output here
  vec out(2, fill::zeros);  // should be nsub
  // Loop
  for (int i=0; i < 2; i++){
    mat W = diagmat(dmvnInt(x, x0.row(i), L));
    mat X_w = X.t()*W*X;
    vec y_w = X.t()*(W*y);
    out(i) = dot(X0.row(i), 
        qrSolverArmadillo(X_w, y_w)
    );
  }
  return out;
}



//[[Rcpp::export]]
// rowvec newnoloop(arma::vec& y, 
//            arma::mat& x0, 
//            arma::mat& X0, 
//            arma::mat& x, 
//            arma::mat& X, 
//            arma::mat& H){
//   // Number of points we are doing local polynomial regression on
//   int nsub = x0.n_rows;
//   // Cholesky decomposition of the matrix L
//   mat L = chol(H, "lower");
//   // store the qr decomposition at every iteration
//   mat Q, R;
//   // Store the final output here
//   vec out(2, fill::zeros);  // should be nsub
//   // Loop
//   int i=0;
//   rowvec rw = sqrt(dmvnInt(x, x0.row(i), L));
//   // mat X1 = rw % X;
//   // qr_econ(Q, R, X1);
//   // mat Q1 = (1 / rw) % Q;
//   // mat Q2 = rw % Q;
//   // mat H = Q1 * Q2.t();
//   // out(i) = H * y;
//   return rw;
// }



// y is m, d matrix where we want to evaluate the kde, 
// basically our x_0
// x is n, d matrix of original samples basically our xi
// H bandwidth matrix
// [[Rcpp::export]]
Rcpp::NumericVector LSfit(arma::vec& y,
                          arma::mat& x0, 
                          arma::mat& xi, 
                          arma::mat& H) {
  unsigned int n = xi.n_rows;
  unsigned int m = x0.n_rows;
  arma::vec out(m, fill::zeros);
  arma::mat cholDec = chol(H, "lower");
  for(int ii = 0; ii < n; ii++){
    out += dmvnInt(x0, xi.row(ii), cholDec);
  }
  out /= n;
  return Rcpp::wrap(out);
}


// y is m, d matrix where we want to evaluate the kde
// x is n, d matrix of original samples
// H bandwidth matrix
// [[Rcpp::export(name = "kdeArma")]]
Rcpp::NumericVector kde_i(arma::mat& y, 
                          arma::mat& x, 
                          arma::mat& H) {
  unsigned int n = x.n_rows;
  unsigned int m = y.n_rows;
  vec out(m, fill::zeros);
  mat cholDec = chol(H, "lower");
  for(int ii = 0; ii < n; ii++){
    out += dmvnInt(y, x.row(ii), cholDec);
  }
  out /= n;
  return Rcpp::wrap(out);
}



















