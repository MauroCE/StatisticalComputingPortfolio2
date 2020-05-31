// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// C H O L E S K Y   S O L V E R
// [[Rcpp::export]]
vec choleskySolverArmadillo(mat& X, 
                            vec& y) {
  // Cholesky decomposition of X^t X
  mat L = chol(X.t() * X);
  // Solve the system using Cholesky
  return solve(L, solve(L.t(), X.t() * y));
}


// Q R    S O L V E R
// [[Rcpp::export]]
arma::vec qrSolverArmadillo(arma::mat& X, 
                            arma::vec& y){
  // QR decomposition. We need to initialize Q and R first
  arma::mat Q, R;
  qr_econ(Q, R, X); 
  // Solve the system
  return solve(R, Q.t() * y);
}



// M U L T I V A R I A T E     G A U S S I A N      D E N S I T Y 
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
      for(ii = 0; ii < irow; ii++){
        acc += z.at(ii) * L.at(irow, ii);
      }
      z.at(irow) = ( x.at(icol, irow) - x0.at(irow) - acc ) / D.at(irow);
    }
    out.at(icol) = sum(square(z));
  }

  out = exp( - 0.5 * out - ( (d / 2.0) * log(2.0 * M_PI) + sum(log(D)) ) );
  
  return out;
}




//[[Rcpp::export]]
vec local_Rcpp(arma::vec& y, 
           arma::mat& x0, 
           arma::mat& X0, 
           arma::mat& x, 
           arma::mat& X, 
           arma::mat& H){
  // Number of points we are doing local polynomial regression on
  unsigned int nsub = x0.n_rows;
  // Cholesky decomposition of the matrix L
  mat L = chol(H, "lower");
  // Store the final output here
  vec out(nsub, fill::zeros); 
  // Store matrices for the QR decomposition
  arma::mat Q, R;
  // Main Loop
  for (int i=0; i < nsub; i++){
    // Square root of the weights
    vec sw = sqrt(dmvnInt(x, x0.row(i), L));
    // QR DECOMPOSITION for weighted least squares
    qr_econ(Q, R, X.each_col() % sw); 
    out(i) = dot(
      X0.row(i), 
      solve(R, Q.t() * (y.each_col() % sw))
    );
  }
  return out;
}


//[[Rcpp::export(name="local_Rcpp2")]]
vec local_Rcpp2(arma::vec& y, 
               arma::mat& x0, 
               arma::mat& X0, 
               arma::mat& x, 
               arma::mat& X, 
               arma::mat& H){
  // Number of points we are doing local polynomial regression on
  unsigned int nsub = x0.n_rows;
  // Cholesky decomposition of the matrix L
  mat L = chol(H, "lower");
  // Store the final output here
  vec out(nsub, fill::zeros); 
  // Store matrices for the QR decomposition
  arma::mat Q, R;
  // Main Loop
  for (int i=0; i < nsub; i++){
    vec w = dmvnInt(x, x0.row(i), L);
    out(i) = dot(X0.row(i), solve((X.each_col() % w).t() * X, X.t() * (w % y)));
  }
  return out;
}



