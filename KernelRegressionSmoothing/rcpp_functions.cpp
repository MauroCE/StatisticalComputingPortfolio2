#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// Kernel Regression Smoothing
// [[Rcpp::export(name="meanKRS_Rcpp")]]
NumericVector krs(const NumericVector x, const NumericVector y,
                  const NumericVector x0, const double lam){
  // Grab the size of the vectors x and x0
  int n  = x.size();
  int n0 = x0.size();
  
  // Create vector where estimates of conditional expectation will be stored
  Rcpp::NumericVector out(n0);
  
  for (int i=0; i < n0; i++){
    // Kernel value
    NumericVector k = dnorm(x, x0[i], lam);
    // Compute KRS and store it in out
    out[i] = sum(k * y) / sum(k);
  }
  return out;
}

// Mean Squared Error
double mean_squared_error(const NumericVector a, const NumericVector b){
  return mean(pow(a - b, 2));
}


// Adaptive Kernel Regression Smoothing
// [[Rcpp::export(name="mean_var_KRS_Rcpp")]]
NumericVector krs_var(const NumericVector x, const NumericVector y,
                      const NumericVector x0, const double lam){
  // Grab the size of the vectors x and x0
  int n  = x.size();
  int n0 = x0.size();
  
  // Create vector where estimates of conditional expectation will be stored
  NumericVector mu(n);
  NumericVector madHat(n0);
  NumericVector out(n0);
  
  // FIRST KRS
  for (int i=0; i < n; i++){
    // Kernel value
    NumericVector k = dnorm(x, x[i], lam);
    // Compute KRS and store it in out
    mu[i] = sum(k * y) / sum(k);
  }
  // RESIDUALS
  NumericVector resAbs = abs(y - mu);
  for (int i=0; i < n0; i++){
    NumericVector k2 = dnorm(x, x0[i], lam);
    madHat[i] = sum(k2 * resAbs) / sum(k2);
  }
  // NORMALIZED WEIGHTS
  NumericVector w_unnorm = 1 / madHat;
  NumericVector w = w_unnorm / mean(w_unnorm);
  for (int i=0; i < n0; i++){
    NumericVector k3 = dnorm(x, x0[i], lam*w[i]);
    out[i] = sum(k3 * y) / sum(k3);
  }
  return out;
}

// Does some pointer stuff
typedef NumericVector (*funcPtr)(const NumericVector x, const NumericVector y, const NumericVector x0, const double lam);

// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "mean")
    return(XPtr<funcPtr>(new funcPtr(&krs)));
  else if (fstr == "var")
    return(XPtr<funcPtr>(new funcPtr(&krs_var)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}


// [[Rcpp::export]]
NumericVector crossvalidation_Rcpp(const NumericVector x, const NumericVector y, 
                       const int k, const NumericVector lambdas, const std::string funname){
  // CHOOSE WHICH FUNCTION TO USE
  XPtr<funcPtr> xpfun = putFunPtrInXPtr(funname);
  funcPtr fun = *xpfun;
  // CREATE SHUFFLED FOLD INDECES
  int n = x.size();
  IntegerVector fold_indeces_ordered = rep_len(seq(0, k-1), n);
  IntegerVector fold_indeces = sample(fold_indeces_ordered, n);

  // STORE BEST MEAN MSE AND THE INDEX OF LAMBDA ACHIEVING SUCH A MEAN MSE
  NumericVector all_mean_mses(lambdas.size());
  
  // LOOP THROUGH LAMBDAS
  for (int lam_index = 0; lam_index < lambdas.size(); lam_index++){
    double lam = lambdas[lam_index];  // Grab current lambda
    double mean_mse = 0.0;            // Stores the mean MSE for this lambda. The mean is over all folds.
    
    // LOOP THROUGH FOLDS
    for (int fold = 0; fold < k; fold++){
      LogicalVector is_in_fold = fold_indeces == fold; // Flag. TRUE if datum in holdout fold, FALSE if in training fold.
      // Holdout data
      NumericVector x_holdout = x[is_in_fold];
      NumericVector y_holdout = y[is_in_fold];
      // Training data
      NumericVector x_kept = x[!is_in_fold];
      NumericVector y_kept = y[!is_in_fold];
      // Kernel Regression Smoothing using holdout as training. Compute MSE against y_holdout
      NumericVector krs_value = fun(x_kept, y_kept, x_holdout, lam);
      double fold_mse = mean_squared_error(krs_value, y_holdout);
      // Iteratively compute the mean
      mean_mse += fold_mse / k; // The mean is over all folds.
    }
    all_mean_mses[lam_index] = mean_mse;
  }
  return all_mean_mses;
  
}




