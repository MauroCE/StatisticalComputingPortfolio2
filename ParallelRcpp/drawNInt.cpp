#include <Rcpp.h>
#include <sitmo.h> 
using namespace Rcpp;

// [[Rcpp::depends(sitmo)]]

// [[Rcpp::export]]
NumericVector drawNInt(unsigned int n, unsigned int seed) {
  // This function simply draws an integer between 0 and sitmo::prng::max()
  NumericVector draws(n);
  
  // Firstly we set up the engine. When this is called, it will generate a random number
  // between 0 and sitmo::prng::max(). We need to feed into the engine a seed, but we need to cast it
  uint32_t casted_seed = static_cast<uint32_t>(seed);
  sitmo::prng engine(casted_seed);
  
  // Draw n times and write onto the "draws" vector
  for (unsigned int i=0; i < n; i++){
    draws(i) = engine();
  }
  return draws;
}