#include <iostream>
#include <random>
#include <cmath>
#include "random_number_generation.h"
#include <limits>


#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

// Define a function that sums two integers and then computes the modulo 7 of the result
// This function will be used to perform reduction. The first argument will be the already 
// reduced value, the second argument will be the new value.
int modulo7_sum(int already_reduced, int newvalue){
  return (already_reduced + newvalue) % 7;
}

// Declare/Define a new reduction
#pragma omp declare reduction (modulo7_summation : int : omp_out = modulo7_sum(omp_out, omp_in)) initializer(omp_priv=0)

int main(int argc, char **argv){
  int n = std::stoi(argv[1]);
  int tot = 0;
  
  #pragma omp parallel reduction(modulo7_summation: tot)
  {
    int one = 0;
    #pragma omp for
    for (int i=0; i < n; i++){
      ++one;
    }
    tot = modulo7_sum(tot, one);
  }
  std::cout << tot << std::endl;
}