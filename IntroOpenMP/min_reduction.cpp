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


int main(int argc, char **argv){
  // Read user inputs. These define how many random integers to generate and between which integers.
  int from = std::stoi(argv[1]);
  int to = std::stoi(argv[2]);
  int n_integers = std::stoi(argv[3]);
  
  // Generate a vector with such characteristics
  std::vector<int> vector1(n_integers);
  std::vector<int> vector2(n_integers);
  #pragma omp parallel
  {
    #pragma omp for
    for (int i=0; i<n_integers; i++){
      vector1[i] = random_int(from, to);
      vector2[i] = random_int(from, to);
    }
  }
  
  int sum_sum = 0;
  int diff_diff = 0;
  int prod_prod = 0;
  int min_max = INT32_MIN;
  
  // Sum vectors and then sum all the elements of this new vector.
  #pragma omp parallel reduction(+ : sum_sum)
  {
    int elementwise_sum = 0;
    
    #pragma omp for
    for (int i=0; i < n_integers; i++){
      elementwise_sum = vector1[i] + vector2[i];  // sum
    }
    sum_sum = sum_sum + elementwise_sum;
  }

  // Subtract vectors and then subtract all the elements of this new vector
  #pragma omp parallel reduction(- : diff_diff)
  {
    int elementwise_diff = 0;
    
    #pragma omp for
    for (int i=0; i < n_integers; i++){
      elementwise_diff = vector1[i] - vector2[i];  // diff
    }
    diff_diff = diff_diff - elementwise_diff;
  }
  
  // Subtract vectors and then subtract all the elements of this new vector
  #pragma omp parallel reduction(* : prod_prod)
  {
    int elementwise_prod = 0;
    
    #pragma omp for
    for (int i=0; i < n_integers; i++){
      elementwise_prod = vector1[i] * vector2[i];  // prod
    }
    prod_prod = prod_prod * elementwise_prod;
  }
  
  // Subtract vectors and then subtract all the elements of this new vector
  #pragma omp parallel reduction(max : min_max)
  {
    int elementwise_min = 0;
    
    #pragma omp for
    for (int i=0; i < n_integers; i++){
      elementwise_min = std::min(vector1[i], vector2[i]);  // min
    }
    min_max = std::max(min_max, elementwise_min);
  }
  
  // Print the results
  std::cout << sum_sum << std::endl;
  std::cout << diff_diff << std::endl;
  std::cout << prod_prod << std::endl;
  std::cout << min_max << std::endl;
}

