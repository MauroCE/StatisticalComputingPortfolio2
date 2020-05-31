#include <iostream>
#include <random>
#include <cmath>
#include "random_number_generation.h"

  
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif



int main(int argc, char **argv){
  // Total number of simulations = n_square
  int n_square = std::stoi(argv[1]);
  int n_circle = 0;

#pragma omp parallel reduction( + : n_circle )
  {
    int private_n_circles = 0;
    
    #pragma omp for
    for (int i=0; i<n_square; i++){
      // Compute x and y coordinates
      double x_coord = 2*random_uniform01() - 1;
      double y_coord = 2*random_uniform01() - 1;
      // Check if it is inside the circle
      if (x_coord*x_coord + y_coord*y_coord <= 1){
        private_n_circles += 1;
      }
    }
    
    #pragma omp critical
    {
      n_circle += private_n_circles;
    }
  }
  double pi = 4*(double)n_circle / (double)n_square;
  std::cout << pi << std::endl;
}