#include <iostream>
#include <random>
#include <cmath>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

int main(int argc, char **argv){
  // Run a simple loop and print stuff
  #pragma omp parallel for
  for (int i=0; i < 10; i++){
    std::cout << "Thread: " << omp_get_thread_num() << "Loop Iteration: " << i << std::endl;
  }
}
   