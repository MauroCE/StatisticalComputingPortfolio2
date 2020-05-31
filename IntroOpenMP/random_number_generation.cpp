#include <iostream>
#include <random>
#include <cmath>
#include "random_number_generation.h"


// Generates a random integer between two integers
int random_int(const int from_integer, const int to_integer){
  std::random_device                  my_device;
  std::mt19937                        my_generator(my_device());
  std::uniform_int_distribution<int>  my_distribution(from_integer, to_integer);
  return my_distribution(my_generator);
}

// Generates a random number between 0 and 1
double random_uniform01(){
  std::random_device                     my_device;
  std::mt19937                           my_generator(my_device());
  std::uniform_real_distribution<double> my_distribution(0, 1);
  return my_distribution(my_generator);
}

// Generates a sample from a Bernoulli distribution with probability p
int random_bernoulli(double p){
  std::random_device                     my_device;
  std::mt19937                           my_generator(my_device());
  std::bernoulli_distribution            my_distribution(p);
  return my_distribution(my_generator);
}


// Generates univariate normal samples N(mean, scale^2)
double random_normal(double mean, double scale){
  std::random_device                     my_device;
  std::mt19937                           my_generator(my_device());
  std::normal_distribution<float>        my_distribution(mean, scale);
  return my_distribution(my_generator);
}

double random_normal_boxmuller(double mean, double scale){
  // Generate two uniform random numbers between 0 and 1
  double uniform1 = random_uniform01();
  double uniform2 = random_uniform01();
  // Transform them using Box, Muller transformation
  double normal = sqrt(-2.0 * log(uniform1)) * cos(2*M_PI*uniform2);
  return mean + scale*normal;
}


std::vector<int> sample_with_replacement(std::vector<double> probabilities, int n){
  // Initiate the sampler
  std::random_device                     my_device;
  std::mt19937                           my_generator(my_device());
  std::discrete_distribution<>           my_distribution(probabilities.begin(), probabilities.end());
  // Create a vector that will store the samples
  std::vector<int> samples(n);
  // Keep sampling and storing the samples in the vector "samples".
  for (int i=0; i < n; i++){
    samples[i] = my_distribution(my_generator);
  }
  return samples;
}

// int main() {
//   std::cout << random_int(10, 20) << std::endl;
//   std::cout << random_uniform01() << std::endl;
//   std::cout << random_bernoulli(0.5) << std::endl;
//   std::cout << random_normal(100.0, 1.0) << std::endl;
//   std::cout << random_normal_boxmuller(-10, 0.1) << std::endl;
// }