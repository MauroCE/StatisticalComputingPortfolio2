#include <iostream>
#include <random>
#include <cmath>
#include "random_number_generation.h"


// double normal_density(double x, double mean, double scale){
//   return (1/sqrt(2*M_PI)*scale) * exp(-std::pow((x - mean)/scale, 2.0) / 2);
// }

double log_normal_density(double x, double mean, double scale){
  return -0.5*log(2*M_PI*pow(scale, 2)) - pow((x - mean) / scale, 2) / 2;
}

std::vector<double> mh(int n_samples, double start, double target_mean, double target_scale,
                       double proposal_scale){
  // Generate a vector where the samples will be stored.
  // Evaluate target at the starting value. Store starting value
  std::vector<double> samples(n_samples);
  double p_current_sample = log_normal_density(start, target_mean, target_scale);
  double current_sample = start;
  
  // Generate uniform random numbers and take the log of them. Will be used for acceptance.
  std::vector<double> log_uniform_samples(n_samples);
  for (int i=0; i < n_samples; i++){
    log_uniform_samples[i] = log(random_uniform01());
  }
  
  // Generate random noise, will be added to the current sample to generate candidates
  std::vector<double> normal_noise(n_samples);
  for (int i=0; i < n_samples; i++){
    normal_noise[i] = random_normal(0.0, proposal_scale);
  }
  
  // Main Loop
  for (int i=0; i < n_samples; i++){
    // Generate candidate and evaluate the target 
    double candidate = current_sample + normal_noise[i];
    double p_candidate = log_normal_density(candidate, target_mean, target_scale);
    // Check whether to accept it or not
    if (log_uniform_samples[i] <= p_candidate - p_current_sample){
      // ACCEPT!
      current_sample = candidate;
      p_current_sample = p_candidate;
    }
    // Store the sample
    samples[i] = current_sample;
  }
  return samples;
}

int main(int n_command_line_arguments_including_filename, char** command_line_arguments){
  // Use stoi to transform command line arguments to integers.
  // Use stod to transform command line arguments to double.
  int n_iterations = std::stoi(command_line_arguments[1]);
  double start = std::stod(command_line_arguments[2]);
  double target_mean = std::stod(command_line_arguments[3]);
  double target_scale = std::stod(command_line_arguments[4]);
  double proposal_scale = std::stod(command_line_arguments[5]);
  
  // generate vector
  std::vector<double> my_vector = mh(n_iterations, start, target_mean, target_scale, proposal_scale);
  int vector_size = my_vector.size();
  // print elements of the vector
  for (int i=0; i < vector_size; i++){
    std::cout << my_vector[i] << std::endl;
  }
  return 0;
}