#include <iostream>
#include <random>
#include <cmath>
#include "random_number_generation.h"

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

// Prior will be N(0, 10)
std::vector<double> sample_prior(int n_particles, double mean, double scale){
  // Create vector that will store the samples
  std::vector<double> prior_samples(n_particles);
  // In a loop set the elements of prior_samples to samples from the normal distribution
  #pragma omp parallel for
  for (int particle_index=0; particle_index < n_particles; particle_index++){
    prior_samples[particle_index] = random_normal(mean, scale);
  }
  return prior_samples;
}

// Transition distribution to generate x_t from x_{t-1}
double transition(double x, int t, double transition_scale){
  double mean = 0.5*x + 25*x/(1+std::pow(x, 2)) + 8*cos(1.2*t);
  return random_normal(mean, transition_scale);
}

// Emission distribution to generate y_t from x_t
double emission(double x, double emission_scale){
  double mean = pow(x, 2)/20;
  return random_normal(mean, emission_scale);
}

// Evaluates the pdf of a normal density. Used by likelihood
double normal_density(double x, double mean, double scale){
  return (1/(sqrt(2*M_PI)*scale))*exp(-pow((x - mean)/scale, 2)/2);
}

// Likelihood (essentially we compute the same density from which we sample from in emission)
double likelihood(double y, double x, double scale){
  return normal_density(y, pow(x, 2)/20, scale);
}

std::vector<std::vector<double> > generate_data(int t_max, double transition_scale, 
                                             double emission_scale, double prior_scale){
  // Generate empty storage matrix
  std::vector<std::vector<double> > data(t_max+1, std::vector<double>(2));
  // Generate the x_0
  data[0][0] = random_normal(0.0, prior_scale);
  for (int i=1; i < t_max+1; i++){
    // Move forward in time with p(x_t | x_{t-1})
    data[i][0] = transition(data[i-1][0], i, transition_scale);
    // Create observation
    data[i][1] = emission(data[i][0], emission_scale);
  }
  return data;
}


int main(int n_command_line_arguments_including_filename, char** command_line_arguments){
  int t_max = std::stoi(command_line_arguments[1]);
  int n_particles = std::stoi(command_line_arguments[2]);
  // Generate dataset
  std::vector<std::vector<double> > data = generate_data(t_max, sqrt(10), sqrt(1.0), sqrt(10));
  // Initialize particle storage matrix
  std::vector<std::vector<double> > particles(t_max+1, std::vector<double>(n_particles));
  // Sample the initial particles
  std::vector<double> prior_sample = sample_prior(n_particles, 0.0, sqrt(10));
  // Set the first row of "particles" to the "prior_sample"
  #pragma omp parallel for
  for (int i=0; i<n_particles; i++){
    particles[0][i] = prior_sample[i];
  }
  // Do bootstrap filter for every time step
  for (int i=1; i<t_max+1; i++){
    // Sample from the transition and store in the new line of particles
    #pragma omp parallel for
    for (int j=0; j < n_particles; j++){
      particles[i][j] = transition(particles[i-1][j], i, sqrt(10));
    }
    // Compute importance weights
    std::vector<double> weights(n_particles);
    #pragma omp parallel for
    for (int j=0; j < n_particles; j++){
      weights[j] = likelihood(data[i][1], particles[i][j], sqrt(1.0));
    }
    // Normalize importance weights
    double weights_sum = 0.0;
    #pragma omp parallel for reduction(+ : weights_sum)
    for (int j=0; j < n_particles; j++){
      weights_sum += weights[j];
    }
    #pragma omp parallel for
    for (int j=0; j < n_particles; j++){
      weights[j] = weights[j] / weights_sum;
    }
    // Resample multinomially
    // std::vector<int> indeces(n_particles);
    // std::iota(indeces.begin(), indeces.end(), 0);  // essentially we fill it with 0, 1, ..., n_particles
    std::vector<int> resampled_indeces = sample_with_replacement(weights, n_particles);
    std::vector<std::vector<double> > resampled_particles(t_max+1, std::vector<double>(n_particles));
    for (int j=0; j < n_particles; j++){ // for every particle
      for (int k=0; k < t_max+1; k++){   // for every time step
        resampled_particles[k][j] = particles[k][resampled_indeces[j]];
      }
    }
    particles = resampled_particles;
  }
  // First, return the dataset
  for (int j=0; j <2; j++){
    for (int i=0; i<t_max+1; i++){
      std::cout << data[i][j] << std::endl;
    }
  }
  // Now return the particles
  for (int j=0; j < n_particles; j++){
    for (int i=0; i <t_max +1; i++){
      std::cout << particles[i][j] << std::endl;
    }
  }
}


  
