#include <Rcpp.h>
#include <sitmo.h> 
using namespace Rcpp;

// [[Rcpp::depends(sitmo)]]


// [[Rcpp::export]]
NumericVector sample_prior(unsigned int n_particles, double mean, double sd){
  return rnorm(n_particles, mean, sd);
}

// [[Rcpp::export]]
double transition(double x, double t, double transition_scale){
  return R::rnorm(0.5*x + 25*x/(1 + pow(x, 2)) + 8*cos(1.2*t), transition_scale);
}

// [[Rcpp::export]]
double likelihood(double y, double x, double scale){
  return R::dnorm(y, pow(x, 2)/20, scale, FALSE);
}

// [[Rcpp::export]]
double emission(double x, double emission_scale){
  return R::rnorm(pow(x, 2)/20, emission_scale);
}

// [[Rcpp::export]]
NumericMatrix generate_data(int t_max, double transition_scale, double emission_scale, double prior_scale){
  NumericMatrix data(t_max+1, 2);  // First column for x, second for y
  data(0, 0) = R::rnorm(0.0, prior_scale);
  for (unsigned int i=1; i < t_max+1; i++){
    data(i, 0) = transition(data(i-1, 0), i, transition_scale);  // move forward in the latent space
    data(i, 1) = emission(data(i, 0), emission_scale);           // generate corresponding observed value
  }
  return data;
}

// [[Rcpp::export]]
NumericVector integer_sequence(unsigned int first, unsigned int last) {
  NumericVector y(last - first + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}


// [[Rcpp::export]]
NumericMatrix particle_filter(unsigned int t_max, unsigned int n_particles, NumericVector y, double prior_mean, double prior_scale,
                              double transition_scale, double likelihood_scale){
  // Initialize particles and sample x_0
  NumericMatrix particles(t_max+1, n_particles);
  particles.row(0) = sample_prior(n_particles, prior_mean, prior_scale);
  // Instantiate the weights vector (will be overwritten)
  NumericVector weights(n_particles);
  // Bootstrap Filter for every time step
  for (unsigned int i=1; i < t_max + 1; i++){
    // Move forward the particles
    for (unsigned int particle_index=0; particle_index < n_particles; particle_index ++){
      particles(i, particle_index) = transition(particles(i-1, particle_index), i, transition_scale);
    }
    
    // Compute importance weights
    double weights_sum = 0.0;
    for (unsigned int particle_index=0; particle_index < n_particles; particle_index++){
      weights(particle_index) = likelihood(y(i), particles(i, particle_index), likelihood_scale);
      weights_sum += weights(particle_index);
    }
    
    // Normalize importance weights
    for (unsigned int particle_index=0; particle_index < n_particles; particle_index++){
      weights(particle_index) = weights(particle_index) / weights_sum;
    }
    // Resample
    NumericVector indices = sample(integer_sequence(0, n_particles-1), n_particles, true, weights);
    NumericMatrix resampled_particles(t_max+1, n_particles);
    for (unsigned int ix=0; ix < n_particles; ix ++){
      resampled_particles.column(ix) = particles.column(indices[ix]);
    }
    for (unsigned int col=0; col < n_particles; col++){
      particles.column(col) = resampled_particles.column(col);
    }
  }
  // Return the particles
  return particles;
}
  
  
  
  
  
  
  
  