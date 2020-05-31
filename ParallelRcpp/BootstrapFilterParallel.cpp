#include <Rcpp.h>
#include <sitmo.h> 
#include <RcppParallel.h>
#include <omp.h>
#include <boost/random/normal_distribution.hpp>
#include <random>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(sitmo)]]


int random_int(const int from_integer, const int to_integer){
  std::random_device                  my_device;
  std::mt19937                        my_generator(my_device());
  std::uniform_int_distribution<int>  my_distribution(from_integer, to_integer);
  return my_distribution(my_generator);
}

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

// Computes mean of the transition function
double meanfunc(double x, unsigned int t){
  return 0.5*x + 25*x/(1 + pow(x, 2)) + 8*cos(1.2*t);
}


// [[Rcpp::export]]
NumericMatrix particle_filter_par(unsigned int t_max, unsigned int n_particles, NumericVector y, double prior_mean, double prior_scale,
                              double transition_scale, double likelihood_scale, unsigned int ncores){
  // Initialize particles and sample x_0
  NumericMatrix particles(t_max+1, n_particles);
  particles.row(0) = sample_prior(n_particles, prior_mean, prior_scale);
  // Initialize temporary container for resampled particles
  NumericMatrix resampled_particles(t_max+1, n_particles);
  // Parallel wrappers for particles and resampled particles
  RcppParallel::RMatrix<double> Pparticles(particles);
  // Instantiate the weights vector (will be overwritten)
  NumericVector weights(n_particles);
  // Set the RNG scope (there's a bug otherwise)
  RNGScope scope;
  // Generate DIFFERENT random seeds, which will generate different engines. These engines will be fed into the distributions
  // from the boost library 
  NumericVector seeds(ncores);
  for (unsigned int ii=0; ii < ncores; ii++){
    seeds[ii] = random_int(100*ii + 1, 100*ii + 100);
  }
  // Bootstrap Filter for every time step
  for (unsigned int i=1; i < t_max + 1; i++){
    // Get the seed and create the engine
    uint32_t coreSeed = static_cast<uint32_t>(seeds[omp_get_thread_num()]);   // grab the corresponding seed and cast it
    sitmo::prng_engine coreEngine = sitmo::prng_engine(coreSeed);    // grab the engine
    // Move forward the particles
    #pragma omp parallel num_threads(ncores)
    {
      #pragma omp for schedule(static)
      for (unsigned int particle_index=0; particle_index < n_particles; particle_index ++){
        // instantiate the distribution
        boost::random::normal_distribution<float> trans(meanfunc(Pparticles(i-1, particle_index), i), transition_scale);
        //particles(i, particle_index) = transition(particles(i-1, particle_index), i, transition_scale); //trans(coreEngine);
        particles(i, particle_index) = trans(coreEngine);
      }
    }
    // PARALLEL : Compute importance weights
    double weights_sum = 0.0;
    #pragma omp parallel num_threads(ncores)
    {
      #pragma omp for schedule(static)
      for (unsigned int particle_index=0; particle_index < n_particles; particle_index++){
        weights(particle_index) = likelihood(y(i), particles(i, particle_index), likelihood_scale);
        weights_sum += weights(particle_index);
      }
    }
    // PARALLEL : normalize importance weights
    #pragma omp parallel num_threads(ncores)
    {
      #pragma omp for schedule(static)
      for (unsigned int particle_index=0; particle_index < n_particles; particle_index++){
        weights(particle_index) = weights(particle_index) / weights_sum;
      }
    }
    // SEQUENTIAL : Resample indices
    NumericVector indices = sample(integer_sequence(0, n_particles-1), n_particles, true, weights);
    // PARALLEL : Select resampled columns
    #pragma omp parallel num_threads(ncores)
    {
      #pragma omp for schedule(static)
      for (unsigned int ix=0; ix < n_particles; ix ++){
        resampled_particles.column(ix) = particles.column(indices[ix]);
      }
    }
    // PARALLEL : Reassigning resampled particles to particles
    #pragma omp parallel num_threads(ncores)
    {
      #pragma omp for schedule(static)
      for (unsigned int col=0; col < n_particles; col++){
        particles.column(col) = resampled_particles.column(col);
      }
    }
  }
  // Return the particles
  return particles;
}







