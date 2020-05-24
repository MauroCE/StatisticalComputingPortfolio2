#ifndef RANDOM_H
#define RANDOM_H

double random_uniform01();
double random_normal(double mean, double scale);
std::vector<int> sample_with_replacement(std::vector<double> probabilities, int n);

#endif