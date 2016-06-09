#pragma once

#include <vector>
#include <cmath>
#include <random>
#include <math.h>
#include <algorithm>
#include <complex>


double l2_norm(std::vector<double> vect);

double one_normal_distribution(double mean, double std_dev);

std::vector<double> N_normal_distribution(int number_of_draws, double mean, double std_dev);

std::vector<double> N_laplace_distribution(int number_of_draws, int mean, int variance);
