#pragma once

#include <vector>
#include <random>
#include <array>
#include <math.h>
#include <algorithm>
#include <complex>
#include <memory>
#include <iostream> //cout
#include "parameter.h"
#include "tools.h"
#include "random_variable.h"
#include "geodesic.h"

// observations in R(N) and time tij
typedef std::vector<std::vector<double>> indiv_observation;

class Subject{
public:
    Subject(int subject_ID, indiv_observation observations, Parameter P, std::shared_ptr<Geodesic> G);
    ~Subject();

    std::vector<double> get_transport_vector(double t); // Here, the time_warp is already taken into account

private:
    int ID;
    int N_dimension;
    int Ns_dimension;
    size_t nb_of_observation;

    // Common Random variables
    double initial_time;
    std::vector<double> propagation_time;


    // Subject-specific random Variables
    double accelerator_factor;            // alpha (proportional to ksi)
    double time_shift;                    // tau
    std::vector<double> laplace_distrib;  // s(i;j)
    //std::vector<double, double> w_k_i;

    // Observation
    indiv_observation observation;

    // Pointer to the geodesic
    std::shared_ptr<Geodesic> geodesic;
};

std::vector<indiv_observation> create_fake_observation(int nb_subjects, int N_dim);

std::vector<Subject> create_fake_subjects(std::vector<indiv_observation> observations, Parameter P, std::shared_ptr<Geodesic> shared_geodesic);