#pragma once

#include <vector>
#include "parameter.h"
#include "tools.h"

class Random_Variable {

public:
    Random_Variable(Parameter P);
    int get_Ns_dim();
    int get_N_dim();
    double get_initial_position();
    double get_initial_time();
    double get_initial_velocity();
    std::vector<double> get_propagation_timing();
    std::vector<double> get_beta();

private:
    int N_dimension;
    int Ns_dimension;
    double initial_position;
    double initial_time;
    double initial_velocity;
    std::vector<double> propagation_timing;
    std::vector<double> beta;

};