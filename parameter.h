#pragma once

#include<vector>
#include "tools.h"

class Parameter{
public:
    Parameter(double position_var, double velocity_var, double time_var, double propagation_time_var, double beta_var);
    void initialize_parameters(int N_dim, int Ns_dim);

    int get_N_dim();
    int get_Ns_dim();
    double get_position_mean();
    double get_position_var();
    double get_time_mean();
    double get_time_var();
    double get_velocity_mean();
    double get_velocity_var();
    std::vector<double> get_propagation_timing_mean();
    double get_propagation_timing_var();
    std::vector<double> get_beta_mean();
    double get_beta_var();
    double get_preacceleration_factor_var();
    double get_time_shift_var();
    double get_uncertainty_var();


// PARAMETERS THAT ARE ESTIMATED BY THE MODEL : POSSIBLE CHANGE
// tau = (...) as described at the end of page 5 of NIPS
private:


    double position_mean;
    double time_mean;
    double velocity_mean;
    std::vector<double> propagation_timing_mean;
    std::vector<double> beta_mean;
    double preacceleration_factor_variance;
    double time_shift_variance;
    double uncertainties_variance;


// PARAMETERS THAT ARE FIXED : NO CHANGE POSSIBLE
private:
    int N_dimension;
    int Ns_dimension;

    double position_variance;
    double velocity_variance;
    double time_variance;
    double propagation_timing_variance;
    double beta_variance;

};