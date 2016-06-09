#include "parameter.h"


Parameter::Parameter(double position_var, double velocity_var, double time_var, double propagation_time_var, double beta_var) {
    position_variance = position_var;
    velocity_variance = velocity_var;
    time_variance = time_var;
    propagation_timing_variance = propagation_time_var;
    beta_variance = beta_var;
}

void Parameter::initialize_parameters(int N_dim, int Ns_dim) {
    // TO BE CHANGED : HOW TO INITIALIZE theta(0)
    N_dimension = N_dim;
    Ns_dimension = Ns_dim;

    position_mean = one_normal_distribution(1,1);
    time_mean = one_normal_distribution(1,1);
    velocity_mean = one_normal_distribution(1,1);
    propagation_timing_mean = N_normal_distribution(N_dim-1, 1,1);
    beta_mean = N_normal_distribution( (N_dim-1)*Ns_dim , 1,1);
    preacceleration_factor_variance = one_normal_distribution(1,1);
    time_shift_variance = one_normal_distribution(1,1);
    uncertainties_variance = one_normal_distribution(1,1);
}

int Parameter::get_N_dim() { return N_dimension; }

int Parameter::get_Ns_dim() { return Ns_dimension; }

double Parameter::get_position_mean() { return position_mean; }

double Parameter::get_position_var() { return position_variance; }

double Parameter::get_time_mean() { return time_mean; }

double Parameter::get_time_var() { return time_variance; }

double Parameter::get_velocity_mean() { return velocity_mean; }

double Parameter::get_velocity_var() { return velocity_variance; }

std::vector<double> Parameter::get_propagation_timing_mean() { return  propagation_timing_mean; }

double Parameter::get_propagation_timing_var() { return propagation_timing_variance; }

std::vector<double> Parameter::get_beta_mean() { return beta_mean; }

double Parameter::get_beta_var() { return beta_variance; }

double Parameter::get_preacceleration_factor_var() { return preacceleration_factor_variance; }

double Parameter::get_time_shift_var() { return time_shift_variance; }

double Parameter::get_uncertainty_var() { return uncertainties_variance; }