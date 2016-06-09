#include "random_variable.h"

Random_Variable::Random_Variable(Parameter P) {
    N_dimension = P.get_N_dim();
    Ns_dimension = P.get_Ns_dim();

    initial_position = one_normal_distribution( P.get_position_mean(), P.get_position_var());
    initial_time = one_normal_distribution(P.get_time_mean(), P.get_time_var());
    initial_velocity = one_normal_distribution(P.get_velocity_mean(), P.get_velocity_var());

    std::vector<double> prop_times;
    for(int i = 0; i < N_dimension-1; ++i) {
        double time = one_normal_distribution(P.get_propagation_timing_mean()[i], P.get_propagation_timing_var());
        prop_times.push_back(time);
    }
    propagation_timing = prop_times;

    std::vector<double> betas;
    for(int i = 0; i< (N_dimension-1)*Ns_dimension ; ++i ) {
        double bet = one_normal_distribution(P.get_beta_mean()[i], P.get_beta_var());
        betas.push_back(bet);
    }
    beta = betas;

}

int Random_Variable::get_Ns_dim() { return Ns_dimension; }

int Random_Variable::get_N_dim() { return N_dimension; }

double Random_Variable::get_initial_position() { return initial_position; }

double Random_Variable::get_initial_time() { return initial_time; }

double Random_Variable::get_initial_velocity() { return initial_velocity; }

std::vector<double> Random_Variable::get_propagation_timing() { return propagation_timing; }

std::vector<double> Random_Variable::get_beta() { return beta; };