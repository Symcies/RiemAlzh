#include "geodesic.h"


Geodesic::Geodesic(ptr_rv Random_Var) {
    random_var  = Random_Var;
    N_dimension = Random_Var->get_N_dim();


    std::vector<geodesic0> Geodesic;
    std::vector<geodesic0> Geodesic_derivative;

    geodesic0 geo0 = create_geodesic;
    geodesic0 geo0_derivative = create_geodesic_derivative;

    for(int i =0; i<N_dimension; ++i) {
        Geodesic.push_back(geo0);
        Geodesic_derivative.push_back(geo0_derivative);
    }
    geodesic = Geodesic;
    geodesic_derivative = Geodesic_derivative;

}

Geodesic::~Geodesic() {

}

ptr_rv Geodesic::get_random_var() { return random_var; }

double Geodesic::get_geodesic_k_at_t(int k, double t) {
    return geodesic[k](t, random_var);

}

double Geodesic::get_geodesic_derivative_k_at_t(int k, double t) {
    return geodesic_derivative[k](t, random_var);
}



double create_geodesic(double t, ptr_rv rv) {
    double p0 = rv->get_initial_position();
    double t0 = rv->get_initial_time();
    double v0 = rv->get_initial_velocity();

    double geodesic = (1 + (1.0/p0 - 1) * exp( - v0 / (p0*(1-p0)) * ( t - t0)));
    return geodesic;
}

double create_geodesic_derivative(double t, ptr_rv rv) {
    double p0 = rv->get_initial_position();
    double t0 = rv->get_initial_time();
    double v0 = rv->get_initial_velocity();
    double exp_value = exp(-v0 / (p0 * (1-p0)) * (t - t0));
    double denom = (p0 + (1-p0) * exp_value) * (p0 + (1-p0) * exp_value);

    return v0 * exp_value / denom;
}


