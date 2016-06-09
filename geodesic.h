#pragma once

#include <iostream>
#include <memory>
#include <bits/shared_ptr.h>
#include "random_variable.h"

typedef std::shared_ptr<Random_Variable> ptr_rv;
typedef double (*geodesic0)(double t, ptr_rv random_var);


class Geodesic {

public:
    Geodesic(std::shared_ptr<Random_Variable> Random_Var);

    ~Geodesic();

    ptr_rv get_random_var();
    double get_geodesic_k_at_t(int k, double t);
    double get_geodesic_derivative_k_at_t(int k, double t);



private:
    double N_dimension;
    ptr_rv random_var;
    std::vector<geodesic0> geodesic;
    std::vector<geodesic0> geodesic_derivative;
};


double create_geodesic(double t, ptr_rv random_var);

double create_geodesic_derivative(double t, ptr_rv random_var);


