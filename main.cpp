#include <iostream>

#include "subject.h"

using namespace std;

int main() {
    ////////////// INITIALIZATION //////////////
    const int N_dim = 3;                // N
    const int Ns_dim = 2;               // Ns
    const int number_of_subjects = 10;  // p


    //////////////// PARAMETERS ////////////////
    // SHARED AND UPDATED : p0_hat, t0_hat, v0_hat, (delta_hat)k, (beta_hat)k, sigma_ksi, sigma_tau, sigma
    // FIXED : sigma_p0, sigma_t0, sigma_vo, sigma_delta, sigma_beta
    Parameter P(1,1,1,1,1);
    P.initialize_parameters(N_dim, Ns_dim);

    ////////////// RANDOM VARIABLES ////////////
    auto shared_rv = std::make_shared<Random_Variable>(P);

    /////////////// GEODESIC ///////////////////
    auto shared_geodesic = std::make_shared<Geodesic>(shared_rv);

    ///////////// FAKE SUBJECTS ////////////////
    std::vector<indiv_observation> observations = create_fake_observation(number_of_subjects, N_dim);
    std::vector<Subject> subjects = create_fake_subjects(observations, P, shared_geodesic);


    //////////////// MCMC SAEM /////////////////


    return 0;
}