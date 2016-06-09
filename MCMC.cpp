
/*

void simulation_step(std::vector<std::vector<double>>& z, std::vector<Subject> S) {


}

std::vector<std::vector<double>> MCMC_SAEM(std::vector<Subject> S) {
    // INITIALIZATION
    std::vector<std::vector<double>> theta;
    std::vector<std::vector<double>> z; // <-- RANDOM (A PRIORI?)
    // S <- 0
    // epsilon(t)
    int max_nb_of_simulation = 10;

    // SIMULATION
    for(int i = 0; i< max_nb_of_simulation; ++i){
        // CHECK CONVERGENCE OR NOT

        // SIMULATION STEP
        simulation_step(z, S, theta);

        // COMPUTE THE SUFFICIENT STATISTICS

        // STOCHASTIC APPROXIMATION STEP

        // MAXIMIZATION STEP

    }

    return theta;
}
*/