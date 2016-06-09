#include "tools.h"

double l2_norm(std::vector<double> vect) {
    double accum = 0;
    for(int i = 0; i < vect.size() ; ++i) {
        accum += vect[i]*vect[i];
    }
    return sqrt(accum);
}

double one_normal_distribution(double mean, double std_dev) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(mean, std_dev);

    return distribution(generator);
}

std::vector<double> N_normal_distribution(int number_of_draws, double mean, double std_dev) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::normal_distribution<double> distribution(mean, std_dev);

    std::vector<double> draws;

    for(int i=0 ; i<number_of_draws ; ++i) {
        double number = distribution(generator);
        draws.push_back(number);
    }

    return draws;
}

std::vector<double> N_laplace_distribution(int number_of_draws, int mean, int variance) {
    std::vector<double> draws;
    for(int i =0; i< number_of_draws; ++i) {
        // TO BE CHANGED TO FIND THE LAPLACE DISTRIBUTION AND UNDERSTAND WHAT IT MEANS THE 1/2, page 4 of NIPS
        double draw = 0.0;
        draws.push_back(draw);
    }
    return draws;
}

