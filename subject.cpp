#include "subject.h"


Subject::Subject(int subject_ID, indiv_observation observations, Parameter P, std::shared_ptr<Geodesic> G) {
    ID = subject_ID;
    N_dimension = P.get_N_dim();
    Ns_dimension = P.get_Ns_dim();
    nb_of_observation = observations.size();

    initial_time = G->get_random_var()->get_initial_time();
    propagation_time = G->get_random_var()->get_propagation_timing();

    double pre_accelerator_factor = one_normal_distribution(0, P.get_preacceleration_factor_var());
    accelerator_factor = exp(pre_accelerator_factor);
    time_shift = one_normal_distribution(0, P.get_time_shift_var());
    laplace_distrib = N_laplace_distribution(Ns_dimension, 0, 1/2);

    observation = observations;
    geodesic = G;
}

Subject::~Subject() {

}

std::vector<double> Subject::get_transport_vector(double t) {
    std::vector<double> transport_vector;
    double FAKE_W_i = 0.5;
    double t0 = geodesic->get_random_var()->get_initial_time();
    double constant_part = accelerator_factor *  (t - t0 - time_shift) + t0;

    for(int k = 0; k<N_dimension; ++k) {
        double delta = propagation_time[k];
        double k_part = FAKE_W_i / geodesic->get_geodesic_derivative_k_at_t(k, t + delta) + delta;
        double time_point = constant_part + k_part;
        double geo_value = geodesic->get_geodesic_k_at_t(k, time_point);
        transport_vector.push_back(geo_value);
    }
    return transport_vector;
}

std::vector<indiv_observation> create_fake_observation(int nb_subjects, int N_dim) {
    std::vector<indiv_observation> observations;

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<int> nb_observ(2,10);


    for(int i=0; i<nb_subjects ; ++i) {
        int nb_observation = nb_observ(generator);
        indiv_observation observation;

        for(int j=0; j<nb_observation ; ++j) {
            std::vector<double> obs;
            for(int k = 0; k<N_dim ; ++k) {
                // TO BE CHANGED
                double coordinate = i+k;
                obs.push_back(coordinate);
            }
        observation.push_back(obs);
        }
        observations.push_back(observation);
    }
    return observations;
}

std::vector<Subject> create_fake_subjects(std::vector<indiv_observation> observations, Parameter P, std::shared_ptr<Geodesic> shared_geodesic) {
    size_t number_of_subjects = observations.size();
    std::vector<Subject> subjects;

    for(int i=0; i<number_of_subjects; ++i) {
        Subject S = Subject(i, observations[i], P, shared_geodesic);
        subjects.push_back(S);
    }

    return subjects;
}
