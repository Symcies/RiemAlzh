#include <iostream>
#include <memory>
#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Algorithm/Algorithm.h"

using namespace std;

typedef std::vector< std::vector< std::pair< std::vector<double>, unsigned int> > > Data;

int main() {
    /// Manifold
    unsigned int NumberDimension = 5;
    unsigned int NumberIndependentComponents = 3;
    //shared_ptr<AbstractManifold> Manifold(new PropagationManifold(NumberDimension, NumberIndependentComponents));
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension);
    Manifold->InitializeRandomVariables();

    /// Model
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents);
    Model->SetManifold(Manifold);
    Model->InitializeRandomVariables();

    /// Data
    shared_ptr<Data> D = make_shared<Data>();

    /// Sampler
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbs>();

    /// Algo
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);


    cout << "Hello, World!" << endl;
    return 0;
}