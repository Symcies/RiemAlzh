#include <iostream>
#include <memory>
#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Algorithm/Algorithm.h"
#include "Samplers/HMWithinGibbsSampler.h"
//#include "itkXMLFile.h"

using namespace std;

typedef std::vector< std::vector< std::pair< std::vector<double>, double> > > Data;

int main() {
    /// Manifold
    unsigned int NumberDimension = 5;
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension);
    Manifold->InitializeRandomVariables();

    /// Model
    unsigned int NumberIndependentComponents = 3;
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents);
    Model->SetManifold(Manifold);
    Model->InitializeFakeRandomVariables();

    /// Data
    shared_ptr<Data> D(Model->SimulateData(15, 2, 5));
    Model->InitializeRandomVariables();

    /// Sampler
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();

    /// Algo
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);


    cout << "Hello, World!" << endl;
    return 0;
}