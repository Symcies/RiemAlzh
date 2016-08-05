#include <iostream>
#include <memory>
#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Algorithm/Algorithm.h"
#include "Samplers/HMWithinGibbsSampler.h"
//#include "itkXMLFile.h"

using namespace std;

typedef std::vector< std::vector< std::pair< std::vector<double>, double> > > Data;
typedef std::map<std::string, std::vector<double>> Realizations;


int main() {
    /// Manifold
    unsigned int NumberDimension = 5;
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension);
    Manifold->InitializeRandomVariables();

    /// Model
    unsigned int NumberIndependentComponents = 3;
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents);
    Model->SetManifold(Manifold);

    /// Data
    Model->InitializeFakeRandomVariables();
    //shared_ptr<Data> D(Model->SimulateData(15, 2, 5));

    // Model
    Model->InitializeRandomVariables();

    /// Sampler
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();

    ///////////////////////////////
    //// Debugging : Unit tests ///
    shared_ptr<Realizations> R(Model->SimulateRealizations(15));
    double T0 = R->at("T0")[0];
    double P0 = R->at("P0")[0];
    std::cout << "T0/P0 : " << T0 << "/" << P0 << endl;
    std::vector<double> Geodesic = Manifold->GetGeodesic(T0, R);
    std::vector<double> SpaceShift;
    SpaceShift.push_back(0.04);
    SpaceShift.push_back(0.02);
    SpaceShift.push_back(0.15);
    SpaceShift.push_back(0.21);
    SpaceShift.push_back(0.11);
    for(double t = 60.0; t < 73.0; t += 0.5)
    {

        std::vector<double> A = Manifold->ComputeParallelCurve(t, SpaceShift, R);
        double Norm = Manifold->ComputeScalarProduct(A, A, Geodesic);
        std::cout << t << " : " << Norm << std::endl;
    }


    /// Algo
    /*
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);
     */


    cout << "Hello, World!" << endl;
    return 0;

    //
    //                       _oo0oo_
    //                      o8888888o
    //                      88" . "88
    //                      (| -_- |)
    //                      0\  =  /0
    //                    ___/`---'\___
    //                  .' \\|     |// '.
    //                 / \\|||  :  |||// \
    //                / _||||| -:- |||||- \
    //               |   | \\\  -  /// |   |
    //               | \_|  ''\---/''  |_/ |
    //               \  .-\__  '-'  ___/-. /
    //             ___'. .'  /--.--\  `. .'___
    //          ."" '<  `.___\_<|>_/___.' >' "".
    //         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
    //         \  \ `_.   \_ __\ /__ _/   .-` /  /
    //     =====`-.____`.___ \_____/___.-`___.-'=====
    //                       `=---='
    //
    //
    //     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //               佛祖保佑         永无BUG
    //



}