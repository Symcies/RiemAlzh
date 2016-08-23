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
    unsigned int NumberDimension = 4;
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension);

    /// Model
    unsigned int NumberIndependentComponents = 2;
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents);
    Model->SetManifold(Manifold);

    /// Data
    Model->InitializeFakeRandomVariables();
    shared_ptr<Data> D(Model->SimulateData(1, 2, 3));

    // Model
    Model->InitializeRandomVariables();

    /// Sampler
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();

    /// Realizations
    shared_ptr<Realizations> R(Model->SimulateRealizations(15));

    /*
    /// Algo
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