#include <iostream>
#include <memory>
#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Algorithm/Algorithm.h"
#include "Samplers/HMWithinGibbsSampler.h"
#include "Manifolds/BaseManifold/LogisticBaseManifold.h"
//#include "itkXMLFile.h"

using namespace std;

typedef std::vector< std::vector< std::pair< std::vector<double>, double> > > Data;
typedef std::map<std::string, std::vector<double>> Realizations;


int main() {

    //// INITIALIZATION ///
    unsigned int NumberDimension = 3;
    unsigned int NumberIndependentComponents = 2;


    /// Base Manifold
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();

    /// Manifold
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension, BaseManifold);

    /// Model
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents, Manifold);

    /// Data
    Model->InitializeFakeRandomVariables();
    std::shared_ptr<Data> D = std::make_shared<Data>( Model->SimulateData(10, 3, 3) );

    // Model
    Model->InitializeRandomVariables();

    /// Sampler
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();


    /// Algo
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);

    /// Realizations
    auto R = std::make_shared<Realizations>( Model->SimulateRealizations(5) );


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