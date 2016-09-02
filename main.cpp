#include <iostream>
#include <fstream>
#include <memory>


#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Algorithm/Algorithm.h"
#include "Samplers/HMWithinGibbsSampler.h"
#include "Manifolds/BaseManifold/LogisticBaseManifold.h"
#include "Outputs/RandomVariableRealizations.h"
//#include "itkXMLFile.h"

using namespace std;

typedef vector< vector< pair< vector<double>, double> > > Data;
typedef map<std::string, vector<double>> Realizations;


int main() {

    //// INITIALIZATION ///
    unsigned int NumberDimension = 3;
    unsigned int NumberIndependentComponents = 2;


    /// Base Manifold, Manifold, Model & Sampler ///
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension, BaseManifold);
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents, Manifold);
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();

    /// DATA GENERATION ///
    Model->InitializeFakeRandomVariables();
    shared_ptr<Data> D = make_shared<Data>( Model->SimulateData(100, 3, 4) );

    // Model
    Model->InitializeRandomVariables();

    /// Algo
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);

    /// Realizations
    //auto R = std::make_shared<Realizations>( Model->SimulateRealizations(5) );

    /// Realizations evolution
    map< string, vector<vector<double>> > Real = Algo->GetRealizationEvolution();
    ofstream File;
    File.open("Realizations.txt");
    RealizationsEvolution(Real, File);


    cout <<  endl << "Hello, World!" << endl;
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