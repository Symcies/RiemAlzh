#include <iostream>
#include <fstream>
#include <memory>


#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Algorithm/Algorithm.h"
#include "Samplers/HMWithinGibbsSampler.h"
#include "Samplers/BlockedGibbsSampler.h"
#include "Manifolds/BaseManifold/LogisticBaseManifold.h"
#include "Tests/TestAssert.h"
#include "Outputs/RandomVariableRealizations.h"
//#include "itkXMLFile.h"


using namespace std;

typedef vector< vector< pair< vector<double>, double> > > Data;
typedef map<std::string, vector<double>> Realizations;


int main() {

    //// INITIALIZATION ///
    unsigned int NumberDimension = 3;
    unsigned int NumberIndependentComponents = 2;
    
    /// TESTS ///
    bool Active = true;
    TestAssert::Init(Active);


    /// Base Manifold, Manifold, Model & Sampler ///
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension, BaseManifold);
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents, Manifold);
    // shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();

    /// DATA GENERATION ///
    Model->InitializeFakeRandomVariables();
    shared_ptr<Data> D = make_shared<Data>( Model->SimulateData(100, 4, 5) );

    /// Model
    Model->InitializeRandomVariables();


    /// Python call
    std::string filename = "/Users/igor.koval/PycharmProjects/RiemAlzh/plotGraph.py";
    std::string command = "python ";
    command += filename;
    //system(command.c_str());

    /// Algo
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);



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