#include <iostream>
#include <fstream>
#include <memory>

typedef double ScalarType;

#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Models/UnivariateModel.h"
#include "Models/NetworkPropagationModel.h"
#include "Algorithm/Algorithm.h"
#include "Samplers/HMWithinGibbsSampler.h"
#include "Samplers/BlockedGibbsSampler.h"
#include "Manifolds/BaseManifold/LogisticBaseManifold.h"
#include "Tests/TestAssert.h"
#include "Outputs/RandomVariableRealizations.h"
#include "LinearAlgebra/LinearAlgebra.h"
#include "Inputs/ReadData.h"

//#include "itkXMLFile.h"


using namespace std;

typedef vector< vector< pair< LinearAlgebra<ScalarType>::VectorType, double> > > Data;
typedef map<std::string, vector<double>> Realizations;



int main() {

    ///////////////////
    /// Python call ///
    ///////////////////
    std::string filename = "/Users/igor.koval/PycharmProjects/RiemAlzh/plotGraph.py";
    std::string command = "python ";
    command += filename;
    //system(command.c_str());
    
    
    ///////////////////////
    //// Initialization ///
    ///////////////////////
    unsigned int NumberDimension = 4;
    unsigned int NumberIndependentComponents = 1;
    clock_t start = clock();
    
    /////////////
    /// Tests ///
    /////////////
    bool Active = true;
    TestAssert::Init(Active);
    
    /////////////////////////
    /// Propagation Model ///
    /////////////////////////
    
    /// Open the files
    std::string KernelMatrixPath ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Kd_toyexample.csv");
    std::string InterpolationMatrixPath ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Kxd_toyexample.csv");
    shared_ptr<NetworkPropagationModel::MatrixType> KernelMatrix = std::make_shared<NetworkPropagationModel::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
    shared_ptr<NetworkPropagationModel::MatrixType> InterpolationMatrix = std::make_shared<NetworkPropagationModel::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
    
    /// Initiate the Manifolds and the model
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(KernelMatrix->rows(), BaseManifold);
    shared_ptr<AbstractModel> Model = make_shared<NetworkPropagationModel>(NumberIndependentComponents, Manifold, KernelMatrix, InterpolationMatrix);
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();
    
    //////////////////////////
    /// Multivariate Model ///
    //////////////////////////
    /*
     * 
    /// Read the data 
    shared_ptr<Data> D = std::make_shared<Data>(ReadData::OpenFilesMultivariate());
     
    /// Initialize the manidolds, model and sampler 
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension, BaseManifold);
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents, Manifold);
    shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();
    //shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    */
     
    ////////////////////////
    /// Univariate Model ///
    ////////////////////////
    /*
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractModel> Model = make_shared<UnivariateModel>(BaseManifold);
    //shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    */
    
    
    ////////////////////////////////////////
    /// Data Generation & Initialization ///
    ////////////////////////////////////////
    
    Model->InitializeFakeRandomVariables();
    shared_ptr<Data> D = make_shared<Data>( Model->SimulateData(200, 6, 8) );
    Model->Initialize();
    
    
    //////////////////////////
    /// Algorithm pipeline ///
    //////////////////////////
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);



    clock_t end = clock();
    cout <<  endl << "Time : " << (end - start) / CLOCKS_PER_SEC << " secondes" << endl;
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
