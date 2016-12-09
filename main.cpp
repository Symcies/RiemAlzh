#include <iostream>
#include <fstream>
#include <memory>

typedef double ScalarType;

#include "Algorithm/Algorithm.h"
#include "Manifolds/ExponentialCurveManifold.h"
#include "Manifolds/PropagationManifold.h"
#include "Manifolds/BaseManifold/LogisticBaseManifold.h"
#include "Models/LongitudinalModel.h"
#include "Models/UnivariateModel.h"
#include "Models/NetworkPropagationModel2.h"
#include "Models/TestModel.h"
#include "Samplers/BlockedGibbsSampler.h"
#include "Tests/TestAssert.h"
#include "Outputs/RandomVariableRealizations.h"
#include "LinearAlgebra/LinearAlgebra.h"
#include "Inputs/ReadData.h"

#include "omp.h"

//#include "itkXMLFile.h"


using namespace std;

typedef vector< vector< pair< LinearAlgebra<ScalarType>::VectorType, double> > > Data;
typedef map<std::string, vector<double>> Realizations;



int main() {
    // TODO : Remplacer les const 'std::shared<T>&' par 'std::shared<const T>'
    // TODO : Change Model->UpdateParameters because it ain't parameters
    

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
    //unsigned int NumberDimension = 10;
    unsigned int NumberIndependentComponents = 2;
    clock_t start = clock();
    
    /////////////
    /// Tests ///
    /////////////
    bool Active = true;
    TestAssert::Init(Active);
    
    /////////////////////////////////
    /// Network Propagation Model ///
    /////////////////////////////////
        
    /// Open the files - toy model
    /*
    std::string KernelMatrixToyPath ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Kd_toyexample.csv");
    std::string InterpolationMatrixToyPath ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Kxd_toyexample.csv");
    shared_ptr<NetworkPropagationModel::MatrixType> KernelMatrix = std::make_shared<NetworkPropagationModel::MatrixType>(ReadData::OpenKernel(KernelMatrixToyPath));
    shared_ptr<NetworkPropagationModel::MatrixType> InterpolationMatrix = std::make_shared<NetworkPropagationModel::MatrixType>(ReadData::OpenKernel(InterpolationMatrixToyPath));
    */
     
    /// Open the files - real example
    
    std::string KernelMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/DataCorticalThickness/invKd.csv");
    auto KernelMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
    std::string InterpolationMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/DataCorticalThickness/Kxd.csv");
    auto InterpolationMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
    shared_ptr<Data> D = std::make_shared<Data>(ReadData::OpenFilesMultivariate());
    
    //TODO : Check the data 
    // Read the initializations
    
    /// Initiate the Manifolds and the model
    shared_ptr<AbstractManifold> Manifold = make_shared<ExponentialCurveManifold>(InterpolationMatrix->rows());
    shared_ptr<AbstractModel> Model = make_shared<NetworkPropagationModel2>(NumberIndependentComponents, Manifold, KernelMatrix, InterpolationMatrix);
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    //Model->InitializeFakeRandomVariables();
    //shared_ptr<Data> D = make_shared<Data>( Model->SimulateData(300, 4, 6) );
    
    
    //////////////////////////
    /// Multivariate Model ///
    //////////////////////////
    /// Read the data 
    /*
    //shared_ptr<Data> D = std::make_shared<Data>(ReadData::OpenFilesMultivariate());
    
    /// Initialize the manifolds, model and sampler 
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension, BaseManifold);
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents, Manifold);
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    */
     
    ////////////////////////
    /// Univariate Model ///
    ////////////////////////
    /*
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractModel> Model = make_shared<UnivariateModel>(BaseManifold);
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    */
    
    //////////////////
    /// Test Model ///
    //////////////////
    /*
    shared_ptr<AbstractModel> Model = make_shared<TestModel>();
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    */
    
    ////////////////////////////////////////
    /// Data Generation & Initialization ///
    ////////////////////////////////////////
    
    
    //////////////////////////
    /// Algorithm pipeline ///
    //////////////////////////
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);
    auto Parameters = Model->GetParameters();
    
    
    for(auto it = Parameters.begin(); it != Parameters.end(); ++it)
    {
        std::cout << it->first << " --> " << it->second << std::endl;
    }



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
