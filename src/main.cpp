#include <iostream>
#include <fstream>
#include <memory>

typedef double ScalarType;

#include "Algorithm.h"

#include "ExponentialCurveManifold.h"
#include "PropagationManifold.h"
#include "LogisticBaseManifold.h"
#include "LinearManifold.h"

#include "LongitudinalModel.h"
#include "UnivariateModel.h"
#include "FastNetworkModel.h"
#include "NetworkPropagationModel.h"
#include "NetworkPropagationModel2.h"
#include "TestModel.h"

#include "BlockedGibbsSampler.h"

#include "LinearAlgebra.h"

#include "ReadData.h"
#include "AlgorithmSettings.h"

#include "omp.h"
#include "tinyxml2.h"


using namespace std;

typedef vector< vector< pair< LinearAlgebra<ScalarType>::VectorType, double> > > Data;

// TODO : Change Model->UpdateModel because it ain't parameters
// TODO : Finish the input classes : Model, Algo and data
// TODO : Create an output file if it does not exist
// TO ADD : Name of the output file
// TO ADD : If writing a file 1 out of n iteration - or only at the end 


int main(int argc, char* argv[]) {
    
    if(argc != 4) {
        std::cout << "Usage: " << " /path/to/executable " << " model_settings.xml " << " algorithm_settings " << "data_settings.xml" << std::endl;
        // return 1;
    }
    
        
    
    ///////////////////////
    //// Initialization ///
    ///////////////////////
    
    //TestAssert::Init(false);
    std::string ModelType = "FastNetwork";
    
    std::string FilePath;
    std::shared_ptr<AbstractModel> Model;
    std::shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    
    /////////////////////////////
    //// Initialize the model ///
    /////////////////////////////
    
    if(ModelType == "Test")
    {
        Model = make_shared<TestModel>();

    }
    else if(ModelType == "Univariate")
    {
        shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
        Model = make_shared<UnivariateModel>(BaseManifold);
    }
    else if(ModelType == "Multivariate")
    {
        FilePath = "/Users/igor.koval/Documents/Work/RiemAlzh/data/CognitiveScores/SimulatedData/";
        
        unsigned int NumberDimensions = 4, NbIndependentComponents = 2;
        
        shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
        shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimensions, BaseManifold);
        Model = make_shared<LongitudinalModel>(NbIndependentComponents, Manifold);
      
    }
    else if(ModelType == "FastNetwork")
    {
        FilePath = "/Users/igor.koval/Documents/Work/RiemAlzh/data/CorticalThickness/MCIconvertAD/MCIconvertAD_";
        
        unsigned int NbIndependentComponents = 5;
        
        std::string KernelMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/data/invKd_16.csv");
        auto KernelMatrix = std::make_shared<FastNetworkModel::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
        std::string InterpolationMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/data/Kxd_16.csv");
        auto InterpolationMatrix = std::make_shared<FastNetworkModel::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
        
        Model = make_shared<FastNetworkModel>(NbIndependentComponents, KernelMatrix, InterpolationMatrix);
    }
    else if (ModelType == "Network")
    {
        FilePath = "/Users/igor.koval/Documents/Work/RiemAlzh/data/CorticalThickness/MCIconvertAD/MCIconvertAD_";
        
        unsigned int NbIndependentComponents = 5;
        
        std::string KernelMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/data/invKd_16.csv");
        auto KernelMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
        std::string InterpolationMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/data/Kxd_16.csv");
        auto InterpolationMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
        
        Model = make_shared<NetworkPropagationModel>(NbIndependentComponents, KernelMatrix, InterpolationMatrix);
    }
    else if(ModelType == "Network2")
    {
        FilePath = "/Users/igor.koval/Documents/Work/RiemAlzh/datatest/CorticalThickness/MCIconvertAD/MCIconvertAD_";
        
        unsigned int NbIndependentComponents = 10;
                
        std::string KernelMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/data/invKd_16.csv");
        auto KernelMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
        std::string InterpolationMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/data/Kxd_16.csv");
        auto InterpolationMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
        
        shared_ptr<AbstractManifold> Manifold = make_shared<ExponentialCurveManifold>(InterpolationMatrix->columns());
        
        Model = make_shared<NetworkPropagationModel2>(NbIndependentComponents, Manifold, KernelMatrix, InterpolationMatrix);
    }
    
    
    /////////////////////////////
    /// Read or Simulate Data ///
    /////////////////////////////
    
    Data D;
    bool ReadData = true;
    if(ReadData)
    {
        int NbMaxOfSubjects = 100;
        D = ReadData::OpenFilesMultivariate(FilePath, NbMaxOfSubjects);
    }
    else
    {
        Model->InitializeFakeRandomVariables();
        D = Model->SimulateData(300, 4, 6);
    }
    
    
    
    if(argc == 3)
    {
        Model->InitializeFakeRandomVariables();
        D = Model->SimulateData(300, 4, 6);
    }
    if(argc == 5)
    {
        int NbMaxOfSubjects = 100;
        D = ReadData::OpenFilesMultivariate(FilePath, NbMaxOfSubjects);
    }
     
    
    
    //////////////////////////
    /// Algorithm pipeline ///
    //////////////////////////
    AlgorithmSettings AS(argv[2]);
    auto Algo = make_shared<Algorithm>(AS);
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);
    
    
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
