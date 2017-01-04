#include <iostream>
#include <fstream>
#include <memory>

typedef double ScalarType;

#include "Algorithm/Algorithm.h"
#include "Manifolds/ExponentialCurveManifold.h"
#include "Manifolds/PropagationManifold.h"
#include "Manifolds/BaseManifold/LogisticBaseManifold.h"
#include "Manifolds/LinearManifold.h"
#include "Models/LongitudinalModel.h"
#include "Models/UnivariateModel.h"
#include "Models/NetworkPropagationModel.h"
#include "Models/NetworkPropagationModel2.h"
#include "Models/TestModel.h"
#include "Samplers/BlockedGibbsSampler.h"
#include "Tests/TestAssert.h"
#include "Outputs/RandomVariableRealizations.h"
#include "LinearAlgebra/LinearAlgebra.h"
#include "Inputs/ReadData.h"

#include "omp.h"


using namespace std;

typedef vector< vector< pair< LinearAlgebra<ScalarType>::VectorType, double> > > Data;
typedef map<std::string, vector<double>> Realizations;
typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;


int main() {
    // TODO : Change Model->UpdateModel because it ain't parameters
    
    ///////////////////////
    //// Initialization ///
    ///////////////////////
    
    TestAssert::Init(true);
    std::string ModelType = "Multivariate";
    bool RealData = true;
    Data D;
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
        unsigned int NumberDimensions = 4, NbIndependentComponents = 2;
        shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
        shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimensions, BaseManifold);
        Model = make_shared<LongitudinalModel>(NbIndependentComponents, Manifold);
    }
    else if (ModelType == "Network")
    {
        unsigned int NbIndependentComponents = 5;
        
        std::string KernelMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/invKd_16.csv");
        auto KernelMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
        std::string InterpolationMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/Kxd_16.csv");
        auto InterpolationMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
        
        Model = make_shared<NetworkPropagationModel>(NbIndependentComponents, KernelMatrix, InterpolationMatrix);
    }
    else if(ModelType == "Network2")
    {
        unsigned int NbIndependentComponents = 5;
                
        std::string KernelMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/invKd_16.csv");
        auto KernelMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(KernelMatrixPath));
        std::string InterpolationMatrixPath("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/Kxd_16.csv");
        auto InterpolationMatrix = std::make_shared<NetworkPropagationModel2::MatrixType>(ReadData::OpenKernel(InterpolationMatrixPath));
        
        shared_ptr<AbstractManifold> Manifold = make_shared<ExponentialCurveManifold>(InterpolationMatrix->rows());
        //shared_ptr<AbstractManifold> Manifold = make_shared<LinearManifold>(InterpolationMatrix->rows());
        
        Model = make_shared<NetworkPropagationModel2>(NbIndependentComponents, Manifold, KernelMatrix, InterpolationMatrix);
    }
    
    
    /////////////////////////////
    /// Read or Simulate Data ///
    /////////////////////////////
    
    if(RealData)
    {
       D = ReadData::OpenFilesMultivariate();
    }
    else
    {
        Model->InitializeFakeRandomVariables();
        D = Model->SimulateData(300, 4, 6);
    }
    
    
    //////////////////////////
    /// Algorithm pipeline ///
    //////////////////////////
    auto Algo = make_shared<Algorithm>();
    Algo->SetModel(Model);
    Algo->SetSampler(Sampler);
    Algo->ComputeMCMCSAEM(D);
    
    
    
       
    ///////////////////
    /// Python call ///
    ///////////////////
    std::string filename = "/Users/igor.koval/PycharmProjects/RiemAlzh/plotGraph.py";
    std::string command = "python ";
    command += filename;
    //system(command.c_str());
    
    
    
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
