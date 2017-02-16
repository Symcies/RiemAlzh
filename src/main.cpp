#include <iostream>
#include <fstream>
#include <memory>

typedef double ScalarType;

#include "ModelSettings.h"
#include "DataSettings.h"
#include "AlgorithmSettings.h"

#include "Algorithm.h"

#include "MeshworkModel.h"
#include "FastNetworkModel.h"

#include "BlockedGibbsSampler.h"


#include "omp.h"
#include "tinyxml2.h"


using namespace std;

typedef vector< vector< pair< LinearAlgebra<ScalarType>::VectorType, double> > > Data;

// TODO : Change Model->UpdateModel because it ain't parameters
// TODO : Finish the input classes : Model, Algo and data
// TODO : Create an output file if it does not exist
// TODO : add a class member to encode the number of samplings to do per MCMC-sampling method
// TO ADD : Name of the output file
// TO ADD : If writing a file 1 out of n iteration - or only at the end 




//

int main(int argc, char* argv[]) {
    
    if(argc != 4) 
    {
        std::cout << "Usage with real data: " << " /path/to/executable " << " model_settings.xml " << " algorithm_settings " << "data_settings.xml" << std::endl;
        return 1;
    }

    /// GOOOD PART
        /// Initialize tests
    TestAssert::Init(false);
    
    /// Load the XML file arguments
    ModelSettings     MS(argv[1]);
    AlgorithmSettings AS(argv[2]);
    DataSettings      DS(argv[3]);
    
    /// Initialize the sampler
    std::shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    
    /// Initialize the model
    std::shared_ptr<AbstractModel> Model;
    //if(MS.GetType() == "Meshwork")  Model = make_shared<MeshworkModel>(MS);
    if(MS.GetType() == "FastNetwork") Model = make_shared<FastNetworkModel>(MS);
    
    /// Initialize the data
    Data D;
    if(DS.RealData())
    {
        int NbMaxOfSubjects = 250;
        D = ReadData::OpenFilesMultivariate(DS, NbMaxOfSubjects);
    }
    else 
    {
        Model->InitializeFakeRandomVariables();
        D = Model->SimulateData(DS);
    }
    /*
    Data D;
    bool ReadData = true;
    if(ReadData)
    {
        int NbMaxOfSubjects = 250;
        std::string FilePath = "/Users/igor.koval/Documents/Work/ok3/RiemAlzh/data/NormalizedThickness/MCIconvertAD/";
        D = ReadData::OpenFilesMultivariate(FilePath, NbMaxOfSubjects);
    }
    else
    {
        Model->InitializeFakeRandomVariables();
        D = Model->SimulateData(150, 3, 5);
    }
    */
    
    
    
    //////////////////////////
    /// Algorithm pipeline ///
    //////////////////////////
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

    #pragma omp parallel for
    for(int i = 0; i < 10; ++i)
    {
        
        //printf("(%d - %d)", i,omp_get_num_threads());
    }

}
