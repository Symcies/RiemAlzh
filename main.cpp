#include <iostream>
#include <fstream>
#include <memory>


#include "Manifolds/PropagationManifold.h"
#include "Models/LongitudinalModel.h"
#include "Models/UnivariateModel.h"
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



Data 
OpenFiles()
{
    Data D;
    std::ifstream IndivID ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Data_Test_SAEM_group.txt");
    std::ifstream DataX ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Data_Test_SAEM_X.txt");
    std::ifstream DataY ("/Users/igor.koval/Documents/Git/RiemAlzh/datatest/Data_Test_SAEM_Y_3.txt");
    
    /// Open the DATA_text_SAEM_group file;
    if(IndivID.is_open())
    {
        int i = 0;
        string line;
        vector< pair< vector<double>, double> > IndivData;
        pair< vector<double>, double> Observations;
        
        while(getline(IndivID, line))
        {
            int j = std::stoi(line);
            if(j == i)
            {
                IndivData.push_back(Observations);
            }
            if(j != i)
            {
                D.push_back(IndivData);
                IndivData.clear();
                IndivData.push_back(Observations);
                i = j;
            }
        }
        D.push_back(IndivData);
    }
    else { cout << "Unable to open indiv id's"; }
    
    /// Open the DATA : timepoints
    if(DataX.is_open())
    {
        string line;
        getline(DataX, line);
        for(auto it = D.begin(); it != D.end(); it++)
        {
            for(auto it2 = it->begin(); it2 != it->end(); ++it2)
            {
                it2->second = stod(line);
                getline(DataX, line);
            }
        }
        std::cout << line << std::endl;
    }
    else { cout << "Unable to open timepoints"; }
    
    /// Open the DATA : observations
    if(DataY.is_open())
    {
        string line;
        getline(DataY, line);
        for(auto it = D.begin(); it != D.end(); ++it)
        {
            for(auto it2 = it->begin(); it2 != it->end(); ++it2)
            {
                std::vector<double> X;
                for(int i = 0; i < 4; ++i)
                {
                    X.push_back(stod(line));
                    getline(DataY, line);
                }
                it2->first = X;
            }
        }
        std::cout << line << std::endl;
    }
    else { cout << "Unable to open the observation values"; }
    
    if( D[0].size() == 0) 
    {
        D.erase(D.begin());
    }
 
    return D;
}

int main() {

    //// INITIALIZATION ///
    unsigned int NumberDimension = 4;
    unsigned int NumberIndependentComponents = 1;
    clock_t start = clock();
    
    /// TESTS ///
    bool Active = true;
    TestAssert::Init(Active);
    
    
    /// Open a text file /// 
    shared_ptr<Data> D = std::make_shared<Data>(OpenFiles());
    


    /// Base Manifold, Manifold, Model & Sampler ///
    
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractManifold> Manifold = make_shared<PropagationManifold>(NumberDimension, BaseManifold);
    shared_ptr<AbstractModel> Model = make_shared<LongitudinalModel>(NumberIndependentComponents, Manifold);
    //shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    
     
    
    /// Univariate Model
    /*
    shared_ptr<AbstractBaseManifold> BaseManifold = make_shared<LogisticBaseManifold>();
    shared_ptr<AbstractModel> Model = make_shared<UnivariateModel>(BaseManifold);
    //shared_ptr<AbstractSampler> Sampler = make_shared<HMWithinGibbsSampler>();
    shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();
    */
    
    
    /// DATA GENERATION ///
    /*
    Model->InitializeFakeRandomVariables();
    shared_ptr<Data> D = make_shared<Data>( Model->SimulateData(200, 6, 8) );
    */

    /// Model
    Model->Initialize();


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
