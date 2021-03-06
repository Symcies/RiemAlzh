#include <iostream>
#include <fstream>

typedef double ScalarType;
using namespace std;

#include "ModelSettings.h"
#include "DataSettings.h"
#include "AlgorithmSettings.h"

#include "Algorithm.h"

#include "MultivariateModel.h"
#include "NetworkModel.h"
#include "MeshworkModel.h"
#include "FastNetworkModel.h"

#include "BlockedGibbsSampler.h"

#include "Observations.h"

//#include "omp.h"
#include "tinyxml2.h"


int main(int argc, char* argv[]) {

  if(argc != 4)
  {
      std::cerr << "Usage with real data: " << " /path/to/executable " << " model_settings.xml " << " algorithm_settings " << "data_settings.xml" << std::endl;
  }

  /// Initialize tests
  TestAssert::Init(false);

  /// Load the XML file arguments
  io::ModelSettings     MS(argv[1]);
  io::AlgorithmSettings AS(argv[2]);
  io::DataSettings      DS(argv[3]);

  /// Initialize the sampler
  std::shared_ptr<AbstractSampler> Sampler = make_shared<BlockedGibbsSampler>();

  /// Initialize the model
  std::shared_ptr<AbstractModel> Model;
  if(MS.GetType() == "Meshwork")     Model = make_shared<MeshworkModel>(MS);
  if(MS.GetType() == "FastNetwork")  Model = make_shared<FastNetworkModel>(MS);
  if(MS.GetType() == "Network")      Model = make_shared<NetworkModel>(MS);
  if(MS.GetType() == "Multivariate") Model = make_shared<MultivariateModel>(MS);

  /// Initialize the data
  Observations Obs;
  if(DS.IsReal())
  {
    Obs = io::ReadData::ReadObservations(DS);
    Obs.InitializeGlobalAttributes();
  }
  else
  {
    Model->InitializeFakeRandomVariables();
    Obs = Model->SimulateData(DS);
  }

  /// Algorithm pipeline
  auto Algo = make_shared<Algorithm>(AS);
  Algo->SetModel(Model);
  Algo->SetSampler(Sampler);
  Algo->ComputeMCMCSAEM(Obs);


  return 0;

#pragma omp parallel for
  //for(int i = 0; i < 10; ++i)
    //printf("(%d - %d)", i,omp_get_num_threads());

  }


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
