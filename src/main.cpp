#include <iostream>
#include <fstream>

typedef double ScalarType;
using namespace std;

#include "ModelSettings.h"
#include "DataSettings.h"
#include "AlgorithmSettings.h"

#include "Algorithm.h"

#include "MultivariateModel.h"
#include "UnivariateModel.h"
//#include "NetworkModel.h"
//#include "MeshworkModel.h"
//#include "FastNetworkModel.h"

#include "BlockedGibbsSampler.h"

#include "Observations.h"

//#include "omp.h"
#include "tinyxml2.h"


int main(int argc, char* argv[]) {

  if(argc != 4)
  {
      std::cerr << "Usage with real data: " << " /path/to/executable " << " model_settings.xml " << " algorithm_settings.xml " << "data_settings.xml" << std::endl;
  }

  /// Load the XML file arguments
  io::ModelSettings     model_settings(argv[1]);
  io::AlgorithmSettings algo_settings(argv[2]);
  io::DataSettings      data_settings(argv[3]);

  /// Initialize the sampler
  std::shared_ptr<AbstractSampler> sampler = make_shared<BlockedGibbsSampler>();

  /// Initialize the model
  std::shared_ptr<AbstractModel> model;
  //if(model_settings.GetType() == "Meshwork")     model = make_shared<MeshworkModel>(model_settings);
  //if(model_settings.GetType() == "FastNetwork")  model = make_shared<FastNetworkModel>(model_settings);
  //if(model_settings.GetType() == "Network")      model = make_shared<NetworkModel>(model_settings);
  if(model_settings.GetType() == "Multivariate") model = make_shared<MultivariateModel>(model_settings);
  if(model_settings.GetType() == "Univariate")   model = make_shared<UnivariateModel>(model_settings);

  /// Initialize the data
  Observations obs;
  if(data_settings.IsReal())
  {
    obs = io::ReadData::ReadObservations(data_settings);
    obs.InitializeGlobalAttributes();
  }
  else
  {
    obs = model->SimulateData(data_settings);
  }

  /// Algorithm pipeline
  auto algo = make_shared<Algorithm>(algo_settings, model, sampler);
  algo->ComputeMCMCSAEM(obs);


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
