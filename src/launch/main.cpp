#include <iostream>
#include <fstream>

using namespace std;

#include "validate.h"
#include "fit.h"


int main(int argc, char* argv[]) {

  if(argc != 6)
  {
      std::cerr << "Usage with real data: " << "fit OR predict " << " /path/to/executable " << " model_settings.xml " << " algorithm_settings.xml " << "data_settings.xml" << std::endl;
  }

  
  if(string(argv[1]) == "fit")             fit(argc, argv);
  else if (string(argv[1]) == "validate")  validate(argc, argv);
  else std::cerr << "Second argument should be 'fit' or 'predict'" << std::endl;
  

  return 0;
/*
 /// Load the XML file arguments
  io::ModelSettings     model_settings(argv[2]);
  io::AlgorithmSettings algo_settings(argv[3]);
  std::shared_ptr<io::DataSettings> data_settings = Builder::BuilderDataSettings(argv[4]);

  /// Initialize the model
  std::shared_ptr<AbstractModel> model;
  //if(model_settings.GetType() == "Meshwork")     model = make_shared<MeshworkModel>(model_settings);
  //if(model_settings.GetType() == "FastNetwork")  model = make_shared<FastNetworkModel>(model_settings);
  //if(model_settings.GetType() == "Network")      model = make_shared<NetworkModel>(model_settings);
  if(model_settings.GetType() == "Multivariate") model = make_shared<MultivariateModel>(model_settings);
  if(model_settings.GetType() == "Univariate")   model = make_shared<UnivariateModel>(model_settings);
  
  /// Initialize the data
  Observations obs;
  if (data_settings->IsReal()) {
    std::shared_ptr<io::RealDataSettings> ds;
    ds = std::dynamic_pointer_cast<io::RealDataSettings>(data_settings);
    obs = io::ReadData::ReadObservations(*ds);
    obs.InitializeGlobalAttributes();
  } else {
    std::shared_ptr<io::SimulatedDataSettings> ds;
    ds = std::dynamic_pointer_cast<io::SimulatedDataSettings>(data_settings);
    obs = model->SimulateData(*ds, true);
  }
  
    /// Algorithm pipeline
  auto algo = make_shared<Algorithm>(algo_settings);
  
  algo.SetModel(model);
  
  /// Initialize the sampler
  std::shared_ptr<AbstractSampler> sampler = make_shared<BlockedGibbsSampler>();
  std::shared_ptr<AbstractSampler> sampler2 = make_shared<BlockedGibbsSampler>();
  
  
  
  algo.AddSampler(sampler, 25);
  algo.AddSampler(sampler2, 25);
 */
  
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
