#include "ModelSettings.h"
#include "AlgorithmSettings.h"
#include "RealDataSettings.h"
#include "SamplerSettings.h"

#include "Observations.h"
#include "ReadData.h"

#include "ModelFactory.h"

#include "Algorithm.h"

#include "fit.h"

void fit(int argc, char* argv[]) {
  
  /// Load the XML file arguments
  io::ModelSettings     model_settings(argv[2]);
  io::AlgorithmSettings  algo_settings(argv[3]);
  io::RealDataSettings   data_settings(argv[4]);
  io::SamplerSettings sampler_settings(argv[5]);
  
  /// Observations
  Observations obs = io::ReadData::ReadObservations(data_settings);
  obs.InitializeGlobalAttributes();
  
  /// Model
  auto model = ModelFactory::NewModel(model_settings);
  auto initial_random_variables = model_settings.GetInitialRandomVariables();
  model->SetRandomVariableParameters(initial_random_variables);
  
  /// Algorithm pipeline
  Algorithm algo(algo_settings);
  algo.SetModel(model);
  algo.AddSamplers(sampler_settings);
  algo.ComputeMCMCSAEM(obs);
  
  /// Eventually simulate data (on option)
  /// Cannot be used ONLy beause the data_settings is reading the parameter <data-type>
  /// Maybe there should be two different
  /// TODO: yes I believe there should
  //io::SimulatedDataSettings simulated_data_settings(argv[4]);
  //model->SimulateData(simulated_data_settings);
}