#include "ModelSettings.h"
#include "AlgorithmSettings.h"
#include "SimulatedDataSettings.h"
#include "SamplerSettings.h"

#include "Observations.h"
#include "ReadData.h"

#include "ModelFactory.h"

#include "Algorithm.h"

#include "validate.h"

void validate(int argc, char* argv[]) {

  /// Load the XML file arguments
  io::ModelSettings        model_settings(argv[2]);
  io::AlgorithmSettings     algo_settings(argv[3]);
  io::SimulatedDataSettings data_settings(argv[4]);
  io::SamplerSettings    sampler_settings(argv[5]);
  
  /// Initialize the model
  auto model = ModelFactory::NewModel(model_settings);
  auto initial_random_variables = model_settings.GetInitialRandomVariables();
  model->SetRandomVariableParameters(initial_random_variables);
  
  /// Simulate the observations
  model->InitializeValidationDataParameters(data_settings, model_settings);
  Observations obs = model->SimulateData(data_settings);
  obs.InitializeGlobalAttributes();
  
  
  /// Initialize the model with other parameters
  auto second_random_variables = model_settings.GetSecondRandomVariables();
  model->SetRandomVariableParameters(second_random_variables);
  
    
  /// Algorithm pipeline
  Algorithm algo(algo_settings);
  algo.SetModel(model);
  algo.AddSamplers(sampler_settings);
  algo.ComputeMCMCSAEM(obs);

}