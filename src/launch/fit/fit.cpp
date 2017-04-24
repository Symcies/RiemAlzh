#include "fit.h"

void fit(int argc, char* argv[]) {
  time_t start, init_comp, end_comp;
  start = time(0);

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
  init_comp = time(0);
  algo.ComputeMCMCSAEM(obs);
  end_comp = time(0);

  std::cout << "Initialisation duration: " << init_comp - start << std::endl;
  std::cout << "MCMCSAEM computations duration: " << end_comp - init_comp << std::endl;


  PythonUtils utils = PythonUtils(argv);
  utils.CallPythonScript(model_settings.GetOutputFileName());

  /// Eventually simulate data (on option)
  /// Cannot be used ONLy beause the data_settings is reading the parameter <data-type>
  /// Maybe there should be two different
  /// TODO: yes I believe there should
  //io::SimulatedDataSettings simulated_data_settings(argv[4]);
  //model->SimulateData(simulated_data_settings);
}