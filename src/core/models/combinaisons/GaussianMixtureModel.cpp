#include "GaussianMixtureModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<AbstractRandomVariable> GaussianMixtureModel::GetRandomVariable(int key) const {
  int class_number = key_to_class_.at(key);
  
  return models_[class_number].GetRandomVariable(key);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianMixtureModel::GaussianMixtureModel(io::ModelSettings &model_settings) {
  /// Initialize the number of classes
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  number_of_classes_ = model_settings.GetNumberOfClasses();
  class_probabilities.set_size(number_of_classes_);
  models_.clear();
  
  for(size_t i = 0; i < number_of_classes_; ++i) {
    models_.push_back( GaussianModel(model_settings) );
    class_probabilities(i) = 1.0/number_of_classes_;
  }
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  
  
}

GaussianMixtureModel::~GaussianMixtureModel() {
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void GaussianMixtureModel::Initialize(const Observations &obs) {
  
  /// Data-related attributes
  manifold_dim_        = obs.GetSubjectObservations(0).GetCognitiveScore(0).size();
  subjects_tot_num_    = obs.GetNumberOfSubjects();
  obs_tot_num_         = obs.GetTotalNumberOfObservations();
  sum_obs_             = obs.GetTotalSumOfCognitiveScores();
  
  last_loglikelihood_.set_size(number_of_classes_, subjects_tot_num_);
  
  for(size_t i = 0; i < number_of_classes_; ++i) {
    InitialRVParameters second_rv_params = rv_params_;
    
    ScalarType mean     = 0;
    ScalarType var      = 0.001;
    ScalarType prop_var = 0.01;
    std::vector<ScalarType> params = {mean, var};
    second_rv_params["Gaussian"] = std::make_pair(params, prop_var);
    
    models_[i].SetRandomVariableParameters(second_rv_params);
    models_[i].Initialize(obs);
  }
  
  proposition_distribution_variance_["Gaussian"] = 0.001;
}

void GaussianMixtureModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                              const io::ModelSettings &model_settings) {
  
  
  int key_count = 0;
  for(size_t i = 0; i < number_of_classes_; ++i) {
    
    InitialRVParameters class_rv_params = rv_params_;
    
    ScalarType mean     = 2*i;
    ScalarType var      = 0.2*i;
    ScalarType prop_var = 0.01;
    std::vector<ScalarType> params = {mean, var};
    class_rv_params["Gaussian"] = std::make_pair(params, prop_var);
    
    models_[i].SetRandomVariableParameters(class_rv_params);
    models_[i].GetRandomVariable().SetKeyCount(key_count);
    models_[i].InitializeValidationDataParameters(data_settings, model_settings);
    
    for(size_t j = key_count; j < models_[i].GetRandomVariable().GetKeyCount() ; ++j)
      key_to_class_[j] = i;
    
    key_count = models_[i].GetRandomVariable().GetKeyCount() + 1; // The 1 is not useful but it prevents to go from class to another!
    
  }
}


void GaussianMixtureModel::UpdateModel(const Realizations &reals,
                                       const MiniBlock &block_info,
                                       const std::vector<std::string, std::allocator<std::string>> names) {
  
  int class_number = GetClassNumber(block_info);
  
  if(class_number != -1) {
    
    Realizations R = GetClassRealizations(class_number, reals);
    models_[class_number].UpdateModel(R, block_info, names);
    
  } else {
    
    for(size_t i = 0; i < number_of_classes_; ++i) {
      Realizations R = GetClassRealizations(i, reals);
      models_[i].UpdateModel(R, block_info, names);
    }
    
  }

}

AbstractModel::SufficientStatisticsVector GaussianMixtureModel::GetSufficientStatistics(const Realizations &reals,
                                                                                        const Observations &obs) {
  SufficientStatisticsVector general_suff_stat;
  
  suff_stat_sizes_.clear();
  for(size_t i = 0; i < number_of_classes_; ++i) {
    Realizations R = GetClassRealizations(i, reals);
    SufficientStatisticsVector s = models_[i].GetSufficientStatistics(R, obs);
    suff_stat_sizes_.push_back(s.size());
    general_suff_stat.insert(general_suff_stat.end(), s.begin(), s.end());
  }
  
  return general_suff_stat;
}

void GaussianMixtureModel::UpdateRandomVariables(const SufficientStatisticsVector &stoch_sufficient_stats) {
  
  unsigned int starting_element = 0;
  
  for(size_t i = 0; i < number_of_classes_; ++i) {
    /// Get the number of sufficient statistics involved in the i-th model
    unsigned int suff_stat_number_in_class = suff_stat_sizes_[i];
    
    /// Retrieve the suff stats of the i-th model
    SufficientStatisticsVector::const_iterator first_element = stoch_sufficient_stats.begin() + starting_element;
    starting_element += suff_stat_number_in_class;
    SufficientStatisticsVector::const_iterator last_element  = stoch_sufficient_stats.begin() +  starting_element;
    
    SufficientStatisticsVector class_ssv(first_element, last_element);
        
    /// Update the random variables of the i-th model
    models_[i].UpdateRandomVariables(class_ssv);
  }
}

Observations GaussianMixtureModel::SimulateData(io::SimulatedDataSettings &data_settings) {
  
  subjects_tot_num_ = data_settings.GetNumberOfSimulatedSubjects();
  
  Observations all_obs;
  unsigned int number_subjects_per_class = data_settings.GetNumberOfSimulatedSubjects() / number_of_classes_;
  data_settings.SetNumberOfSimulatedSubjects(number_subjects_per_class);
  
  for(size_t i = 0; i < number_of_classes_; ++i) {
    Observations class_obs = models_[i].SimulateData(data_settings);
    all_obs.AddObservations(class_obs);
  }
  
  /// Initialize the observation and model attributes
  all_obs.InitializeGlobalAttributes();
  
  return all_obs;
}

std::vector<AbstractModel::MiniBlock> GaussianMixtureModel::GetSamplerBlocks() const {
  
  std::vector<MiniBlock> blocks;
  
  for(size_t i = 0; i < number_of_classes_; ++i) {
     
    /// Individual blocks
    for(size_t j = 0; j < subjects_tot_num_; ++j) {
      MiniBlock indiv_block;
      for(size_t k = 0; k < manifold_dim_; ++k) {
        indiv_block.push_back(std::tuple<int, std::string, int>(0, "Gaussian#" + std::to_string(k) + "#" + std::to_string(i), j));
      }
      blocks.push_back(indiv_block);
    }
  }
  
  
  return blocks;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void GaussianMixtureModel::InitializeLogLikelihood(const Observations &obs) {
  
}

AbstractModel::VectorType GaussianMixtureModel::ComputeLogLikelihood(const Observations &obs,
                                                                     const MiniBlock &block_info) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  int class_number = GetClassNumber(block_info);
  
  return models_[class_number].ComputeLogLikelihood(obs, block_info);
}

ScalarType GaussianMixtureModel::ComputeIndividualLogLikelihood(const IndividualObservations &obs,
                                                                const int subject_num) {
  std::cerr << "GaussianMixtureModel > ComputeIndividualLogLikelihood - eThis function should not be computed for the Mixture Model --> Maybe overwritten" << std::endl; 
}

ScalarType GaussianMixtureModel::GetPreviousLogLikelihood(const MiniBlock &block_info) { 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  int class_number = GetClassNumber(block_info);
  return models_[class_number].GetPreviousLogLikelihood(block_info);
}

void GaussianMixtureModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
                                                    const MiniBlock &block_info) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  int class_number = GetClassNumber(block_info);
  models_[class_number].SetPreviousLogLikelihood(log_likelihood, block_info);
}




////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////


void GaussianMixtureModel::DisplayOutputs(const Realizations &reals) {
  
}

void GaussianMixtureModel::SaveCurrentState(unsigned int iteration_number, const Realizations &reals) {
  
}

void GaussianMixtureModel::SaveFinalState(const Realizations &reals, const Observations& obs) {}

void GaussianMixtureModel::SavePopulationFile() {}

void GaussianMixtureModel::SaveIndividualsFile(const Realizations &reals, const Observations& obs) {}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Overwriten method(s):
////////////////////////////////////////////////////////////////////////////////////////////////////

Realizations GaussianMixtureModel::SimulateRealizations() {
  
  Realizations reals;
  
  for(size_t i = 0; i < number_of_classes_; ++i) {
    Realizations class_realizations = models_[i].SimulateRealizations();
    
    for(auto it = class_realizations.begin(); it != class_realizations.end(); ++it) {
      
      int realization_key = it->first;
      std::string name = class_realizations.ReverseKeyToName(realization_key) + "#" + std::to_string(i);
      VectorType vector_of_reals = it->second;
      
      reals.AddRealizations(name, realization_key, vector_of_reals);
    }
    
  }
 
 
  return reals;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

int GaussianMixtureModel::GetClassNumber(const MiniBlock &block_info) {

  // TODO TODO TODO TODO TODO TODO TODO : return UNSIGNED INT or INT ?
  
  /// Check if it is the same class in the block_info ; if not --> Think about it
  int class_number = std::get<0>(block_info[0]);
  
  for(auto it = block_info.begin(); it != block_info.end(); ++it) {
    
    if(std::get<0>(*it) == class_number) 
      continue;
    else
      std::cerr << " Whoua ! Either allow only the same class within each miniblock; otherwise think how to do it. GAussianMixtureModel > UpdateModel";
      
  }
  
  if(class_number == -1) {
    int a = 0;
  }
  
  return class_number;
}


Realizations GaussianMixtureModel::GetClassRealizations(const unsigned int class_number, const Realizations &reals) {
  Realizations R;
  
  for(auto it = reals.begin(); it != reals.end(); ++it) {
    int realization_key = it->first;
    std::string real_name = reals.ReverseKeyToName(realization_key);
    
    std::string realization_name = real_name.substr(0, real_name.find_last_of("#"));
    int realization_number = stoi(real_name.substr(real_name.find_last_of("#") + 1)); 
    
    if(realization_number == class_number) {
      VectorType realizations = it->second;
      R.AddRealizations(realization_name, realization_key, realizations);
    }
  }
  
  return R;
}