#include "GaussianModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

GaussianModel::GaussianModel(io::ModelSettings &model_settings) {
  
}

GaussianModel::GaussianModel(const GaussianModel &) {
  
}

GaussianModel::~GaussianModel() {
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void GaussianModel::Initialize(const Observations &obs) {
  
  /// Data-related attributes
  manifold_dim_        = obs.GetSubjectObservations(0).GetCognitiveScore(0).size();
  subjects_tot_num_    = obs.GetNumberOfSubjects();
  individual_obs_date_ = obs.GetObservations();
  obs_tot_num_         = obs.GetTotalNumberOfObservations();
  sum_obs_             = obs.GetTotalSumOfCognitiveScores();
  
  last_loglikelihood_.set_size(subjects_tot_num_);
  
  /// Population variables
  noise_ = std::make_shared<GaussianRandomVariable>(rv_params_.at("noise").first[0], rv_params_.at("noise").first[1]);
  
  for(size_t i = 0; i < manifold_dim_; ++i) {
    std::string name = "Gaussian#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Gaussian").first);
    asso_num_real_per_rand_var_[name] = subjects_tot_num_;
  }

  proposition_distribution_variance_["Gaussian"] = rv_params_.at("Gaussian").second;
}

void GaussianModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                       const io::ModelSettings &model_settings) {
  
  /// Initialize the model
  manifold_dim_ = data_settings.GetDimensionOfSimulatedObservations();
  
  /// Population variables
  noise_ = std::make_shared<GaussianRandomVariable>(rv_params_.at("noise").first[0], rv_params_.at("noise").first[1]);
  
  for(size_t i = 0; i < manifold_dim_; ++i) {
    std::string name = "Gaussian#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Gaussian").first);
    asso_num_real_per_rand_var_[name] = subjects_tot_num_;
  }
}

void GaussianModel::UpdateModel(const Realizations &reals,
                                const MiniBlock &block_info,
                                const std::vector<std::string, std::allocator<std::string>> names) {
  
  int type = GetType(block_info);
  
  ComputeGaussianRealizations(reals, type);
  
}

AbstractModel::SufficientStatisticsVector GaussianModel::GetSufficientStatistics(const Realizations &reals,
                                                                                 const Observations &obs) {
  
  
  /// s1 <- y_ij * eta_ij    &    s2 <- eta_ij * eta_ij
  VectorType s1(obs_tot_num_), s2(obs_tot_num_);
  auto it_s1 = s1.begin(), it_s2 = s2.begin();
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      for(size_t j = 0; j < obs.GetNumberOfTimePoints(i); ++j)
      {
          VectorType parallel_curve = ComputeParallelCurve(i, j);
          *it_s1 = dot_product(parallel_curve, obs.GetSubjectCognitiveScore(i, j));
          *it_s2 = parallel_curve.squared_magnitude();
          ++it_s1, ++it_s2;
      }
  }
  
  SufficientStatisticsVector ssv = {s1, s2};
  
  /// s3 <- mu_i    &    s4 <- mu_i * mu_i
  for(size_t i = 0; i < manifold_dim_; ++i) {
    std::string name = "Gaussian#" + std::to_string(i);
    
    VectorType s3 = reals.at(name);
    VectorType s4 = reals.at(name) % reals.at(name);    
    
    ssv.push_back(s3);
    ssv.push_back(s4);
  }

  
  return ssv;
}

void GaussianModel::UpdateRandomVariables(const SufficientStatisticsVector &stoch_sufficient_stats) {
  /// Update the noise variance, sigma
  ScalarType noise_variance = sum_obs_;
  const ScalarType * it_s1 = stoch_sufficient_stats[0].memptr();
  const ScalarType * it_s2 = stoch_sufficient_stats[1].memptr();
  for(size_t i = 0; i < stoch_sufficient_stats[0].size(); ++i)
      noise_variance += - 2 * it_s1[i] + it_s2[i];

  noise_variance /= obs_tot_num_ * manifold_dim_;
  noise_->SetVariance(noise_variance);
  
  /// Update gaussian
  for(size_t i = 0; i < manifold_dim_; ++i) {
    int mean_number = (i+1)*2;
    int var_number = (i+1)*2 + 1;
    
    ScalarType gaussian_mean = 0.0, gaussian_variance = 0.0;
    
    const ScalarType * it_s3 = stoch_sufficient_stats[mean_number].memptr();
    const ScalarType * it_s4 = stoch_sufficient_stats[var_number].memptr();
    
    for(size_t j = 0; j < subjects_tot_num_; ++j) {
      gaussian_mean     += it_s3[j];
      gaussian_variance += it_s4[j];
    }
    
    gaussian_mean     /= subjects_tot_num_;
    gaussian_variance -= subjects_tot_num_ * gaussian_mean * gaussian_mean;
    gaussian_variance /= subjects_tot_num_;
    
    rand_var_.UpdateRandomVariable("Gaussian#" + std::to_string(i), {{"Mean", gaussian_mean}, {"Variance", gaussian_variance}});
  }

}

Observations GaussianModel::SimulateData(io::SimulatedDataSettings &data_settings) {
  individual_obs_date_.clear();
  
  subjects_tot_num_ = data_settings.GetNumberOfSimulatedSubjects();
  for(size_t i = 0; i < manifold_dim_; ++i) {
   asso_num_real_per_rand_var_["Gaussian#" + std::to_string(i)] = subjects_tot_num_; 
  }
  
  GaussianRealizations.resize(manifold_dim_);
  
  auto reals = SimulateRealizations();
  ComputeGaussianRealizations(reals, -1);
  
  ScalarType mean = GaussianRealizations[0].sum() / subjects_tot_num_;
  ScalarType var = 0;
  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    var += (GaussianRealizations[0](i) - mean) * (GaussianRealizations[0](i) - mean);
  }
  var /= subjects_tot_num_;
  std::cout << mean << " & " << var << std::endl;
  
   
  /// Simulate the data
  std::random_device rand_device;
  std::mt19937 rand_num_gen( GV::TEST_RUN ? 1 : rand_device());
  
  std::uniform_int_distribution<int> uni(data_settings.GetMinimumNumberOfObservations(), data_settings.GetMaximumNumberOfObservations());
  UniformRandomVariable ran_time_points_num(60, 95);
  GaussianRandomVariable noise(0, noise_->GetVariance());
  
  /// Simulate the data
  Observations obs;
  for(size_t i = 0; i < subjects_tot_num_; ++i) {
  
    VectorType time_points = ran_time_points_num.Samples(uni(rand_num_gen));
    time_points.sort();
    individual_obs_date_.push_back(time_points);
    
    IndividualObservations indiv_obs(time_points, i);
    std::vector<VectorType> cognitive_scores; 
    
    for(size_t j = 0; j < time_points.size(); ++j) {
      VectorType parallel_curve = ComputeParallelCurve(i, j);
      VectorType noise_sample = noise.Samples(manifold_dim_);
      cognitive_scores.push_back(parallel_curve + noise_sample);
    }
    
    indiv_obs.AddCognitiveScores(cognitive_scores);
    obs.AddIndividualData(indiv_obs);
  }
  
  /// Initialize the observation and model attributes
  obs.InitializeGlobalAttributes();

  return obs;
  
}

std::vector<AbstractModel::MiniBlock> GaussianModel::GetSamplerBlocks() const {
  std::vector<MiniBlock> blocks;
  
  /// Blocks
  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    MiniBlock indiv_block;
    for(size_t j = 0; j < manifold_dim_; ++j) {
      indiv_block.push_back(std::tuple<int, std::string, int>(0, "Gaussian#" + std::to_string(j), i));
    }
    blocks.push_back(indiv_block);
  }
  
  return blocks;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void GaussianModel::InitializeLogLikelihood(const Observations &obs) {
  last_loglikelihood_.set_size(subjects_tot_num_);
  ScalarType * l = last_loglikelihood_.memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i) 
    l[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);
}


AbstractModel::VectorType GaussianModel::ComputeLogLikelihood(const Observations &obs,
                                                              const MiniBlock &block_info) {
    /// It computes the likelihood of the model. For each subject i, it sums its likelihood, namely the distance,
  /// for each time t_ij, between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve
  
  int type = GetType(block_info);
  
  if(type == -1) {
    VectorType loglikelihood(subjects_tot_num_);
    ScalarType *l_ptr = loglikelihood.memptr();
    for (size_t i = 0; i < subjects_tot_num_; ++i)
      l_ptr[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);

    return loglikelihood;
  } else {
    return VectorType(1, ComputeIndividualLogLikelihood(obs.GetSubjectObservations(type), type));
  }
}

ScalarType GaussianModel::ComputeIndividualLogLikelihood(const IndividualObservations &obs,
                                                         const int subject_num) {
    /// Given a particular subject i, it computes its likelihood, namely the distance, for each observation t_ij,
  /// between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve

  ScalarType log_likelihood = 0;
  auto time_points = obs.GetNumberOfTimePoints();

  /// For each timepoints of the particular subject
  for(size_t i = 0; i < time_points; ++i)
  {
    auto& indiv_cog_scores = obs.GetCognitiveScore(i);
    VectorType parallel_curve = ComputeParallelCurve(subject_num, i);
    log_likelihood += (indiv_cog_scores - parallel_curve).squared_magnitude();
    //std::cout << indiv_cog_scores(0) << " - " << parallel_curve(0) << " = " << indiv_cog_scores(0) - parallel_curve(0) << std::endl;
  }

  log_likelihood /= - 2 * noise_->GetVariance();
  log_likelihood -= time_points * log(2 * noise_->GetVariance() * M_PI) / 2.0;

  return log_likelihood;
}

ScalarType  GaussianModel::GetPreviousLogLikelihood(const MiniBlock &block_info) {
  int type = GetType(block_info);
  
  if (type == -1) {
    return last_loglikelihood_.sum();
  }
  else if (type >= 0 && type <= last_loglikelihood_.size()) {
    return last_loglikelihood_(type);
  }
  else {
    std::cerr << "there is something wrong with the type";
  }
}

void GaussianModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
                                             const MiniBlock &block_info) {
  int type = GetType(block_info);
  
  if (type == -1) {
    last_loglikelihood_ = log_likelihood;
  }  else if (type >= 0 && type <= last_loglikelihood_.size()) {
    last_loglikelihood_(type) = log_likelihood.sum();
  }
  else {
    std::cerr << "there is something wrong with the type";
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////


void GaussianModel::DisplayOutputs(const Realizations &reals) {
  auto gaussian = rand_var_.GetRandomVariable("Gaussian#0");
  
  std::cout << "noise: " << noise_->GetVariance();
  std::cout << " - Gaussian mean: " << gaussian->GetParameter("Mean");
  std::cout << " - Gaussian variance: " << gaussian->GetParameter("Variance") << std::endl;
  
  
  ScalarType mean = GaussianRealizations[0].sum() / subjects_tot_num_;
  ScalarType var = 0;
  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    var += (GaussianRealizations[0](i) - mean) * (GaussianRealizations[0](i) - mean);
  }
  var /= subjects_tot_num_;
  //std::cout << mean << " & " << var << std::endl;
  
}



void GaussianModel::SaveCurrentState(unsigned int IterationNumber, const Realizations &reals) {
  // TODO TODO TODO TODO TODO TODO 
  // We'll see later on what we need
}

void GaussianModel::SaveFinalState(const Realizations &reals, const Observations& obs) {}

void GaussianModel::SavePopulationFile() {}

void GaussianModel::SaveIndividualsFile(const Realizations &reals, const Observations& obs) {}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void GaussianModel::ComputeGaussianRealizations(const Realizations &reals, const int type) {
  
  if(type == -1) {
    for(size_t i = 0; i < manifold_dim_; ++i) {
      GaussianRealizations[i] = reals.at("Gaussian#" + std::to_string(i));
    }
  } else {
    for(size_t i = 0; i < manifold_dim_; ++i) {
      GaussianRealizations[i](type) = reals.at("Gaussian#" + std::to_string(i))(type);
    }
  }
}



AbstractModel::VectorType GaussianModel::ComputeParallelCurve(int subjects_num, int obs_num) {
  
  VectorType parallel_curve(manifold_dim_);
  ScalarType * p = parallel_curve.memptr();
  ScalarType time = individual_obs_date_[subjects_num](obs_num);
   
  
  for(size_t i = 0; i < manifold_dim_; ++i) {
    ScalarType individual_variable = GaussianRealizations[i](subjects_num);
    p[i] = individual_variable * time;
  }
  
  return parallel_curve;
}


int GaussianModel::GetType(const MiniBlock &block_info) {
  int type = std::get<2>(block_info[0]);
  
  for(auto it = block_info.begin(); it != block_info.end(); ++it) {
    int class_number = std::get<0>(*it);
    std::string real_name = std::get<1>(*it);
    real_name = real_name.substr(0, real_name.find_first_of("#"));
    int real_number = std::get<2>(*it);
    
    if(real_name != "Gaussian" && real_name != "All") {
      std::cerr << "What is going on in Gaussian Model? > GetType method with realization : " << real_name <<  std::endl;
    }
    if(real_name == "All") 
      return -1;
    if(type != real_number) {
      return -1;
    }
  }
  
  return type;
}
