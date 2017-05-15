#include "MultivariateModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


MultivariateModel::MultivariateModel(io::ModelSettings &model_settings)
{
  /// Initialize the data dimension and the number of sources
  indep_sources_num_ = model_settings.GetIndependentSourcesNumber();
  output_file_name_ = model_settings.GetOutputFileName();
  std::remove((GV::BUILD_DIR + output_file_name_ ).c_str());
  std::remove((GV::BUILD_DIR + "LastRealizationOf" + output_file_name_ + "_pop").c_str());
  std::remove((GV::BUILD_DIR + "LastRealizationOf" + output_file_name_ + "_ind").c_str());
  acceptance_ratio_to_display_ = model_settings.GetAcceptanceRatioToDisplay();
}

MultivariateModel::~MultivariateModel()
{

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void MultivariateModel::Initialize(const Observations& obs)
{
  /// The function initialize different attributes of the model
  /// As well as the specific random variables and their realizations used by the model
  // TODO : Pitfall : How to have a smart initialization, that may come out of a .txt / .csv file
  // TODO : instead of using the same initialization or modifiying it in the code
  // TODO : Some cases may fall in between the two cases (default values needed)


  /// Data-related attributes
  manifold_dim_           = obs.GetSubjectObservations(0).GetCognitiveScore(0).size();
  subjects_tot_num_       = obs.GetNumberOfSubjects();
  individual_obs_date_    = obs.GetObservations();
  individual_time_points_ = obs.GetObservations();
  obs_tot_num_            = obs.GetTotalNumberOfObservations();
  sum_obs_                = obs.GetTotalSumOfCognitiveScores();


  /// Initialize the size of some parameters
  deltas_.set_size(manifold_dim_);
  block_.set_size(manifold_dim_);
  last_loglikelihood_.set_size(subjects_tot_num_);

  /// Population variables
  /// (Initialization of the random variable m_Random Variable
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  noise_ = std::make_shared<GaussianRandomVariable>(rv_params_.at("noise").first[0], rv_params_.at("noise").first[1]);

  rand_var_.AddRandomVariable("G", "Gaussian", rv_params_.at("G").first);
  asso_num_real_per_rand_var_["G"] = 1;
  proposition_distribution_variance_["G"] = rv_params_.at("G").second;


  for(size_t i = 1; i < manifold_dim_; ++i)
  {
    std::string name = "Delta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Delta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  proposition_distribution_variance_["Delta"] = rv_params_.at("Delta").second;

  for(int i = 0; i < indep_sources_num_*(manifold_dim_ - 1); ++i)
  {
    std::string name = "Beta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Beta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  proposition_distribution_variance_["Beta"] = rv_params_.at("Beta").second;

  /// Individual realizations
  /// (Initialization of the random variable m_Random Variable
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  rand_var_.AddRandomVariable("Ksi", "Gaussian", rv_params_.at("Ksi").first);
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  proposition_distribution_variance_["Ksi"] = rv_params_.at("Ksi").second;
  

  rand_var_.AddRandomVariable("Tau", "Gaussian", rv_params_.at("Tau").first);
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
  proposition_distribution_variance_["Tau"] = rv_params_.at("Tau").second;
  
  
  for(int i = 0; i < indep_sources_num_; ++i)
  {
    std::string name = "S#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("S").first);
    asso_num_real_per_rand_var_[name] = subjects_tot_num_;
  }
  proposition_distribution_variance_["S"] = rv_params_.at("S").second;
}

void MultivariateModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                           const io::ModelSettings &model_settings) {
  
  /// Initialize the model
  manifold_dim_ = data_settings.GetDimensionOfSimulatedObservations();
  indep_sources_num_ = model_settings.GetIndependentSourcesNumber();
  
  deltas_.set_size(manifold_dim_);
  block_.set_size(manifold_dim_);

  /// Population variables
  noise_ = std::make_shared<GaussianRandomVariable>(rv_params_.at("noise").first[0], rv_params_.at("noise").first[1]);

  rand_var_.AddRandomVariable("G", "Gaussian", rv_params_.at("G").first);
  asso_num_real_per_rand_var_["G"] = 1;

  for(size_t i = 1; i < manifold_dim_; ++i)
  {
    std::string name = "Delta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Delta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }

  for(int i = 0; i < indep_sources_num_*(manifold_dim_ - 1); ++i)
  {
    std::string name = "Beta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Beta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }

  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", rv_params_.at("Ksi").first);
  rand_var_.AddRandomVariable("Tau", "Gaussian", rv_params_.at("Tau").first);
  for(size_t i = 0; i < indep_sources_num_; ++i){
    rand_var_.AddRandomVariable("S#" + std::to_string(i), "Gaussian", rv_params_.at("S").first);
  }
}

void MultivariateModel::UpdateModel(const Realizations &reals, const MiniBlock& block_info,
            const std::vector<std::string, std::allocator<std::string>> names)
{
  /// Given a list of names (which in fact corresponds to the variables that have potentially changed),
  /// the function updates the parameters associated to these names

  int type = GetType(block_info);

  /// Possible parameters to update, depending on the name being in "vect<> names"
  bool compute_g       = false;
  bool compute_delta   = false;
  bool compute_basis   = false;
  bool compute_a       = false;
  bool compute_space_shift = false;
  bool compute_block   = false;
  bool individual_only = (type > -1);

  /// Parameters to update, depending on the names called
  //TODO: switch case?
  for(auto it = names.begin(); it != names.end(); ++it)
  {
    std::string name = it->substr(0, it->find_first_of("#"));
    if(name == "None")
    {
      continue;
    }
    else if(name == "Ksi" or name == "Tau")
    {
      continue;
    }
    else if("G" == name)
    {
      individual_only = false;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block = true;
    }
    else if("Delta" == name)
    {
      individual_only = false;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block = true;
    }
    else if("Beta" == name)
    {
      individual_only = false;
      compute_a = true;
      compute_space_shift = true;
    }
    else if("S" == name)
    {
      compute_space_shift = true;
    }
    else if("All" == name)
    {
      ComputeSubjectTimePoint(reals, -1);
      individual_only = false;
      compute_g = true;
      compute_delta = true;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block = true;
    }
    else
    {
      std::cerr << "The realization does not exist in the multivariate model > update model" << std::endl;
    }
  }


  // TODO : To parse it even faster, update just the coordinates within the names
  if(individual_only)     ComputeSubjectTimePoint(reals, type);
  if(compute_g)           g_ = exp(reals.at("G", 0));
  if(compute_delta)       ComputeDeltas(reals);
  if(compute_basis)       ComputeOrthonormalBasis();
  if(compute_a)           ComputeAMatrix(reals);
  if(compute_space_shift) ComputeSpaceShifts(reals);
  if(compute_block)       ComputeBlock(reals);
}

AbstractModel::SufficientStatisticsVector MultivariateModel::GetSufficientStatistics(const Realizations &reals, const Observations& obs)
{
  /// Computation of the suffisient statistics of the model
  /// Basically, it is a vector (named s1 to s9) of vectors

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

  /// s3 <- Ksi_i   &    s4 <- Ksi_i * Ksi_i
  VectorType s3 = reals.at("Ksi");
  VectorType s4 = reals.at("Ksi") % reals.at("Ksi");

  /// s5 <- Tau_i   &    s6 <- Tau_i * Tau_i
  VectorType s5 = reals.at("Tau");
  VectorType s6 = reals.at("Tau") % reals.at("Tau");

  /// s7 <- G
  VectorType s7(1, reals.at("G", 0));

  /// s8 <- beta_k
  VectorType s8((manifold_dim_-1) * indep_sources_num_);
  ScalarType * it_s8 = s8.memptr();
  for(size_t i = 0; i < s8.size(); ++i)
      it_s8[i] = reals.at("Beta#" + std::to_string(i), 0);

  /// s8 <- delta_k
  VectorType s9(manifold_dim_ - 1);
  ScalarType * it_s9 = s9.memptr();
  for(size_t i = 1; i < s9.size(); ++i)
      it_s9[i] = reals.at("Delta#" + std::to_string(i), 0);

  return {s1, s2, s3, s4, s5, s6, s7, s8, s9};

}


void MultivariateModel::UpdateRandomVariables(const SufficientStatisticsVector &sufficient_stats_vector)
{
  /// This function updates the random variables of the model rand_var_
  /// According to the sufficient statistic vector sufficient_stats_vector
  /// See documentation for further information about the update process

  /// Update the noise variance, sigma
  ScalarType noise_variance = sum_obs_;
  const ScalarType * it_s1 = sufficient_stats_vector[0].memptr();
  const ScalarType * it_s2 = sufficient_stats_vector[1].memptr();
  for(size_t i = 0; i < sufficient_stats_vector[0].size(); ++i)
      noise_variance += - 2 * it_s1[i] + it_s2[i];

  noise_variance /= obs_tot_num_ * manifold_dim_;
  noise_->SetVariance(noise_variance);

  /// Update Ksi : Mean & Variance
  ScalarType ksi_mean = 0.0, ksi_variance = 0.0;
  const ScalarType * it_s3 = sufficient_stats_vector[2].memptr();
  const ScalarType * it_s4 = sufficient_stats_vector[3].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    ksi_mean     += it_s3[i];
    ksi_variance += it_s4[i];
  }

  ksi_mean     /= subjects_tot_num_;
  ksi_variance -= subjects_tot_num_ * ksi_mean * ksi_mean;
  ksi_variance /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Ksi", {{"Mean", ksi_mean}, { "Variance", ksi_variance}});

  /// Update Tau : Mean & Variance
  ScalarType tau_mean = 0.0, tau_variance = 0.0;
  const ScalarType * it_s5 = sufficient_stats_vector[4].memptr();
  const ScalarType * it_s6 = sufficient_stats_vector[5].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
    tau_mean     += it_s5[i];
    tau_variance += it_s6[i];
  }

  tau_mean     /= subjects_tot_num_;
  tau_variance -= subjects_tot_num_ * tau_mean * tau_mean;
  tau_variance /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Tau", {{"Mean", tau_mean}, {"Variance", tau_variance}});

  /// Update G : Mean
  rand_var_.UpdateRandomVariable("G", {{"Mean", sufficient_stats_vector[6](0)}});

  /// Update Beta_k : Mean
  const ScalarType * it_s8 = sufficient_stats_vector[7].memptr();
  for(size_t i = 0; i < sufficient_stats_vector[7].size(); ++i)
    rand_var_.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", it_s8[i]}});

  /// Update Delta_k : Mean
  const ScalarType * it_s9 = sufficient_stats_vector[8].memptr();
  for(size_t i = 1; i < sufficient_stats_vector[8].size(); ++i)
    rand_var_.UpdateRandomVariable("Delta#" + std::to_string(i), {{"Mean", it_s9[i]}});

}

Observations MultivariateModel::SimulateData(io::SimulatedDataSettings &data_settings)
{
  /// This function simulates observations (Patients and their measurements y_ij at different time points t_ij)
  /// according to the model, with a given noise level e_ij, such that y_ij = f(t_ij) + e_ij
  /// Their is two option:
  /// 1) The model is not initialized (neither random variables of number of realizations) and it has to be
  /// This case corresponds to simulated data used to test the model, meaning if it can recover the random variables
  /// used to simulate the data
  /// 2) The model is already initialized, thus it should rely on its current state (random variables, realizations, ...)
  /// This case corresponds to new data, for a dataset augmentation for instance

  // TODO :
  /// PITFALL  : As for now, only the first option is implemented
  /// TODO: PITFALL2 : Take a closer look at / merge with InitializeFakeRandomVariables
  
  individual_obs_date_.clear();
  individual_time_points_.clear();
  
  subjects_tot_num_  = data_settings.GetNumberOfSimulatedSubjects();

  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;

  for(int i = 0; i < indep_sources_num_; ++i)
    asso_num_real_per_rand_var_["S#" + std::to_string(i)] = subjects_tot_num_;

  auto reals = SimulateRealizations();

  /// Update the model
  // NB: Cannot use UpdateModel because the subjectTimePoints cannot be updated yet
  g_ = exp(reals.at("G", 0));
  ComputeDeltas(reals);
  ComputeOrthonormalBasis();
  ComputeAMatrix(reals);
  ComputeSpaceShifts(reals);
  ComputeBlock(reals);

  /// Simulate the data
  std::random_device rand_device;
  std::mt19937 rand_num_gen( GV::TEST_RUN ? 1 : rand_device());

  std::uniform_int_distribution<int> uni(data_settings.GetMinimumNumberOfObservations(), data_settings.GetMaximumNumberOfObservations());
  UniformRandomVariable ran_time_points_num(60, 95);
  GaussianRandomVariable noise(0, noise_->GetVariance());

  /// Simulate the data
  Observations obs;
  for(int i = 0; i < subjects_tot_num_; ++i)
  {
    /// Get a random number of timepoints and sort them
    VectorType time_points = ran_time_points_num.Samples(uni(rand_num_gen));
    time_points.sort();
    individual_obs_date_.push_back(time_points);
    VectorType transformed_time_points = exp(reals.at("Ksi")(i)) * (time_points - reals.at("Tau")(i));
    individual_time_points_.push_back(transformed_time_points);

    /// Simulate the data base on the time-points
    IndividualObservations indiv_obs(time_points, i);
    std::vector<VectorType> cognitive_scores;
    for(size_t j = 0; j < time_points.size(); ++j)
    {
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

std::vector<AbstractModel::MiniBlock> MultivariateModel::GetSamplerBlocks() const
{

  int population_type = -1;
  std::vector<MiniBlock> blocks;

  /// Block G
  MiniBlock g_block;
  g_block.push_back(std::tuple<int, std::string, int>(population_type, "G", 0));
  blocks.push_back(g_block);

  /// Block Delta
  MiniBlock delta_block;
  for(size_t i = 1; i < manifold_dim_; ++i)
    delta_block.push_back(std::tuple<int, std::string, int>(population_type, "Delta#" + std::to_string(i), 0));
  blocks.push_back(delta_block);

  /// Block beta
  MiniBlock beta_block;
  for(size_t i = 0; i < indep_sources_num_*(manifold_dim_ - 1); ++i)
    beta_block.push_back(std::tuple<int, std::string, int>(population_type, "Beta#" + std::to_string(i), 0));
  blocks.push_back(beta_block);

  /// Individual variables --> Each block corresponds to one individual with all his/her associated random variables
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
    MiniBlock indiv_block;
    indiv_block.push_back(std::tuple<int, std::string, int>(i, "Ksi",i));
    indiv_block.push_back(std::tuple<int, std::string, int>(i, "Tau",i));
    for(size_t j = 0; j < indep_sources_num_; ++j)
      indiv_block.push_back(std::tuple<int, std::string, int>(i, "S#" + std::to_string(j), i));

    blocks.push_back(indiv_block);
  }

  return blocks;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void MultivariateModel::InitializeLogLikelihood(const Observations &obs) {
  last_loglikelihood_.set_size(subjects_tot_num_);
  ScalarType * l = last_loglikelihood_.memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i) 
    l[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);
}


AbstractModel::VectorType MultivariateModel::ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info)
{
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


ScalarType MultivariateModel::ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int subject_num)
{
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
  }

  log_likelihood /= - 2 * noise_->GetVariance();
  log_likelihood -= time_points * log(2 * noise_->GetVariance() * M_PI) / 2.0;

  return log_likelihood;
}

ScalarType MultivariateModel::GetPreviousLogLikelihood(const MiniBlock& block_info) {
  
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

void MultivariateModel::SetPreviousLogLikelihood(VectorType& log_likelihood, const MiniBlock& block_info) {
  
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

void MultivariateModel::DisplayOutputs(const Realizations &AR)
{
  /// It defines the outputs displayed on the terminal

  auto g = rand_var_.GetRandomVariable("G")->GetParameter("Mean");
  auto tau = rand_var_.GetRandomVariable("Tau");
  auto ksi = rand_var_.GetRandomVariable("Ksi");


  std::cout << "noise: "  << noise_->GetVariance();
  std::cout << " - G: "   << g;
  std::cout << " - T0: "  << tau->GetParameter("Mean") << " - Var(Tau): " << tau->GetParameter("Variance");
  std::cout << " - Ksi: " << ksi->GetParameter("Mean") << " - Var(Ksi): " << ksi->GetParameter("Variance") << std::endl;
}

void MultivariateModel::SaveCurrentState(unsigned int iter_num, const Realizations &reals)
{
  /// It saves the random variables / realizations / whatever model parameters
  /// Mainly needed for post processing
  
  std::ofstream log_file;
  
  if(!GV::TEST_RUN) {
    log_file.open(GV::BUILD_DIR + output_file_name_ , std::ofstream::out | std::ofstream::app);
    auto g = rand_var_.GetRandomVariable("G")->GetParameter("Mean");
    auto tau = rand_var_.GetRandomVariable("Tau");
    auto ksi = rand_var_.GetRandomVariable("Ksi");
    if (iter_num == 0){
      log_file << "iter Noise G TauMean TauVar KsiMean KsiVar" << std::endl;
    }
    // This part should be tuned by a xml file
    log_file << iter_num << " " << noise_->GetVariance() << " "
             << g << " " << tau->GetParameter("Mean") << " "
             << tau->GetParameter("Variance") << " "
             << ksi->GetParameter("Mean") << " "
             << ksi->GetParameter("Variance") << std::endl;
    log_file.close();
  }
  else {
    log_file.open(GV::TEST_DIR + "log_multivariate_file.txt", std::ofstream::out | std::ofstream::app);
    /// Save all the random variables parameters
    for(auto it = rand_var_.begin(); it != rand_var_.end(); ++it) {
      log_file << it->second->GetParameter(0) << " " << it->second->GetParameter(1) << " ";
    }
    
    /// Save all the realizations
    for(auto it = reals.begin(); it != reals.end(); ++it) {
      for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        log_file << *it2 << " ";
      }
    }
  }
}



void MultivariateModel::SavePopulationFile(){
  std::ofstream log_file_pop;

  log_file_pop.open(GV::BUILD_DIR + "LastRealizationOf" + output_file_name_ + "_pop", std::ofstream::out | std::ofstream::app);

  log_file_pop << "Noise " << noise_->GetVariance() << std::endl;
  auto block_g = rand_var_.GetRandomVariable("G")->GetParameter("Mean");
  log_file_pop << "G " << block_g << std::endl;
  auto tau = rand_var_.GetRandomVariable("Tau")->GetParameter("Mean");
  log_file_pop << "TauMean " << tau << std::endl;
  auto ksi = rand_var_.GetRandomVariable("Ksi")->GetParameter("Mean");
  log_file_pop << "KsiMean " << ksi << std::endl;
  log_file_pop << "Deltas ";
  for (int i = 0; i < deltas_.size(); i++){
    log_file_pop << deltas_[i] << " ";
  }
  log_file_pop << std::endl;

}

void MultivariateModel::SaveIndividualsFile(const Realizations &reals, const Observations &obs){
  std::ofstream log_file_ind;

  log_file_ind.open(GV::BUILD_DIR + "LastRealizationOf" + output_file_name_ + "_ind", std::ofstream::out | std::ofstream::app);

  // First part of file
  log_file_ind << "id Tau Ksi W" << std::endl;
  log_file_ind << "1 1 1 " << deltas_.size() << std::endl;

  // Second part of file
/*  // Label management
  log_file_ind << "id Tau Ksi ";
  for (int i = 0; i < deltas_.size(); i++) {
    log_file_ind << "W" << i << " ";
  }
  log_file_ind << std::endl;*/

  std::vector<std::string> val_realization = {"Tau", "Ksi"};

  int num_col = val_realization.size() + deltas_.size() + 1;
  int num_row = reals.at(val_realization[0]).size();
  std::cout << num_col << "  " << num_row;
  double realisations[num_row][num_col];

  // Id
  for (int j = 0; j < obs.GetNumberOfSubjects(); j++){
    realisations[j][0] = obs.GetId(j);
  }

  // Paramètres de realisation
  for (int i = 0; i < num_col - deltas_.size() - 1; i++) {
    std::string var = val_realization[i];
    for (int j = 0; j < reals.at(var).size(); j++) {
      realisations[j][1 + i] = reals.at(var)[j];
    }
  }

  // Space shifts
  for (int i = 0; i < deltas_.size(); i++) {
    for (int j = 0; j < space_shifts_.get_row(i).size(); j++) {
      realisations[j][1 + val_realization.size() + i] = space_shifts_(i, j);
    }
  }

  // Inversion: Going from row to columns
  for (int i = 0; i < num_row; i++) {
    for (int j = 0; j < num_col; j++) {
      log_file_ind << realisations[i][j] << " ";
    }
    log_file_ind << std::endl;
  }

  log_file_ind.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void MultivariateModel::ComputeSubjectTimePoint(const Realizations &reals, const int subject_num)
{
  /// The model introduces a time-warp for each subject, namely a time reparametrization
  /// For a given subject i, t_ij becomes alpha_i * ( t_ij - tau_i)
  /// alpha_i is the pace of disease propagation of patient i (Faster/Slower than average)
  /// tau_i is the disease onset of patient i (Beginning of the disease before/after the average time of conversion)

  /// If the subject number is -1, then the function recalculates the reparametrization for all the subjects
  if(subject_num != -1)
  {
      double acc_factor = exp(reals.at("Ksi", subject_num));
      double time_shift = reals.at("Tau", subject_num);
      individual_time_points_[subject_num] = acc_factor * (individual_obs_date_[subject_num] - time_shift);
  }
  else
  {

      for(size_t i = 0; i < subjects_tot_num_; ++i)
      {
          double acc_factor = exp(reals.at("Ksi")(i));
          double time_shift = reals.at("Tau")(i);

          individual_time_points_[i] = acc_factor * (individual_obs_date_[i] - time_shift);
      }
  }

}

void MultivariateModel::ComputeDeltas(const Realizations &reals)
{
  /// Create the vector of time shifts delta accross all coordinates
  /// delta = (delta_1, delta_2, ..., delta_N) where delta_1 = 0

  deltas_(0) = 0.0;
  ScalarType * d = deltas_.memptr();

  for(size_t i = 1; i < manifold_dim_; ++i) {
    d[i] = reals.at("Delta#" + std::to_string(i), 0);
  }
}

//TODO: change names, not clear at all
void MultivariateModel::ComputeOrthonormalBasis()
{
  /// It computes a basis of vector orthogonal to the geodesic gamma_derivative at t0, with respect to
  /// the scalar product defined on the riemannian manifold
  /// Further information about the mathematical operation in the documentation

  ScalarType v0 = - g_ / (g_ + 1) * exp(rand_var_.GetRandomVariable("Ksi")->GetParameter("Mean"));
  VectorType u_vector(manifold_dim_);
  ScalarType * u = u_vector.memptr();
  ScalarType * d = deltas_.memptr();

  for(size_t i = 0; i < manifold_dim_; ++i){
    u[i] = v0 * (1.0/g_ + exp(d[i]));
  }

  /// Compute the initial pivot vector u_vector
  double norm = u_vector.magnitude();
  u_vector(0) += copysign(1, -u_vector(0)) * norm;

  // Check the vectorise_row_wise function of the ArmadilloMatrixWrapper
  MatrixType u2(u_vector);
  double norm_u2 = u_vector.squared_magnitude();
  MatrixType final_matrix2 = (-2.0/norm_u2) * u2*u2.transpose();
  for(size_t i = 0; i < manifold_dim_; ++i){
    final_matrix2(i, i ) += 1;
  }

  orthog_basis_ = final_matrix2;
}

void MultivariateModel::ComputeAMatrix(const Realizations &reals)
{
  /// It computes the A matrix (ref documentation) based on the orthonormal basis B and the beta coefficients

  MatrixType new_a_matrix(manifold_dim_, indep_sources_num_);

  for(int i = 0; i < indep_sources_num_; ++i)
  {
    VectorType beta(manifold_dim_, 0.0);
    for(size_t j = 0; j < manifold_dim_ - 1; ++j)
    {
      std::string number = std::to_string(int(j + i*(manifold_dim_ - 1)));
      beta(j) = reals.at( "Beta#" + number, 0);
    }

    new_a_matrix.set_column(i, orthog_basis_ * beta);
  }

  a_matrix_ = new_a_matrix;
}

void MultivariateModel::ComputeSpaceShifts(const Realizations &reals)
{
  /// It computes the individual space shifts w_i (ref documentation or NIPS paper)
  MatrixType sufficient_stats_vector(indep_sources_num_, subjects_tot_num_);

  for(int i = 0; i < indep_sources_num_; ++i) {
    sufficient_stats_vector.set_row(i, reals.at("S#" + std::to_string(i)));
  }

  space_shifts_ = a_matrix_ * sufficient_stats_vector;
}

void MultivariateModel::ComputeBlock(const Realizations &reals)
{
  /// It compute a block used to optimise the number of computation within the likelihood function
  ScalarType * b = block_.memptr();
  ScalarType * d = deltas_.memptr();

  for(size_t i = 0; i < manifold_dim_; ++i)
  {
    ScalarType q = g_ * exp(-d[i]);
    b[i] = (q + 1) * (q + 1) / q;
  }
}


int MultivariateModel::GetType(const MiniBlock &block_info) {
  
  int type = std::get<0>(block_info[0]);
  
  for(auto it = block_info.begin(); it != block_info.end(); ++it) {
    int class_number = std::get<0>(*it);
    std::string real_name = std::get<1>(*it);
    int real_number = std::get<2>(*it);
    
    if(real_name == "G" || real_name == "Delta" || real_name == "Beta") {
      return -1;
    }
    if(type != real_number) {
      return -1;
    }
  }
  
  return type;
}


AbstractModel::VectorType MultivariateModel::ComputeParallelCurve(int subject_num, int obs_num)
{
  /// Given a subject i = subject_num and a timepoint t_ij where j = Observation number,
  /// it computes the f(t_i) value corresponding to the current model

  VectorType parallel_curve(manifold_dim_);
  ScalarType * p = parallel_curve.memptr();

  double time = individual_time_points_[subject_num](obs_num);
  ScalarType * d = deltas_.memptr();
  ScalarType * b = block_.memptr();
  ScalarType * w = space_shifts_.get_column(subject_num).memptr();

  for(size_t i = 0; i < manifold_dim_; ++i){
    p[i] = 1.0/(1.0 + g_ * exp(-d[i] - w[i]*b[i] - time));
  }

  return parallel_curve;
}

