#include "FastNetworkModel.h"

FastNetworkModel::FastNetworkModel(io::ModelSettings &MS)
{
  indep_components_nb_ = MS.GetIndependentSourcesNumber();
  acceptance_ratio_to_display_ = MS.GetAcceptanceRatioToDisplay();
  
  std::string kernel_matrix_path = MS.GetInvertKernelPath();
  std::string interp_matrix_path = MS.GetInterpolationKernelPath();

  invert_kernel_matrix_ = io::ReadData::OpenKernel(kernel_matrix_path).transpose();
  interpolation_matrix_ = io::ReadData::OpenKernel(interp_matrix_path);

  control_points_nb_ = invert_kernel_matrix_.columns();
  manifold_dim_ = interpolation_matrix_.rows();
  nus_.set_size(manifold_dim_);
  deltas_.set_size(manifold_dim_);
  block1_.set_size(manifold_dim_);
  block2_.set_size(manifold_dim_);
}

FastNetworkModel::~FastNetworkModel()
{



}

void FastNetworkModel::Initialize(const Observations& obs)
{
  /// Data-related attributes
  subjects_tot_num_       = obs.GetNumberOfSubjects();
  individual_obs_date_    = obs.GetObservations();
  individual_time_points_ = obs.GetObservations();
  obs_tot_num_            = obs.GetTotalNumberOfObservations();
  sum_obs_                = obs.GetTotalSumOfLandmarks();

   /// Population variables
  noise_ = std::make_shared<GaussianRandomVariable>( 0.0, 0.000001 );

  rand_var_.AddRandomVariable("P", "Gaussian", rv_params_.at("P").first);
  asso_num_real_per_rand_var_["P"] = 1;
  proposition_distribution_variance_["P"] = rv_params_.at("P").second;

  for(int i = 1; i < control_points_nb_; ++i)
  {
    std::string name = "Delta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Delta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  proposition_distribution_variance_["Delta"] = rv_params_.at("Delta").second;


  rand_var_.AddRandomVariable("Nu", "Gaussian", rv_params_.at("Nu").first);
  asso_num_real_per_rand_var_["Nu"] = control_points_nb_;
  proposition_distribution_variance_["Nu"] = rv_params_.at("Nu").second;

  for(int i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
  {
    std::string name = "Beta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Beta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  proposition_distribution_variance_["Beta"] = rv_params_.at("Beta").second;


  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", rv_params_.at("Ksi").first);
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  proposition_distribution_variance_["Ksi"] = rv_params_.at("Ksi").second;

  rand_var_.AddRandomVariable("Tau", "Gaussian", rv_params_.at("Tau").first);
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
  proposition_distribution_variance_["Tau"] = rv_params_.at("Tau").second;
  
  
  for(int i = 0; i < indep_components_nb_; ++i)
  {
    std::string name = "S#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("S").first);
    asso_num_real_per_rand_var_[name] = subjects_tot_num_;
  }
  proposition_distribution_variance_["S"] = rv_params_.at("S").second;

}

void FastNetworkModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                          const io::ModelSettings &model_settings) {
  /// Initialize the model
  // TODO : Check which one has to be updated --> Some (a lot) are initialized in the constructor
  // DO AN ASSERT : manifold_dim_ = data_settings.GetDimensionOfSimulatedObservations();
  indep_components_nb_ = model_settings.GetIndependentSourcesNumber();
  last_loglikelihood_.set_size(subjects_tot_num_);
  
  /// Population variables
  rand_var_.AddRandomVariable("P", "Gaussian", rv_params_.at("P").first);
  asso_num_real_per_rand_var_["P"] = 1;

  for(int i = 1; i < control_points_nb_; ++i)
  {
    std::string name = "Delta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Delta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  
  rand_var_.AddRandomVariable("Nu", "Gaussian", rv_params_.at("Nu").first);
  asso_num_real_per_rand_var_["Nu"] = control_points_nb_;
  
  for(int i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
  {
    std::string name = "Beta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Beta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  

  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", rv_params_.at("Ksi").first);
  rand_var_.AddRandomVariable("Tau", "Gaussian", rv_params_.at("Tau").first); 
  
  for(int i = 0; i < indep_components_nb_; ++i)
  {
    std::string name = "S#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("S").first);
    
  }
   
  
}


void FastNetworkModel::UpdateModel(const Realizations &reals, const MiniBlock &block_info,
                                  const std::vector<std::string, std::allocator<std::string>> names)
{

  int type = GetType(block_info);

  bool compute_position = false;
  bool compute_delta = false;
  bool compute_nu = false;
  bool compute_basis = false;
  bool compute_a = false;
  bool compute_space_shift = false;
  bool compute_block1 = false;
  bool compute_block2 = false;

  bool individual_only = (type > -1);
  
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
    else if(name == "Nu")
    {
      individual_only = false;
      compute_nu = true;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block2 = true;
      continue;
    }
    else if(name == "Delta")
    {
      individual_only = false;
      compute_delta = true;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block1 = true;
      continue;
    }
    else if(name == "P")
    {
      individual_only = false;
      compute_position = true;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block1 = true;
      compute_block2 = true;
      continue;
    }
    else if(name == "Beta")
    {
      individual_only = false;
      compute_a = true;
      compute_space_shift = true;
      continue;
    }
    else if(name == "S")
    {
      compute_space_shift = true;
      continue;
    }
    else if(name == "All")
    {
      ComputeSubjectTimePoint(reals, -1);
      individual_only = false;
      compute_position = true;
      compute_delta = true;
      compute_nu = true;
      compute_basis = true;
      compute_a = true;
      compute_space_shift = true;
      compute_block1 = true;
      compute_block2 = true;
      break;
    }
    else
    {
      std::cerr << "PROBLEM WITH FAST NETWORK MODEL" << std::endl;
    }
  }

  // TODO : To parse it even faster, update just the coordinates within the names
  if(individual_only) ComputeSubjectTimePoint(reals, type);

  if(compute_position)    init_pos_ = exp(reals.at("P", 0));
  if(compute_delta)       ComputeDeltas(reals);
  if(compute_nu)          ComputeNus(reals);
  if(compute_basis)       ComputeOrthonormalBasis();
  if(compute_a)           ComputeAMatrix(reals);
  if(compute_space_shift) ComputeSpaceShifts(reals);
  if(compute_block1)      ComputeBlock1();
  if(compute_block2)      ComputeBlock2();

}


FastNetworkModel::SufficientStatisticsVector FastNetworkModel::GetSufficientStatistics(const Realizations& simulated_real,  const Observations& obs)
{

  /// s1 <- y_ij * eta_ij    &    s2 <- eta_ij * eta_ij
  VectorType s1(obs_tot_num_), s2(obs_tot_num_);
  auto it_s1 = s1.begin(), it_s2 = s2.begin();
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      for(size_t j = 0; j < obs.GetNumberOfTimePoints(i); ++j)
      {
          VectorType parallel_curve = ComputeParallelCurve(i, j);
          *it_s1 = dot_product(parallel_curve, obs.GetSubjectLandmark(i, j));
          *it_s2 = parallel_curve.squared_magnitude();
          ++it_s1, ++it_s2;
      }
  }

  /// s3 <- ksi_i * ksi_i
  VectorType s3 = simulated_real.at("Ksi") % simulated_real.at("Ksi");

  /// s4 <- Tau_i   &    s5 <- Tau_i * Tau_i
  VectorType s4 = simulated_real.at("Tau");
  VectorType s5 = simulated_real.at("Tau") % simulated_real.at("Tau");

  /// s6 <- p0
  VectorType s6(1, simulated_real.at("P", 0));

  /// s7 <- beta_k
  VectorType s7((manifold_dim_-1) * indep_components_nb_);
  ScalarType * it_s7 = s7.memptr();
  for(size_t i = 0; i < s7.size(); ++i)  {
      it_s7[i] = simulated_real.at("Beta#" + std::to_string(i), 0);
  }

  /// s8 <- delta_k
  VectorType s8(control_points_nb_ - 1);
  ScalarType * it_s8 = s8.memptr();
  for(size_t i = 1; i < s8.size(); ++i)  {
      it_s8[i] = simulated_real.at("Delta#" + std::to_string(i), 0);
  }

  /// s9 <- nu_k, s10 = nu_k * nu_k
  VectorType s9 = simulated_real.at("Nu");
  VectorType s10 = simulated_real.at("Nu") % simulated_real.at("Nu");


  SufficientStatisticsVector sufficient_stat_vec = {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10};
  return sufficient_stat_vec;
}

void FastNetworkModel::UpdateRandomVariables(const SufficientStatisticsVector &sufficient_stats_vector)
{
  /// Update the noise variance, sigma
  ScalarType noise_variance = sum_obs_;
  const ScalarType * it_s1 = sufficient_stats_vector[0].memptr();
  const ScalarType * it_s2 = sufficient_stats_vector[1].memptr();
  for(size_t i = 0; i < sufficient_stats_vector[0].size(); ++i)
      noise_variance += - 2 * it_s1[i] + it_s2[i];

  noise_variance /= obs_tot_num_ * manifold_dim_;
  noise_->SetVariance(noise_variance);

  /// Update ksi and sigma_ksi
  ScalarType ksi_variance = 0.0;
  const ScalarType * it_s3 = sufficient_stats_vector[2].memptr();
  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    ksi_variance += it_s3[i];
  }
  ksi_variance /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Ksi", {{"Variance", ksi_variance}});

  /// Update tau and sigma_tau
  ScalarType tau_mean = 0.0, tau_variance = 0.0;
  const ScalarType * it_s4 = sufficient_stats_vector[3].memptr();
  const ScalarType * it_s5 = sufficient_stats_vector[4].memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    tau_mean     += it_s4[i];
    tau_variance += it_s5[i];
  }

  tau_mean     /= subjects_tot_num_;
  tau_variance -= subjects_tot_num_ * tau_mean * tau_mean;
  tau_variance /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Tau", {{"Mean", tau_mean}, {"Variance", tau_variance}});

  /// Update p0
  rand_var_.UpdateRandomVariable("P", {{"Mean", sufficient_stats_vector[5](0)}});

  /// Update beta_k
  const ScalarType * it_s7 = sufficient_stats_vector[6].memptr();
  for(size_t i = 0; i < (manifold_dim_-1) * indep_components_nb_; ++i)
    rand_var_.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", it_s7[i]}});

  /// Update delta_k
  const ScalarType * it_s8 = sufficient_stats_vector[7].memptr();
  for(size_t i = 1; i < control_points_nb_ - 1; ++i)
    rand_var_.UpdateRandomVariable("Delta#" + std::to_string(i), {{"Mean", it_s8[i]}});
  

  /// Update nu_k = v0 and sigma_nu
  ScalarType v0 = 0.0, nu_variance = 0.0;
  const ScalarType * it_s9  = sufficient_stats_vector[8].memptr();
  const ScalarType * it_s10 = sufficient_stats_vector[9].memptr();
  
  for(size_t i = 0; i < sufficient_stats_vector[8].size(); ++i) {
    v0          +=  it_s9[i];
    nu_variance += it_s10[i];
  }
  
  v0          /= subjects_tot_num_;
  nu_variance -= subjects_tot_num_ * v0 * v0;
  nu_variance /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Nu", {{"Mean", v0}, {"Variance", nu_variance}});
}



Observations FastNetworkModel::SimulateData(io::SimulatedDataSettings& data_settings)
{

  individual_obs_date_.clear();
  individual_time_points_.clear();
  
  subjects_tot_num_ = data_settings.GetNumberOfSimulatedSubjects();
  
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;

  for(int i = 0; i < indep_components_nb_; ++i)
    asso_num_real_per_rand_var_["S#" + std::to_string(i)] = subjects_tot_num_;

  auto reals = SimulateRealizations();


  /// Update the model
  init_pos_ = exp(reals.at("P")(0));
  ComputeDeltas(reals);
  ComputeNus(reals);
  ComputeOrthonormalBasis();
  ComputeAMatrix(reals);
  ComputeSpaceShifts(reals);
  ComputeBlock1();
  ComputeBlock2();


  /// Simulate the data
  std::random_device rand_device;
  std::mt19937 rand_num_gen(GV::TEST_RUN ? 1 : rand_device());
  std::uniform_int_distribution<int> unif_distrib(data_settings.GetMinimumNumberOfObservations(), data_settings.GetMaximumNumberOfObservations());
  UniformRandomVariable time_points_num(60, 95);
  
  /// Simulate the data
  Observations obs;
  for(int i = 0; i < subjects_tot_num_; ++i)
  {
    /// Get a random number of timepoints and sort them
    VectorType time_points = time_points_num.Samples(unif_distrib(rand_num_gen));
    time_points.sort();
    individual_obs_date_.push_back(time_points);
    VectorType transformed_time_points = exp(reals.at("Ksi")(i)) * (time_points - reals.at("Tau")(i) );
    individual_time_points_.push_back(transformed_time_points);
    
    /// Simulate the data base on the time-points
    IndividualObservations indiv_obs(time_points, i);
    std::vector<VectorType> landmarks;
    for(size_t j = 0; j < time_points.size(); ++j)
    {
      landmarks.push_back(ComputeParallelCurve(i, j) + noise_->Samples(manifold_dim_));
    }

    indiv_obs.AddLandmarks(landmarks);
    obs.AddIndividualData(indiv_obs);
  }

  /// Initialize the observation and model attributes
  obs.InitializeGlobalAttributes();

  return obs;

}

std::vector<AbstractModel::MiniBlock> FastNetworkModel::GetSamplerBlocks() const
{
  int population_type = -1;
  std::vector<MiniBlock> blocks;

  
  /// Insert p0;
  MiniBlock block_p;
  block_p.push_back(std::tuple<int, std::string, int>(population_type, "P", 0));
  blocks.push_back(block_p);

  /// Insert Beta_k
  MiniBlock block_beta;
  for(size_t i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
    block_beta.push_back(std::tuple<int, std::string, int>(population_type, "Beta#" + std::to_string(i), 0));
  blocks.push_back(block_beta);


  /// Insert Delta_k
  MiniBlock block_delta;
  for(size_t i = 1; i < control_points_nb_; ++i)
      block_delta.push_back(std::tuple<int, std::string, int>(population_type, "Delta#" + std::to_string(i), 0));
  blocks.push_back(block_delta);

  /// Insert Nu_k
  MiniBlock block_nu;
  for(size_t i = 0; i < control_points_nb_; ++i)
      block_nu.push_back(std::tuple<int, std::string, int>(population_type, "Nu", i));
  blocks.push_back(block_nu);

  /// Individual variables

  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      MiniBlock indiv_block;
      indiv_block.push_back(std::tuple<int, std::string, int>(i, "Ksi", i));
      indiv_block.push_back(std::tuple<int, std::string, int>(i, "Tau", i));
      for(size_t j = 0; j < indep_components_nb_; ++j)
            indiv_block.push_back(std::tuple<int, std::string, int>(i, "S#" + std::to_string(j), i));

      blocks.push_back(indiv_block);
  }
  
  return blocks;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void FastNetworkModel::InitializeLogLikelihood(const Observations &obs) {
  last_loglikelihood_.set_size(subjects_tot_num_);
  ScalarType * l = last_loglikelihood_.memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i) 
    l[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);
}

AbstractModel::VectorType FastNetworkModel::ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info)
{
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

ScalarType FastNetworkModel::ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int subject_num)
{
  /// Get the data
  ScalarType log_likelihood = 0;
  auto time_points_num = obs.GetNumberOfTimePoints();

#pragma omp parallel for reduction(+:log_likelihood)
  for(size_t i = 0; i < time_points_num; ++i)
  {
      auto& it = obs.GetLandmark(i);
      VectorType parallel_curve2 = ComputeParallelCurve(subject_num, i);
      log_likelihood += (it - parallel_curve2).squared_magnitude();
  }

  log_likelihood /= -2*noise_->GetVariance();
  log_likelihood -= time_points_num * log(2 * noise_->GetVariance() * M_PI) / 2.0;

  return log_likelihood;
}

ScalarType FastNetworkModel::GetPreviousLogLikelihood(const MiniBlock &block_info) {

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


void FastNetworkModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
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

void FastNetworkModel::DisplayOutputs(const Realizations& simulated_real)
{
  auto p0 = rand_var_.GetRandomVariable("P")->GetParameter("Mean");
  auto tau = rand_var_.GetRandomVariable("Tau");

  auto nu = rand_var_.GetRandomVariable("Nu");

  double nu_max = simulated_real.at("Nu").max_value();
  double nu_min = simulated_real.at("Nu").min_value();


  double delta_min = simulated_real.at("Delta#1", 0);
  double delta_max = delta_min;
  for(size_t i = 1; i < 258; ++i)
  {
      double delta_k = simulated_real.at("Delta#" + std::to_string(i), 0);
      delta_max = std::max(delta_max, delta_k);
      delta_min = std::min(delta_min, delta_k);
  }

  std::cout << "noise: " << noise_->GetVariance();
  std::cout << " - p0: " << exp(p0) << " - t0: " << tau->GetParameter("Mean") << " - S_tau:" << tau->GetParameter("Variance");
  std::cout << " - v0: " << nu->GetParameter("Mean") << " - S_nu:" << nu->GetParameter("Variance");
  std::cout << " - MaxNu: " << nu_max << " - MinNu: " << nu_min;
  std::cout << " - MaxDelta: " << delta_max << " - MinDelta: " << delta_min << std::endl;

}

void FastNetworkModel::SaveCurrentState(unsigned int iter_num, const Realizations& simulated_real)
{
  /*
  unsigned int indiv_tot_num = simulated_real.at("Tau").size();
  std::ofstream outputs;
  std::string file_name = "/Users/igor.koval/Documents/Work/RiemAlzh/src/io/outputs/Simulated/Parameters" + std::to_string(iter_num) + ".txt";
  outputs.open(file_name, std::ofstream::out | std::ofstream::trunc);

  /// Save the final noise variance
  outputs << noise_->GetVariance() << std::endl;

 /// Save number of subject, Dimensions, number of Sources, number of control points
  outputs << indiv_tot_num << ", " << manifold_dim_ << ", " << indep_components_nb_  << ", " << control_points_nb_ << std::endl;

  /// Save p0, P0_mean and P0_var
  auto p0 = rand_var_.GetRandomVariable("P");
  outputs << simulated_real.at("P")(0) << ", " << p0->GetParameter("Mean") << ", " << p0->GetParameter("Variance") << std::endl;
  // std::cout  << simulated_real.at("P")(0) << ", " << p0->GetParameter("Mean") << ", " << p0->GetParameter("Variance") << std::endl;

  /// Save (ksi_i)
  for(size_t i = 0; i < indiv_tot_num; ++i)
  {
      outputs << simulated_real.at("Ksi")(i);
      if(i != indiv_tot_num - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save ksi_mean and ksi_Var
  auto ksi = rand_var_.GetRandomVariable("Ksi");
  outputs << ksi->GetParameter("Mean") << ", " << ksi->GetParameter("Variance") << std::endl;

  /// Save v0
  auto nu = rand_var_.GetRandomVariable("Nu");
  outputs << nu->GetParameter("Mean") << std::endl;

  /// Save (Tau_i)
  for(size_t i = 0; i < indiv_tot_num; ++i)
  {
      outputs << simulated_real.at("Tau")(i);
      if(i != indiv_tot_num - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save Tau_Mean and Tau_Var
  auto tau = rand_var_.GetRandomVariable("Tau");
  outputs << tau->GetParameter("Mean") << ", " << tau->GetParameter("Variance") << std::endl;

  /// Save (Delta_tilde_k)
  outputs << 0 << ", ";
  for(size_t i = 1; i < control_points_nb_; ++i)
  {
      outputs << simulated_real.at("Delta#" + std::to_string(i))(0);
      if(i != control_points_nb_ - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save (Delta_k)
  outputs << 0 << ", ";
  for(size_t i = 1; i < control_points_nb_; ++i)
  {
      std::string name = "Delta#" + std::to_string(i);
      outputs << rand_var_.GetRandomVariable(name)->GetParameter("Mean");
      if(i != control_points_nb_ - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save (Nu_k)
  for(size_t i = 0; i < control_points_nb_; ++i)
  {
      outputs << simulated_real.at("Nu", i);
      if(i != control_points_nb_ - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save (Nu_k_mean)
  for(size_t i = 0; i < control_points_nb_; ++i)
  {
      outputs << rand_var_.GetRandomVariable("Nu")->GetParameter("Mean");
      if(i != control_points_nb_ - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save (S_i)
  for(size_t i = 0; i < indiv_tot_num; ++i)
  {
      for(size_t j = 0; j < indep_components_nb_; ++j)
      {
          outputs << simulated_real.at("S#" + std::to_string(j))(i);
          if(i != indep_components_nb_ - 1) { outputs << ", "; }
      }
      outputs << std::endl;
  }

  /// Save (W_i)
  auto size_w = subjects_tot_num_;
  for(size_t i = 0; i < indiv_tot_num; ++i)
  {
      VectorType w_vec = space_shifts_.get_column(i);
      for(auto it = w_vec.begin(); it != w_vec.end(); ++it)
      {
          outputs << *it;
          if(i != size_w - 1) { outputs << ", "; }
      }
      outputs << std::endl;
  }
  */
}

void FastNetworkModel::SaveFinalState(const Realizations &reals, const Observations& obs) {}

void FastNetworkModel::SavePopulationFile() {}

void FastNetworkModel::SaveIndividualsFile(const Realizations &reals, const Observations& obs) {}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void FastNetworkModel::ComputeSubjectTimePoint(const Realizations &reals, const int indiv_num)
{
  if(indiv_num != -1)
  {
      double acc_factor = exp(reals.at("Ksi", indiv_num));
      double time_shift = reals.at("Tau", indiv_num);
      individual_time_points_[indiv_num] = acc_factor * (individual_obs_date_[indiv_num] - time_shift);
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


void FastNetworkModel::ComputeDeltas(const Realizations& reals)
{
  VectorType delta(control_points_nb_);
  delta(0) = 0.0;
  ScalarType * d = delta.memptr();

#pragma omp parallel for
  for(size_t i = 1; i < control_points_nb_; ++i)
      d[i] = reals.at("Delta#" + std::to_string(i), 0);
  

  auto interp_coeff = invert_kernel_matrix_ * delta;
  deltas_ = interpolation_matrix_ * interp_coeff;
}

void
FastNetworkModel
::ComputeNus(const Realizations& simulated_real)
{
  nus_ = interpolation_matrix_ * invert_kernel_matrix_ * simulated_real.at("Nu");

}

void
FastNetworkModel
::ComputeOrthonormalBasis()
{
  /// Get the data
  auto v0 = rand_var_.GetRandomVariable("Nu")->GetParameter("Mean");
  v0 = exp(v0);
  
  VectorType u_vec = (v0/ (init_pos_*init_pos_)) * nus_ % deltas_.exp();

  /// Compute the initial pivot vector u_vec
  double norm = u_vec.magnitude();
  u_vec(0) += copysign(1, -u_vec(0)) * norm;

  // Check the vectorise_row_wise function of the ArmadilloMatrixWrapper
  MatrixType u2(u_vec);
  double norm_u2 = u_vec.squared_magnitude();
  MatrixType final_matrix2 = (-2.0/norm_u2) * u2*u2.transpose();
  for(size_t i = 0; i < manifold_dim_; ++i)
    final_matrix2(i, i ) += 1;


  orthog_basis_ = final_matrix2;

}

void FastNetworkModel::ComputeAMatrix(const Realizations& simulated_real)
{

  MatrixType new_a(manifold_dim_, indep_components_nb_);

  for(int i = 0; i < indep_components_nb_; ++i)
  {
    VectorType beta(manifold_dim_, 0.0);
    for(size_t j = 0; j < manifold_dim_ - 1; ++j)
    {
      std::string number = std::to_string(int(j + i*(manifold_dim_ - 1)));
      beta(j) = simulated_real.at( "Beta#" + number, 0);
    }

    new_a.set_column(i, orthog_basis_ * beta);
  }

  a_matrix_ = new_a;

}

void FastNetworkModel::ComputeSpaceShifts(const Realizations& simulated_real)
{
  MatrixType sufficient_stat_vec(indep_components_nb_, subjects_tot_num_);

  for(int i = 0; i < indep_components_nb_; ++i)
    sufficient_stat_vec.set_row(i, simulated_real.at("S#" + std::to_string(i)));

  space_shifts_ = a_matrix_ * sufficient_stat_vec;
}


void FastNetworkModel::ComputeBlock1()
{
  ScalarType * d = deltas_.memptr();
  ScalarType *b = block1_.memptr();

#pragma omp simd
  for(size_t i = 0; i < manifold_dim_; ++i)
      b[i] = 1 / (init_pos_ * exp(d[i]));
}


void FastNetworkModel::ComputeBlock2()
{
  ScalarType * n = nus_.memptr();
  ScalarType * b = block2_.memptr();

#pragma omp simd
  for(size_t i = 0; i < manifold_dim_; ++i)
      b[i] = n[i] / init_pos_;
}


int FastNetworkModel::GetType(const MiniBlock &block_info) { 
  int type = std::get<0>(block_info[0]);
  
  for(auto it = block_info.begin(); it != block_info.end(); ++it) {
    
    int class_number = std::get<0>(*it);
    std::string real_name = std::get<1>(*it);
    int real_number = std::get<2>(*it);
    
    if(real_name == "P" || real_name == "Delta" || real_name == "Beta" || real_name == "Nu") {
      return -1;
    }
    if(type != real_number) {
      return -1;
    }
    
  }
  
  return type;
}

FastNetworkModel::VectorType FastNetworkModel::ComputeParallelCurve(int indiv_num, int ObservationNumber)
{
  double TimePoint = individual_time_points_[indiv_num](ObservationNumber);

  VectorType parallel_curve(manifold_dim_);


  ScalarType * d = deltas_.memptr();
  ScalarType * b1 = block1_.memptr();
  ScalarType * n = nus_.memptr();
  ScalarType * b2 = block2_.memptr();
  ScalarType * w = space_shifts_.get_column(indiv_num).memptr();
  ScalarType * p = parallel_curve.memptr();

//#pragma omp simd
  for(size_t i = 0; i < manifold_dim_; ++i)
  {
      p[i] = init_pos_ * exp(d[i] + w[i] * b1[i]  - b2[i] * TimePoint);
      //std::cout << "Observation at : " << indiv_obs_date_[indiv_num](ObservationNumber) << std::endl;
      //std::cout << p[i] << " = " << init_pos_ << " , " << d[i] << " , " << w[i] << " , " << n[i] << " , " << TimePoint << std::endl;
      //int a = 0;
  }


  //parallel_curve = init_pos_ * (deltas_ + space_shifts_.get_column(indiv_num) % block1_ - TimePoint * block2_).exp();

  return parallel_curve;
}
