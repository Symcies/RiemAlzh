#include "NetworkModel.h"

NetworkModel::NetworkModel(io::ModelSettings &model_settings)
{
  indep_components_nb_= model_settings.GetIndependentSourcesNumber();
  acceptance_ratio_to_display_ = model_settings.GetAcceptanceRatioToDisplay();
  output_file_name_ = model_settings.GetOutputFileName();
  
  std::string kernel_mat_path = model_settings.GetInvertKernelPath();
  std::string interp_mat_path = model_settings.GetInterpolationKernelPath();

  invert_kernel_matrix_ = io::ReadData::OpenKernel(kernel_mat_path).transpose();
  interpolation_matrix_ = io::ReadData::OpenKernel(interp_mat_path);

  control_points_nb_ = invert_kernel_matrix_.columns();
  manifold_dim_ = interpolation_matrix_.rows();
  thickenesses_.set_size(manifold_dim_);
  nus_.set_size(manifold_dim_);
  block1_.set_size(manifold_dim_);
  
  std::remove((GV::BUILD_DIR + output_file_name_ + "pop_params.txt").c_str());
  std::remove((GV::BUILD_DIR + output_file_name_ + "indiv_params.txt").c_str());
  std::remove((GV::BUILD_DIR + output_file_name_ + "convergence_params.txt").c_str());
}

NetworkModel::~NetworkModel()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void NetworkModel::Initialize(const Observations& obs)
{
  /// Data-related attributes
  subjects_tot_num_       = obs.GetNumberOfSubjects();
  individual_obs_date_    = obs.GetObservations();
  individual_time_points_ = obs.GetObservations();
  obs_tot_num_            = obs.GetTotalNumberOfObservations();
  sum_obs_                = obs.GetTotalSumOfLandmarks();
  
  /// noise
  noise_ = std::make_shared<GaussianRandomVariable>( rv_params_.at("noise").first[0], rv_params_.at("noise").first[1] );

  /// Population variables
  rand_var_.AddRandomVariable("P", "Gaussian", rv_params_.at("P").first);
  asso_num_real_per_rand_var_["P"] = control_points_nb_;
  proposition_distribution_variance_["P"] = rv_params_.at("P").second;

  rand_var_.AddRandomVariable("Nu", "Gaussian", rv_params_.at("Nu").first);
  asso_num_real_per_rand_var_["Nu"] = control_points_nb_;
  proposition_distribution_variance_["Nu"] = rv_params_.at("Nu").second;

  for(int i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i) {
      std::string name = "Beta#" + std::to_string(i);
      rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Beta").first);
      asso_num_real_per_rand_var_.insert({name, 1});
  }
  proposition_distribution_variance_["Beta"] = rv_params_.at("Beta").second;


  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", rv_params_.at("Ksi").first);
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  proposition_distribution_variance_["Ksi"] = rv_params_.at("Ksi").second;

  rand_var_.AddRandomVariable("Tau", "Gaussian", rv_params_.at("Tau").first);
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
  proposition_distribution_variance_["Tau"] = rv_params_.at("Tau").second;

  for(int i = 0; i < indep_components_nb_; ++i) {
      std::string name = "S#" + std::to_string(i);
      rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("S").first);
      asso_num_real_per_rand_var_[name] = subjects_tot_num_;
  }
  proposition_distribution_variance_["S"] = rv_params_.at("S").second;
}

void NetworkModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                      const io::ModelSettings &model_settings) {
  /// Initialize the model
  // TODO : Check which one has to be updated --> Some (a lot) are initialized in the constructor
  // DO AN ASSERT : manifold_dim_ = data_settings.GetDimensionOfSimulatedObservations();
  indep_components_nb_ = model_settings.GetIndependentSourcesNumber();
  last_loglikelihood_.set_size(subjects_tot_num_);
  
  /// Population variables
  rand_var_.AddRandomVariable("P", "Gaussian", rv_params_.at("P").first);
  asso_num_real_per_rand_var_["P"] = 1;
  
  rand_var_.AddRandomVariable("Nu", "Gaussian", rv_params_.at("Nu").first);
  asso_num_real_per_rand_var_["Nu"] = control_points_nb_;
  
  for(int i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i) {
    std::string name = "Beta#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("Beta").first);
    asso_num_real_per_rand_var_[name] = 1;
  }
  

  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", rv_params_.at("Ksi").first);
  rand_var_.AddRandomVariable("Tau", "Gaussian", rv_params_.at("Tau").first); 
  
  for(int i = 0; i < indep_components_nb_; ++i) {
    std::string name = "S#" + std::to_string(i);
    rand_var_.AddRandomVariable(name, "Gaussian", rv_params_.at("S").first);
  }
}

void NetworkModel::UpdateModel(const Realizations& reals, const MiniBlock& block_info, const std::vector<std::string> names)
{
  
  int type = GetType(block_info);  
  
  bool compute_thickness = false;
  bool compute_nu = false;
  bool compute_block1 = false;
  bool compute_basis = false;
  bool compute_a = false;
  bool compute_space = false;

  bool indiv_only = (type > -1);

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
          indiv_only = false;
          compute_nu = true;
          compute_block1 = true;
          compute_basis = true;
          compute_a = true;
          compute_space = true;
      }
      else if(name == "P")
      {
          indiv_only = false;
          compute_thickness = true;
          compute_block1 = true;
          compute_basis = true;
          compute_a = true;
          compute_space = true;
      }
      else if(name == "Beta")
      {
          indiv_only = false;
          compute_a = true;
          compute_space = true;
          continue;
      }
      else if(name == "S")
      {
          compute_space= true;
      }
      else if(name == "All")
      {
          indiv_only = false;
          ComputeSubjectTimePoint(reals, -1);

          compute_thickness = true;
          compute_nu = true;
          compute_basis = true;
          compute_a = true;
          compute_space = true;
          compute_block1 = true;
          break;
      }
      else
      {
          std::cerr << "The random variable name " << name << "is unknown to the meshwork model" << std::endl;
      }
  }


  if(indiv_only)        ComputeSubjectTimePoint(reals, type);
  if(compute_thickness) ComputeThicknesses(reals);
  if(compute_nu)        ComputeNus(reals);
  if(compute_block1)    ComputeBlock();
  if(compute_basis)     ComputeOrthonormalBasis();
  if(compute_a)         ComputeAMatrix(reals);
  if(compute_space)     ComputeSpaceShifts(reals);
}

AbstractModel::SufficientStatisticsVector NetworkModel::GetSufficientStatistics(const Realizations &reals, const Observations& obs)
{
  /// s1 <- y_ij * eta_ij    &    s2 <- eta_ij * eta_ij
  VectorType s1(obs_tot_num_), s2(obs_tot_num_);
  auto iter_s1 = s1.begin(), iter_s2 = s2.begin();
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      for(size_t j = 0; j < obs.GetNumberOfTimePoints(i); ++j)
      {
          VectorType parallel_curve = ComputeParallelCurve(i, j);
          *iter_s1 = dot_product(parallel_curve, obs.GetSubjectLandmark(i, j));
          *iter_s2 = parallel_curve.squared_magnitude();
          ++iter_s1, ++iter_s2;
      }
  }

  /// Sufficient Statistic P_k and P_k * P_k
  VectorType s3 = reals.at("P");
  VectorType s4 = reals.at("P") % reals.at("P");

  /// Sufficient Statistic Nu_k and Nu_k*Nu_k
  VectorType s5 = reals.at("Nu");
  VectorType s6 = reals.at("Nu") % reals.at("Nu");

  /// Sufficient statistic beta_k
  VectorType s7((manifold_dim_-1) * indep_components_nb_);
  ScalarType * iter_s7 = s7.memptr();
  for(size_t i = 0; i < ((manifold_dim_-1) * indep_components_nb_); ++i)
      iter_s7[i] = reals.at("Beta#" + std::to_string(i), 0);

  /// Sufficient statistic Ksi_i * Ksi_i
  VectorType s8 = reals.at("Ksi") % reals.at("Ksi");
  
  /// Sufficient statistic Tau_i and Tau_i * Tau_i
  VectorType s9 = reals.at("Tau");
  VectorType s10 = reals.at("Tau") % reals.at("Tau");

  return {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10};
}

void NetworkModel::UpdateRandomVariables(const SufficientStatisticsVector &stoch_sufficient_stats)
{
  /// Update the noise variance, sigma
  ScalarType noise_var = sum_obs_;
  const ScalarType * iter_s1 = stoch_sufficient_stats[0].memptr();
  const ScalarType * iter_s2 = stoch_sufficient_stats[1].memptr();
  for(size_t i = 0; i < stoch_sufficient_stats[0].size(); ++i)
      noise_var += - 2 * iter_s1[i] + iter_s2[i];

  noise_var /= obs_tot_num_ * manifold_dim_;
  noise_->SetVariance(noise_var);


  /// Update P_k : Mean and Var
  ScalarType p_mean = 0.0, p_var = 0.0;
  const ScalarType * iter_s3 = stoch_sufficient_stats[2].memptr();
  const ScalarType * iter_s4 = stoch_sufficient_stats[3].memptr();

  for(size_t i = 0; i < control_points_nb_; ++i)
  {
    p_mean += iter_s3[i];
    p_var  += iter_s4[i];
  }

  p_mean /= control_points_nb_;
  p_var  -= control_points_nb_ * p_mean * p_mean;
  p_var  /= control_points_nb_;

  rand_var_.UpdateRandomVariable("P", {{"Mean", p_mean}, {"Variance", p_var}});

  /// Update Nu_k : Mean and Var
  ScalarType nu_mean = 0.0, nu_var = 0.0;
  const ScalarType * iter_s5 = stoch_sufficient_stats[4].memptr();
  const ScalarType * iter_s6 = stoch_sufficient_stats[5].memptr();

  for(size_t i = 0; i < control_points_nb_; ++i)
  {
    nu_mean += iter_s5[i];
    nu_var  += iter_s6[i];
  }

  nu_mean /= control_points_nb_;
  nu_var  -= control_points_nb_ * nu_mean * nu_mean;
  nu_var  /= control_points_nb_;

  rand_var_.UpdateRandomVariable("Nu", {{"Mean", nu_mean}, {"Variance", nu_var}});

  /// Update Beta_k : Mean
  const ScalarType * iter_s7 = stoch_sufficient_stats[6].memptr();
  for(size_t i = 0; i < stoch_sufficient_stats[6].size(); ++i){
    rand_var_.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", iter_s7[i]}});
  }

  /// Update Ksi : Mean and Variance
  ScalarType ksi_var = 0.0;
  const ScalarType * iter_s8 = stoch_sufficient_stats[7].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i) {
    ksi_var += iter_s4[i];
  }

  ksi_var /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Ksi", {{"Variance", ksi_var}});


  /// Update Tau : Mean and Variance
  ScalarType tau_mean = 0.0, tau_var = 0.0;
  const ScalarType * iter_s9  = stoch_sufficient_stats[8].memptr();
  const ScalarType * itS10 = stoch_sufficient_stats[9].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      tau_mean += iter_s9[i];
      tau_var  += itS10[i];
  }

  tau_mean /= subjects_tot_num_;
  tau_var  -= subjects_tot_num_ * tau_mean * tau_mean;
  tau_var  /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Tau", {{"Mean", tau_mean}, {"Variance", tau_var}});

}


Observations NetworkModel::SimulateData(io::SimulatedDataSettings &data_settings)
{
  individual_obs_date_.clear();
  individual_time_points_.clear();
  
  subjects_tot_num_ = data_settings.GetNumberOfSimulatedSubjects();

  /// Initialize the realizations and simulate them
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
  for(int i = 0; i < indep_components_nb_; ++i)
      asso_num_real_per_rand_var_["S#" + std::to_string(i)] = subjects_tot_num_;

  auto reals = SimulateRealizations();

  /// Update the model
  ComputeThicknesses(reals);
  ComputeNus(reals);
  ComputeBlock();
  ComputeOrthonormalBasis();
  ComputeAMatrix(reals);
  ComputeSpaceShifts(reals);

  /// Simulate the data
  std::random_device random_device;
  std::mt19937 rand_num_gen(random_device());
  std::uniform_int_distribution<int> uniform_dist(data_settings.GetMinimumNumberOfObservations(), data_settings.GetMaximumNumberOfObservations());
  UniformRandomVariable time_points_num(60, 95);

  /// Simulate the data
  Observations obs;
  for(int i = 0; i < subjects_tot_num_; ++i)
  {
    /// Get a random number of timepoints and sort them
    VectorType time_points = time_points_num.Samples(uniform_dist(rand_num_gen));
    time_points.sort();
    individual_time_points_.push_back(time_points);
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

std::vector<AbstractModel::MiniBlock> NetworkModel::GetSamplerBlocks(unsigned int blocks_number) const
{
  int population_type = -1;

  std::vector<MiniBlock> blocks;

  /// Insert P
  MiniBlock block_p;
  for(size_t i = 0; i < control_points_nb_; ++i)
      block_p.push_back(std::tuple<int, std::string, int>(population_type, "P", i));
  blocks.push_back(block_p);

  /// Insert Nu
  MiniBlock block_nu;
  for(size_t i = 0; i < control_points_nb_; ++i)
      block_nu.push_back(std::tuple<int, std::string, int>(population_type, "Nu", i));
  blocks.push_back(block_nu);

  /// Insert Beta
  MiniBlock block_beta;
  for(size_t i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
    block_beta.push_back(std::tuple<int, std::string, int>(population_type, "Beta#" + std::to_string(i), 0));
  blocks.push_back(block_beta);

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

void NetworkModel::InitializeLogLikelihood(const Observations &obs) {
  
  last_loglikelihood_.set_size(subjects_tot_num_);
  ScalarType * l = last_loglikelihood_.memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i) 
    l[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);
  
}


AbstractModel::VectorType NetworkModel::ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info)
{
  int type = GetType(block_info);
  
  if(type == -1) {
    VectorType loglikelihood(subjects_tot_num_);
    ScalarType * l_ptr = loglikelihood.memptr();
    
    for (size_t i = 0; i < subjects_tot_num_; ++i)
      l_ptr[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);

    return loglikelihood;
  } else {
    return VectorType(1, ComputeIndividualLogLikelihood(obs.GetSubjectObservations(type), type));
  }
}

ScalarType NetworkModel::ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int subject_num)
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

ScalarType NetworkModel::GetPreviousLogLikelihood(const MiniBlock &block_info) {
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


void NetworkModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
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
/// outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void NetworkModel::DisplayOutputs(const Realizations &reals)
{
  auto p = rand_var_.GetRandomVariable("P");
  auto p_real = reals.at("P");

  auto nu = rand_var_.GetRandomVariable("Nu");
  auto nu_real = reals.at("Nu");

  auto tau = rand_var_.GetRandomVariable("Tau");
  auto ksi = rand_var_.GetRandomVariable("Ksi");

  std::cout << "noise: " << noise_->GetVariance() << " - p_mean: " << exp(p->GetParameter("Mean"));
  std::cout << " - PVar: " << p->GetParameter("Variance") << " - PMin: " << exp(p_real.min_value()) << " - PMax: " << exp(p_real.max_value());
  std::cout << " - T0: " << tau->GetParameter("Mean") << " - TauVar: " << tau->GetParameter("Variance");
  std::cout << " - KsiVar: " << ksi->GetParameter("Variance");
  std::cout << " - nu_mean/V0: " << nu->GetParameter("Mean") << " - nu_var: " << nu->GetParameter("Variance");
  std::cout << " - NuMin: " << nu_real.min_value() << " - NuMax: " << nu_real.max_value() << std::endl;

}

void NetworkModel::SaveCurrentState(unsigned int iter_num, const Realizations &reals)
{
  
    std::ofstream log_file;
  
  if(!GV::TEST_RUN) {
    log_file.open(GV::BUILD_DIR + output_file_name_ + "convergence_params.txt" , std::ofstream::out | std::ofstream::app );
    
    auto p     = rand_var_.GetRandomVariable("P");
    auto nu    = rand_var_.GetRandomVariable("Nu");
    auto tau   = rand_var_.GetRandomVariable("Tau");
    auto ksi   = rand_var_.GetRandomVariable("Ksi");
    auto beta  = rand_var_.GetRandomVariable("Beta#1")->GetParameter("Mean");
    auto s     = rand_var_.GetRandomVariable("S#1")->GetParameter("Mean");
    
    
    
    if (iter_num == 0){
      log_file << "iter Noise PMean PVar NuMean NuVar TauMean TauVar KsiMean KsiVar Beta#1 S#1" << std::endl;
    }
    
    log_file << iter_num << " " << noise_->GetVariance() << " " 
             << p->GetParameter("Mean")   << " " << p->GetParameter("Variance")   << " "
             << nu->GetParameter("Mean")  << " " << nu->GetParameter("Variance")  << " "
             << tau->GetParameter("Mean") << " " << tau->GetParameter("Variance") << " "
             << ksi->GetParameter("Mean") << " " << ksi->GetParameter("Variance") << " "
             << beta << " " << s << std::endl;
             
    
    log_file.close();
  }
  else {
    // TODO TODO TODO TODO TODO
    // TODO TODO TODO TODO TODO
    // TODO TODO TODO TODO TODO
  }

}

void NetworkModel::SavePopulationFile() {
  std::ofstream log_file;
  log_file.open(GV::BUILD_DIR + output_file_name_ + "pop_params.txt", std::ofstream::out | std::ofstream::app);
  
  log_file << "Noise " << noise_->GetVariance() << std::endl;
  
  log_file << "P ";
  for (size_t i = 0; i < manifold_dim_; ++i) {
    log_file << thickenesses_[i] << " ";
  }
  log_file << std::endl;
  
  log_file << "Nu ";
  for (size_t i = 0; i < manifold_dim_; ++i) {
    log_file << nus_[i] << " ";
  }
  log_file << std::endl;
  
  log_file << "V0 " << exp(rand_var_.GetRandomVariable("Ksi")->GetParameter("Mean")) << std::endl;
  
  log_file << "T0 " << rand_var_.GetRandomVariable("Tau")->GetParameter("Mean") << std::endl;
  
  log_file.close();
}


void NetworkModel::SaveIndividualsFile(const Realizations &reals, const Observations &obs) {
  std::ofstream log_file;
  
  log_file.open(GV::BUILD_DIR + output_file_name_ + "indiv_params.txt", std::ofstream::out | std::ofstream::app);
  
  log_file << "ID Tau Ksi W" << std::endl;
  log_file << "1 1 1 " << manifold_dim_ << std::endl;
  
  for (size_t i = 0; i < subjects_tot_num_; ++i) {
    log_file << obs.GetId(i) << " " << reals.at("Tau")(i) << " " << reals.at("Ksi")(i) << " ";
    for (size_t j = 0; j < manifold_dim_; ++j) {
      log_file << space_shifts_.get_column(i)(j) << " ";
    }
    log_file << std::endl;
  }
  
  log_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void NetworkModel::ComputeSubjectTimePoint(const Realizations &reals, const int indiv_num)
{
  if(indiv_num != -1)
  {
      ScalarType acc_factor = exp(reals.at("Ksi", indiv_num));
      ScalarType time_shift = reals.at("Tau", indiv_num);
      individual_time_points_[indiv_num] = acc_factor * (individual_obs_date_[indiv_num] - time_shift);
  }
  else
  {

      for(size_t i = 0; i < subjects_tot_num_; ++i)
      {
          ScalarType acc_factor = exp(reals.at("Ksi")(i));
          ScalarType time_shift = reals.at("Tau")(i);

          individual_time_points_[i] = acc_factor * (individual_obs_date_[i] - time_shift);
      }
  }
}


void NetworkModel::ComputeNus(const Realizations &reals)
{
  nus_ = interpolation_matrix_ * invert_kernel_matrix_ * reals.at("Nu");
}

void NetworkModel::ComputeThicknesses(const Realizations &reals)
{
  thickenesses_ = interpolation_matrix_ * invert_kernel_matrix_ * reals.at("P").exp();
}

void NetworkModel::ComputeOrthonormalBasis()
{

  VectorType u_vec(manifold_dim_);
  ScalarType * u = u_vec.memptr();
  ScalarType * t = thickenesses_.memptr();
  ScalarType * n = nus_.memptr();
  
  for(size_t i = 0; i < manifold_dim_; ++i){
      u[i] = n[i] / (t[i] * t[i]);
  }

  /// Compute the initial pivot vector u
  double norm = u_vec.magnitude();
  u_vec(0) += copysign(1, -u_vec(0)) * norm;

  // Check the vectorise_row_wise function of the ArmadilloMatrixWrapper
  MatrixType u2_mat(u_vec);
  double norm_u2 = u_vec.squared_magnitude();
  MatrixType final_mat2 = (-2.0/norm_u2) * u2_mat*u2_mat.transpose();
  for(size_t i = 0; i < manifold_dim_; ++i){
      final_mat2(i, i ) += 1;
  }


  orthog_basis_ = final_mat2;
}


void NetworkModel::ComputeAMatrix(const Realizations &reals)
{
  MatrixType new_a_matrix(manifold_dim_, indep_components_nb_);

  for(int i = 0; i < indep_components_nb_; ++i)
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


void NetworkModel::ComputeSpaceShifts(const Realizations &reals)
{
  MatrixType stoch_sufficient_stats(indep_components_nb_, subjects_tot_num_);
  for(int i = 0; i < indep_components_nb_; ++i)
  {
      stoch_sufficient_stats.set_row(i, reals.at("S#" + std::to_string(i)));
  }
  space_shifts_ = a_matrix_ * stoch_sufficient_stats;
}

void NetworkModel::ComputeBlock()
{
  ScalarType * b = block1_.memptr();
  ScalarType * n = nus_.memptr();
  ScalarType * t = thickenesses_.memptr();

  for(size_t i = 0; i < manifold_dim_; ++i) {
      b[i] = n[i] / t[i] ;
  }
}

int NetworkModel::GetType(const MiniBlock &block_info) {
  int type = std::get<2>(block_info[0]);
  
  for(auto it = block_info.begin(); it != block_info.end(); ++it) {
    
    int class_number = std::get<0>(*it);
    std::string real_name = std::get<1>(*it);
    int real_number = std::get<2>(*it);
    
    if(real_name == "P" or real_name == "Nu" or real_name == "Beta") {
      return -1;
    }
    if(type != real_number) {
      return -1;
    }
    
  }
  
  return type;
}

AbstractModel::VectorType NetworkModel::ComputeParallelCurve(int indiv_num, int obs_num)
{
  double time_point = individual_time_points_[indiv_num](obs_num);
  
  VectorType parallel_curve(manifold_dim_);
  
  ScalarType * p = parallel_curve.memptr();
  ScalarType * t = thickenesses_.memptr();
  ScalarType * w = space_shifts_.get_column(indiv_num).memptr();
  ScalarType * b = block1_.memptr();

  for(size_t i = 0; i < manifold_dim_; ++i)
      p[i] = t[i] * exp(w[i] / t[i] - b[i] * time_point);
  
  return parallel_curve;
}
