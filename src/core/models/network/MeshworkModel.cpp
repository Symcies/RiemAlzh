#include "MeshworkModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

MeshworkModel::MeshworkModel(io::ModelSettings &model_settings)
{
  indep_components_nb_= model_settings.GetIndependentSourcesNumber();

  std::string kernel_mat_path = model_settings.GetInvertKernelPath();
  std::string interp_mat_path = model_settings.GetInterpolationKernelPath();

  invert_kernel_matrix_ = io::ReadData::OpenKernel(kernel_mat_path).transpose();
  interpolation_matrix_ = io::ReadData::OpenKernel(interp_mat_path);

  control_points_nb_ = invert_kernel_matrix_.columns();
  manifold_dim_ = interpolation_matrix_.rows();
  thickenesses_.set_size(manifold_dim_);
  deltas_.set_size(manifold_dim_);
  block1_.set_size(manifold_dim_);

}

MeshworkModel::~MeshworkModel()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void MeshworkModel::Initialize(const Observations& obs)
{
  /// Data-related attributes
  subjects_tot_num_  = obs.GetNumberOfSubjects();
  indiv_obs_date_    = obs.GetObservations();
  indiv_time_points_ = obs.GetObservations();
  obs_tot_num_       = obs.GetTotalNumberOfObservations();
  sum_obs_           = obs.GetTotalSumOfLandmarks();

  /// Population variables
  noise_ = std::make_shared<GaussianRandomVariable>( 0.0, 0.01 );

  rand_var_.AddRandomVariable("P", "Gaussian", {0.13, 0.00005 * 0.00005});
  asso_num_real_per_rand_var_["P"] = control_points_nb_;

  for(size_t i = 1; i < control_points_nb_; ++i)
  {
      std::string name = "Delta#" + std::to_string(i);
      rand_var_.AddRandomVariable(name, "Gaussian", {0, 0.003 * 0.003});
      asso_num_real_per_rand_var_[name] = 1;
  }

  for(size_t i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
  {
      std::string name = "Beta#" + std::to_string(i);
      rand_var_.AddRandomVariable(name, "Gaussian", {0, 0.001 * 0.001});
      asso_num_real_per_rand_var_[name] = 1;
  }

  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", {-3.1971, 0.000000004});
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;

  rand_var_.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;

  for(int i = 0; i < indep_components_nb_; ++i)
  {
      std::string name = "S#" + std::to_string(i);
      rand_var_.AddRandomVariable(name, "Gaussian", {0.0, 1});
      asso_num_real_per_rand_var_[name] =  subjects_tot_num_;
  }

}

void MeshworkModel::InitializeValidationDataParameters(const io::SimulatedDataSettings &data_settings,
                                                       const io::ModelSettings &model_settings) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
}


void UnivariateModel::UpdateModel(const Realizations &reals, const MiniBlock& block_info, const std::vector<std::string> names)
{
  bool compute_thickness = false;
  bool compute_delta = false;
  bool compute_basis = false;
  bool compute_a = false;
  bool compute_space_shift = false;
  bool compute_block = false;

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
      else if(name == "P")
      {
          indiv_only = false;
          compute_thickness = true;
          compute_basis = true;
          compute_a = true;
          compute_space_shift = true;
          compute_block = true;
          continue;
      }
      else if(name == "Delta")
      {
          indiv_only = false;
          compute_delta = true;
          compute_basis = true;
          compute_a = true;
          compute_space_shift = true;
          compute_block = true;
          continue;
      }
      else if(name == "Beta")
      {
          indiv_only = false;
          compute_a = true;
          compute_space_shift = true;
          continue;
      }
      else if(name == "S")
      {
          compute_space_shift = true;
      }
      else if(name == "All")
      {
          indiv_only = false;
          ComputeSubjectTimePoint(reals, -1);

          compute_thickness = true;
          compute_delta = true;
          compute_basis = true;
          compute_a = true;
          compute_space_shift = true;
          compute_block = true;
          break;
      }
      else
      {
          std::cerr << "The random variable name " << name << "is unknown to the meshwork model" << std::endl;
      }
  }

  if(indiv_only)          ComputeSubjectTimePoint(reals, type);
  if(compute_thickness)   ComputeThicknesses(reals);
  if(compute_delta)       ComputeDeltas(reals);
  if(compute_basis)       ComputeOrthonormalBasis();
  if(compute_a)           ComputeAMatrix(reals);
  if(compute_space_shift) ComputeSpaceShifts(reals);
  if(compute_block)       ComputeBlock();
}


MeshworkModel::SufficientStatisticsVector MeshworkModel::GetSufficientStatistics(
  const Realizations &reals, const Observations& obs)
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

  /// Sufficient Statistic Ksi and Ksi*Ksi
  VectorType s3 = reals.at("Ksi");
  VectorType s4 = reals.at("Ksi") % reals.at("Ksi");

  /// Sufficient statistic Tau and Tau*Tau
  VectorType s5 = reals.at("Tau");
  VectorType s6 = reals.at("Tau") % reals.at("Tau");

  /// Sufficient statistic beta_k
  VectorType s7((manifold_dim_-1) * indep_components_nb_);
  ScalarType * iter_s7 = s7.memptr();
  for(size_t i = 0; i < ((manifold_dim_-1) * indep_components_nb_); ++i){
      iter_s7[i] = reals.at("Beta#" + std::to_string(i), 0);
  }

  /// Sufficient statistic p_k and p_k*p_k
  VectorType s8 = reals.at("P");
  VectorType s9 = reals.at("P") % reals.at("P");

  /// Sufficient statistic delta_k
  VectorType s10(control_points_nb_ - 1);
  ScalarType * iter_s10 = s10.memptr();
  for(size_t i = 1; i < control_points_nb_; ++i){
      iter_s10[i-1] = reals.at("Delta#" + std::to_string(i), 0);
  }


  /// Return the sufficient statistic vector
  return {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10};
}

void MeshworkModel::UpdateRandomVariables(const SufficientStatisticsVector &stoch_sufficient_stats)
{
  /// Update the noise variance, sigma
  ScalarType noise_var = sum_obs_;
  const ScalarType * iter_s1 = stoch_sufficient_stats[0].memptr();
  const ScalarType * iter_s2 = stoch_sufficient_stats[1].memptr();
  for(size_t i = 0; i < stoch_sufficient_stats[0].size(); ++i){
      noise_var += - 2 * iter_s1[i] + iter_s2[i];
  }

  noise_var /= obs_tot_num_ * manifold_dim_;
  noise_->SetVariance(noise_var);

  /// Update Ksi : Mean and Variance
  ScalarType ksi_mean = 0.0, ksi_var = 0.0;
  const ScalarType * itS3 = stoch_sufficient_stats[2].memptr();
  const ScalarType * itS4 = stoch_sufficient_stats[3].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      ksi_mean     += itS3[i];
      ksi_var += itS4[i];
  }

  ksi_mean     /= subjects_tot_num_;
  ksi_var -= subjects_tot_num_ * ksi_mean * ksi_mean;
  ksi_var /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Ksi", {{"Mean", ksi_mean}, {"Variance", ksi_var}});


  /// Update Tau : Mean and Variance
  ScalarType tau_mean = 0.0, tau_var = 0.0;
  const ScalarType * itS5 = stoch_sufficient_stats[4].memptr();
  const ScalarType * itS6 = stoch_sufficient_stats[5].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      tau_mean     += itS5[i];
      tau_var += itS6[i];
  }

  tau_mean     /= subjects_tot_num_;
  tau_var -= subjects_tot_num_ * tau_mean * tau_mean;
  tau_var /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Tau", {{"Mean", tau_mean}, {"Variance", tau_var}});


  /// Update Beta_k : Mean
  const ScalarType * iter_s7 = stoch_sufficient_stats[6].memptr();
  for(size_t i = 0; i < stoch_sufficient_stats[6].size(); ++i){
      rand_var_.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", iter_s7[i]}});
  }

  /// Update P_k : Mean and Var
  ScalarType p_mean = 0.0, p_var = 0.0;
  const ScalarType * iter_s8 = stoch_sufficient_stats[7].memptr();
  const ScalarType * iter_s9 = stoch_sufficient_stats[8].memptr();

  for(size_t i = 0; i < control_points_nb_; ++i)
  {
      p_mean     += iter_s8[i];
      p_var += iter_s9[i];
  }

  p_mean     /= control_points_nb_;
  p_var -= control_points_nb_ * p_mean * p_mean;
  p_var /= control_points_nb_;

  rand_var_.UpdateRandomVariable("P", {{"Mean", p_mean}, {"Variance", p_var}});

  /// Update Delta_k : Mean
  const ScalarType * iter_s10 = stoch_sufficient_stats[9].memptr();
  for(size_t i = 0; i < control_points_nb_ - 1; ++i){
      rand_var_.UpdateRandomVariable("Delta#" + std::to_string(i+1), {{"Mean", iter_s10[i]}});
  }
}

Observations MeshworkModel::SimulateData(io::DataSettings &data_settings)
{
  typedef std::vector< std::pair< VectorType, double> > IndividualData;

  subjects_tot_num_ = data_settings.GetNumberOfSimulatedSubjects();

  /// Initialize the realizations and simulate them
  asso_num_real_per_rand_var_["P"] = control_points_nb_;

  for(int i = 1; i < control_points_nb_; ++i)
      asso_num_real_per_rand_var_["Delta#" + std::to_string(i)] = 1;


  for(size_t i = 0; i <  indep_components_nb_*(manifold_dim_ - 1); ++i)
      asso_num_real_per_rand_var_["Beta#" + std::to_string(i)] = 1;

  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;

  for(int i = 0; i < indep_components_nb_; ++i)
      asso_num_real_per_rand_var_["S#" + std::to_string(i)] = subjects_tot_num_;

  auto reals = SimulateRealizations();

  /// Update the model
  ComputeDeltas(reals);
  ComputeThicknesses(reals);
  ComputeOrthonormalBasis();
  ComputeAMatrix(reals);
  ComputeSpaceShifts(reals);

  ComputeBlock();

  /// Simulate the data
  std::random_device RD;
  std::mt19937 rand_num_gen(RD());
  std::uniform_int_distribution<int> Uni(data_settings.GetMinimumNumberOfObservations(), data_settings.GetMaximumNumberOfObservations());
  UniformRandomVariable time_points_num(60, 95);
  GaussianRandomVariable noise(0, noise_->GetVariance());

  /// Simulate the data
  Observations obs;
  for(int i = 0; i < subjects_tot_num_; ++i)
  {
    /// Get a random number of timepoints and sort them
    VectorType timepoints_vec = time_points_num.Samples(Uni(rand_num_gen));
    timepoints_vec.sort();
    indiv_time_points_.push_back(timepoints_vec);

    /// Simulate the data base on the time-points
    IndividualObservations indiv_obs(timepoints_vec);
    std::vector<VectorType> landmarks;
    for(size_t j = 0; j < timepoints_vec.size(); ++j)
    {
      landmarks.push_back(ComputeParallelCurve(i, j) + noise.Samples(manifold_dim_));
    }

    indiv_obs.AddLandmarks(landmarks);
    obs.AddIndividualData(indiv_obs);
  }

  /// Initialize the observation and model attributes
  obs.InitializeGlobalAttributes();
  indiv_obs_date_ = obs.GetObservations();
  sum_obs_           = obs.GetTotalSumOfLandmarks();
  obs_tot_num_     = obs.GetTotalNumberOfObservations();

  return obs;
}


std::vector<AbstractModel::MiniBlock> MeshworkModel::GetSamplerBlocks()
const
{
  int pop_type = -1;
  int nb_beta = 3;
  int nb_delta = 0;
  int nb_p = 5;


  std::vector<SamplerBlock> blocks;


  /// Insert P;
  MiniBlock block_p;

  for(size_t i = 0; i < control_points_nb_; ++i)
  {
      block_p.push_back(std::make_pair("P", i));
      //P.clear();
  }
  blocks.push_back(std::make_pair(pop_type, block_p));


  /// Insert Beta_k
  /*
  MiniBlock Beta;
  int BetaModulo = (int)indep_components_nb_*(manifold_dim_ - 1) /nb_beta;
  for(size_t i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
  {
      Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
      bool HingeCondition = (i%BetaModulo == 0 && i != 0);
      bool FinalCondition = (i == indep_components_nb_*(manifold_dim_ - 1) - 1);
      if(FinalCondition || HingeCondition)
      {
          blocks.push_back(std::make_pair(pop_type, Beta));
          Beta.clear();
      }
  }
   */
  MiniBlock block_beta;
  for(size_t i = 0; i < indep_components_nb_*(manifold_dim_ - 1); ++i)
  {
      block_beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
  }
  blocks.push_back(std::make_pair(pop_type, block_beta));

  /// Insert Delta_k
  MiniBlock block_delta;
  for(size_t i = 1; i < control_points_nb_; ++i)
  {
      block_delta.push_back(std::make_pair("Delta#" + std::to_string(i), 0));
  }
  blocks.push_back(std::make_pair(pop_type, block_delta));

  /// Individual variables
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      MiniBlock indiv_block;
      indiv_block.push_back(std::make_pair("Ksi", i));
      indiv_block.push_back(std::make_pair("Tau", i));
      for(size_t j = 0; j < indep_components_nb_; ++j)
          indiv_block.push_back(std::make_pair("S#" + std::to_string(j), i));

      blocks.push_back(std::make_pair(i, indiv_block));
  }


  return blocks;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
  

AbstractModel::VectorType MeshworkModel::ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info)
{
  int type = std::get<0>(block_info[0]);
  
  if(type == -1) {
    VectorType loglikelihood(subjects_tot_num_);
    ScalarType *l_ptr = ok.memptr();
    for (size_t i = 0; i < subjects_tot_num_; ++i)
      l_ptr[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);

    return loglikelihood;
  } else {
    return VectorType(1, ComputeIndividualLogLikelihood(obs.GetSubjectObservations(type), type));
  }
}

ScalarType MeshworkModel::ComputeIndividualLogLikelihood(const IndividualObservations& obs, const int subject_num)
{
  /// Get the data
  double log_likelihood = 0;
  auto num_time_points = obs.GetNumberOfTimePoints();

#pragma omp parallel for reduction(+:log_likelihood)
  for(size_t i = 0; i < num_time_points; ++i)
  {
      auto& it = obs.GetLandmark(i);
      VectorType parallel_curve = ComputeParallelCurve(indiv_num, i);
      log_likelihood += (it - parallel_curve).squared_magnitude();
  }

  log_likelihood /= -2*noise_->GetVariance();
  log_likelihood -= num_time_points * log(2 * noise_->GetVariance() * M_PI) / 2.0;

  return log_likelihood;
}

ScalarType MeshworkModel::GetPreviousLogLikelihood(const MiniBlock &block_info) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
}


void MeshworkModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
                                            const MiniBlock &block_info) {
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// outputs
////////////////////////////////////////////////////////////////////////////////////////////////////


void MeshworkModel::DisplayOutputs(const Realizations &reals)
{
  auto p = rand_var_.GetRandomVariable("P");
  auto p_real = reals.at("P");

  auto tau = rand_var_.GetRandomVariable("Tau");
  auto ksi = rand_var_.GetRandomVariable("Ksi");

  double delta_min = reals.at("Delta#1", 0);
  double delta_max = delta_min;
  for(size_t i = 1; i < 258; ++i)
  {
      double delta_k = reals.at("Delta#" + std::to_string(i), 0);
      delta_max = std::max(delta_max, delta_k);
      delta_min = std::min(delta_min, delta_k);
  }


  std::cout << "noise: " << noise_->GetVariance() << " - p_mean: " << exp(p->GetParameter("Mean"));
  std::cout << " - PVar: " << p->GetParameter("Variance") << " - PMin: " << exp(p_real.min_value()) << " - PMax: " << exp(p_real.max_value());
  std::cout << " - T0: " << tau->GetParameter("Mean") << " - TauVar: " << tau->GetParameter("Variance");
  std::cout << " - v0: " << exp(ksi->GetParameter("Mean")) << " - KsiVar: " << ksi->GetParameter("Variance");
  std::cout << " - MinDelta: " << delta_min << " - MaxDelta: " << delta_max << std::endl;
}


void MeshworkModel::SaveData(unsigned int iter_num, const Realizations &reals)
{

  std::ofstream outputs;
  std::string file_name = "/Users/igor.koval/Documents/Work/RiemAlzh/src/io/outputs/Meshwork/Parameters" + std::to_string(iter_num) + ".txt";
  outputs.open(file_name, std::ofstream::out | std::ofstream::trunc);

  /// Save the final noise variance
  outputs << noise_->GetVariance() << std::endl;

  /// Save the number of subjects, the manifold dimension, the number of sources, and, the number of control points
  outputs << subjects_tot_num_ << ", " << manifold_dim_ << ", " << indep_components_nb_ << ", " << control_points_nb_ << std::endl;

  /// Save the delta_mean -> First one being equal to 0 as the reference
  outputs << 0 << ", ";
  for(size_t i = 1; i < control_points_nb_; ++i)
  {
      outputs << rand_var_.GetRandomVariable("Delta#" + std::to_string(i))->GetParameter("Mean");
      if(i != control_points_nb_ - 1) { outputs << ", "; }
  }
  outputs << std::endl;

  /// Save the thicknesses
  for(size_t i = 0; i < control_points_nb_; ++i)
  {
      outputs << reals.at("P", i);
      if(i != control_points_nb_ - 1) { outputs << ", ";}
  }
  outputs << std::endl;

  /// Save the tau
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      outputs << reals.at("Tau", i) ;
      if(i != subjects_tot_num_ - 1) { outputs << ", ";}
  }
  outputs << std::endl;

  /// Save the ksi
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      outputs << reals.at("Ksi", i) ;
      if(i != subjects_tot_num_ - 1) { outputs << ", ";}
  }

      /// Save (S_i)
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      for(size_t j = 0; j < indep_components_nb_; ++j)
      {
          outputs << reals.at("S#" + std::to_string(j))(i);
          if(i != indep_components_nb_ - 1) { outputs << ", "; }
      }
      outputs << std::endl;
  }

  /// Save (W_i)
  auto size_w = subjects_tot_num_;
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      VectorType w = space_shifts_.get_column(i);
      for(auto it = w.begin(); it != w.end(); ++it)
      {
          outputs << *it;
          if(i != size_w - 1) { outputs << ", "; }
      }
      outputs << std::endl;
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void MeshworkModel::ComputeSubjectTimePoint(const Realizations &reals, const int indiv_num)
{
  if(indiv_num != -1)
  {
      double acc_factor = exp(reals.at("Ksi", indiv_num));
      double time_shift = reals.at("Tau", indiv_num);
      indiv_time_points_[indiv_num] = acc_factor * (indiv_obs_date_[indiv_num] - time_shift);
  }
  else
  {

      for(size_t i = 0; i < subjects_tot_num_; ++i)
      {
          double acc_factor = exp(reals.at("Ksi")(i));
          double time_shift = reals.at("Tau")(i);

          indiv_time_points_[i] = acc_factor * (indiv_obs_date_[i] - time_shift);
      }
  }
}

void MeshworkModel::ComputeDeltas(const Realizations &reals)
{
  VectorType delta(control_points_nb_);
  delta(0) = 0.0;
  ScalarType * d = delta.memptr();

#pragma omp parallel for
  for(size_t i = 1; i < control_points_nb_; ++i)
  {
      d[i] = reals.at("Delta#" + std::to_string(i), 0);
  }

  deltas_ = interpolation_matrix_ * invert_kernel_matrix_ * delta;
}


void MeshworkModel::ComputeThicknesses(const Realizations &reals)
{
  thickenesses_ = interpolation_matrix_ * invert_kernel_matrix_ * reals.at("P").exp();
}

void MeshworkModel::ComputeOrthonormalBasis()
{
      /// Get the data
  auto v0 = rand_var_.GetRandomVariable("Ksi")->GetParameter("Mean");
  v0 = exp(v0);


  VectorType u_vec(manifold_dim_);
  ScalarType * t = thickenesses_.memptr();
  ScalarType * d = deltas_.memptr();
  ScalarType * u = u_vec.memptr();

#pragma omp simd
  for(int i = 0; i < manifold_dim_; ++i)
  {
      u[i] = v0 / (t[i] * t[i])* exp(-d[i]);
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

void MeshworkModel::ComputeAMatrix(const Realizations &reals)
{

  MatrixType new_a(manifold_dim_, indep_components_nb_);

  for(int i = 0; i < indep_components_nb_; ++i)
  {
      VectorType beta(manifold_dim_, 0.0);
      for(size_t j = 0; j < manifold_dim_ - 1; ++j)
      {
          std::string number = std::to_string(int(j + i*(manifold_dim_ - 1)));
          beta(j) = reals.at( "Beta#" + number, 0);
      }

      new_a.set_column(i, orthog_basis_ * beta);
  }

  a_matrix_ = new_a;
}

void MeshworkModel::ComputeSpaceShifts(const Realizations &reals)
{

  MatrixType stoch_sufficient_stats(indep_components_nb_, subjects_tot_num_);
  for(int i = 0; i < indep_components_nb_; ++i){
      stoch_sufficient_stats.set_row(i, reals.at("S#" + std::to_string(i)));
  }
  space_shifts_ = a_matrix_ * stoch_sufficient_stats;
}

void MeshworkModel::ComputeBlock()
{
  ScalarType * t = thickenesses_.memptr();
  ScalarType * d = deltas_.memptr();
  ScalarType * b = block1_.memptr();

#pragma omp simd
  for(size_t i = 0; i < manifold_dim_; ++i){
      b[i] = 1 / (t[i] * exp(d[i]));
  }
}


AbstractModel::VectorType MeshworkModel::ComputeParallelCurve(int indiv_num, int obs_num)
{
  double time_point = indiv_time_points_[indiv_num](obs_num);

  VectorType parallel_curve(manifold_dim_);

  ScalarType * p = parallel_curve.memptr();
  ScalarType * t = thickenesses_.memptr();
  ScalarType * d = deltas_.memptr();
  ScalarType * b = block1_.memptr();
  ScalarType * w = space_shifts_.get_column(indiv_num).memptr();

  for(size_t i = 0; i < manifold_dim_; ++i){
      p[i] = t[i] * exp( w[i] / (t[i]*exp(d[i])) + d[i] - time_point/t[i]);
  }

  return parallel_curve;
}
