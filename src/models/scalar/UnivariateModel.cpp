#include "UnivariateModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

//TODO(igor): why empty and param?
UnivariateModel::UnivariateModel(io::ModelSettings &model_settings)
{
}

UnivariateModel::~UnivariateModel()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void UnivariateModel::Initialize(const Observations &obs)
{
  /// Data-related attributes
  manifold_dim_        = obs.GetSubjectObservations(0).GetCognitiveScore(0).size();
  subjects_tot_num_    = obs.GetNumberOfSubjects();
  individual_obs_date_ = obs.GetObservations();
  subj_time_points_    = obs.GetObservations();
  obs_tot_num_         = obs.GetTotalNumberOfObservations();
  sum_obs_             = obs.GetTotalSumOfCognitiveScores();

  /// Population variables
  /// (Initialization of the random variable m_Random Variable
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  noise_ = std::make_shared<GaussianRandomVariable>(0.0, 0.00001);
  last_loglikelihood_.set_size(subjects_tot_num_);

  rand_var_.AddRandomVariable("P", "Gaussian", {0.02, 0.0001*0.0001});
  asso_num_real_per_rand_var_["P"] = 1;

  /// Individual realizations
  /// (Initialization of the random variable m_Random Variable
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  rand_var_.AddRandomVariable("Ksi", "Gaussian", {0.13, 0.0004});
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;

  rand_var_.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
}


ScalarType UnivariateModel::InitializePropositionDistributionVariance(std::string name) const
{
  name = name.substr(0, name.find_first_of("#"));

  /// It returns the variance of the proposition distribution
  /// of a given random variable of the model
  if("P" == name)
    return 0.00001;
  if("Ksi" == name)
    return 0.0001;
  if("Tau" == name)
    return 0.5;
}

void UnivariateModel::UpdateModel(const Realizations &reals, const MiniBlock& block_info, const std::vector<std::string> names)
{
  /// Given a list of names (which in fact corresponds to the variables that have potentially changed),
  /// the function updates the parameters associated to these names

  int type = std::get<0>(block_info[0]);

  /// Possible parameters to update, depending on the name being in "vect<> names"
  //TODO: init parameters?
  bool compute_position;
  bool compute_individual = (type > -1);

  /// Parameters to update, depending on the names called
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
      compute_individual = false;
      compute_position = true;
    }
    else if(name == "All")
    {
      ComputeSubjectTimePoint(reals, -1);
      compute_individual = false;
      compute_position = true;
    }
    else
    {
      std::cerr << "The realization does not exist in the multivariate model > update model" << std::endl;
    }
  }
  
  if(compute_individual) {
    ComputeSubjectTimePoint(reals, type);
  }
  if(compute_position) {
    position_ = 1.0 / (1 + exp(-reals.at("P", 0)));
  }
}

AbstractModel::SufficientStatisticsVector UnivariateModel::GetSufficientStatistics(const Realizations& reals, const Observations &obs)
{
  /// Computation of the suffisient statistics of the model
  /// Basically, it is a vector (named s1 to S9) of vectors

  /// s1 <- y_ij * eta_ij    &    s2 <- eta_ij * eta_ij
  VectorType s1(obs_tot_num_), s2(obs_tot_num_);
  auto it_s1 = s1.begin(), it_s2 = s2.begin();
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
    for(size_t j = 0; j < obs.GetNumberOfTimePoints(i); ++j)
    {
      VectorType parallel_curve = ComputeParallelCurve(i, j);
      *it_s1 = parallel_curve(0) * obs.GetSubjectCognitiveScore(i, j)(0);
      *it_s2 = parallel_curve(0) * parallel_curve(0);
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
  VectorType s7(1, reals.at("P", 0));

  return {s1, s2, s3, s4, s5, s6, s7};
}

void UnivariateModel::UpdateRandomVariables(const SufficientStatisticsVector &stoch_sufficient_stats)
{
    /// This function updates the random variables of the model rand_var_
  /// According to the sufficient statistic vector stoch_sufficient_stats
  /// See documentation for further information about the update process

  /// Update the noise variance, sigma
  ScalarType NoiseVariance = sum_obs_;
  const ScalarType * it_s1 = stoch_sufficient_stats[0].memptr();
  const ScalarType * it_s2 = stoch_sufficient_stats[1].memptr();
  for(size_t i = 0; i < stoch_sufficient_stats[0].size(); ++i) {
    NoiseVariance += -2 * it_s1[i] + it_s2[i];
  }

  NoiseVariance /= obs_tot_num_ * manifold_dim_;
  noise_->SetVariance(NoiseVariance);

  /// Update Ksi : Mean & Variance
  ScalarType ksi_mean = 0.0, ksi_variance = 0.0;
  const ScalarType * it_s3 = stoch_sufficient_stats[2].memptr();
  const ScalarType * it_s4 = stoch_sufficient_stats[3].memptr();

  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
    ksi_mean     += it_s3[i];
    ksi_variance += it_s4[i];
  }

  ksi_mean     /= subjects_tot_num_;
  ksi_variance -= subjects_tot_num_ * ksi_mean * ksi_mean;
  ksi_variance /= subjects_tot_num_;

  rand_var_.UpdateRandomVariable("Ksi", {{"Mean", ksi_mean}, { "Variance", ksi_variance}});

  /// Update Tau : Mean & Variance
  ScalarType tau_mean = 0.0, tau_variance = 0.0;
  const ScalarType * it_s5 = stoch_sufficient_stats[4].memptr();
  const ScalarType * it_s6 = stoch_sufficient_stats[5].memptr();

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
  rand_var_.UpdateRandomVariable("P", {{"Mean", stoch_sufficient_stats[6](0)}});
}


Observations UnivariateModel::SimulateData(io::SimulatedDataSettings &data_settings, bool need_data)
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
  /// PITFALL2 : Take a closer look at / merge with InitializeFakeRandomVariables


  /// Initialize the model
  manifold_dim_ = data_settings.GetDimensionOfSimulatedObservations();
  subjects_tot_num_  = data_settings.GetNumberOfSimulatedSubjects();

  if (need_data){
    InitializeFakeRandomVariables();
  }

  asso_num_real_per_rand_var_["P"] = 1;
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;

  auto reals = SimulateRealizations();

  /// Update the model
  position_ = 1.0 / (1 + exp(-reals.at("P", 0)));


  /// Simulate the data
  //Find a way to initialize depending on the test case
  std::random_device rand_device;
  std::mt19937 rand_num_gen( GV::TEST_RUN ? 1 : rand_device());

  std::uniform_int_distribution<int> uni_distrib(data_settings.GetMinimumNumberOfObservations(), data_settings.GetMaximumNumberOfObservations());
  UniformRandomVariable time_points_num(60, 95);
  GaussianRandomVariable noise(0, noise_->GetVariance());

  /// Simulate the data
  Observations obs;
  for(int i = 0; i < subjects_tot_num_; ++i)
  {
    /// Get a random number of timepoints and sort them
    VectorType time_points = time_points_num.Samples(uni_distrib(rand_num_gen));
    time_points.sort();
    subj_time_points_.push_back(time_points);

    /// Simulate the data base on the time-points
    IndividualObservations indiv_obs(time_points);
    std::vector<VectorType> cog_scores;
    for(size_t j = 0; j < time_points.size(); ++j)
    {
      VectorType parallel_curve = ComputeParallelCurve(i, j);
      cog_scores.push_back(parallel_curve + noise.Samples(manifold_dim_));
    }

    indiv_obs.AddCognitiveScores(cog_scores);
    obs.AddIndividualData(indiv_obs);
  }

  /// Initialize the observation and model attributes
  obs.InitializeGlobalAttributes();
  individual_obs_date_ = obs.GetObservations();
  sum_obs_             = obs.GetTotalSumOfCognitiveScores();
  obs_tot_num_         = obs.GetTotalNumberOfObservations();

  return obs;

}

std::vector<AbstractModel::MiniBlock> UnivariateModel::GetSamplerBlocks() const
{
  std::vector<MiniBlock> blocks;

  /// Block P
  MiniBlock block_p;
  block_p.push_back(std::make_tuple<int, std::string, int>(-1, "P", 0));
  blocks.push_back(block_p);

  /// Individual blocks
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
    MiniBlock indiv_block;
    indiv_block.push_back(std::make_tuple<int, std::string, int>(i, "Ksi",i));
    indiv_block.push_back(std::make_tuple<int, std::string, int>(i, "Tau",i));

    blocks.push_back(indiv_block);
  }

  return blocks;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Log-likelihood related method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
 

AbstractModel::VectorType UnivariateModel::ComputeLogLikelihood(const Observations &obs, const MiniBlock& block_info)
{
  /// It computes the likelihood of the model. For each subject i, it sums its likelihood, namely the distance,
  /// for each time t_ij, between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve
  int type = std::get<0>(block_info[0]);
  
  if(type == -1) {
    VectorType ok(subjects_tot_num_);
    ScalarType *ok2 = ok.memptr();
    
    for (size_t i = 0; i < subjects_tot_num_; ++i)
      ok2[i] = ComputeIndividualLogLikelihood(obs.GetSubjectObservations(i), i);

    return ok;
  } else {
    return VectorType(1, ComputeIndividualLogLikelihood(obs.GetSubjectObservations(type), type));
  }
}

//TODO(clem): adapt to use in previous func
ScalarType UnivariateModel::ComputeIndividualLogLikelihood(const IndividualObservations &obs, const int subjects_tot_num_)
{
  /// Given a particular subject i, it computes its likelihood, namely the distance, for each observation t_ij,
  /// between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve

  ScalarType log_likelihood = 0;
  auto num_time_points = obs.GetNumberOfTimePoints();


  /// For each timepoints of the particular subject
  for(size_t i = 0; i < num_time_points; ++i)
  {
    VectorType subject_cog_scores = obs.GetCognitiveScore(i);
    VectorType parallel_curve = ComputeParallelCurve(subjects_tot_num_, i);
    log_likelihood += (subject_cog_scores - parallel_curve).squared_magnitude();
  }

  log_likelihood /= - 2 * noise_->GetVariance();
  log_likelihood -= num_time_points * log(2 * noise_->GetVariance() * M_PI) / 2.0;

  return log_likelihood;
}


ScalarType UnivariateModel::GetPreviousLogLikelihood(const MiniBlock &block_info) {
  
  int type = std::get<0>(block_info[0]);
  
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

void UnivariateModel::SetPreviousLogLikelihood(VectorType &log_likelihood,
                                               const MiniBlock &block_info) {
  
  int type = std::get<0>(block_info[0]);
  
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

//TODO: what is the use of the parameter?
void UnivariateModel::DisplayOutputs(const Realizations &AR)
{
  auto block_p = rand_var_.GetRandomVariable("P")->GetParameter("Mean");
  auto Tau = rand_var_.GetRandomVariable("Tau");
  auto Ksi = rand_var_.GetRandomVariable("Ksi");


  std::cout << "noise: " << noise_->GetVariance();
  std::cout << " - P: " << block_p;
  std::cout << " - T0: " << Tau->GetParameter("Mean") << " - Var(Tau): " << Tau->GetParameter("Variance");
  std::cout << " - Ksi: " << Ksi->GetParameter("Mean") << " - Var(Ksi): " << Ksi->GetParameter("Variance") << std::endl;

}

void UnivariateModel::SaveData(unsigned int iter_num, const Realizations &reals)
{
  /// It saves the random variables / realizations / whatever model parameters
  /// Mainly needed for post processingGV::TEST_DIR + "log_univariate_file.txt"


  std::ofstream log_file;

  if(!GV::TEST_RUN) {
    log_file.open(GV::BUILD_DIR + "log_univariate_file.txt", std::ofstream::out | std::ofstream::app);
    auto g = rand_var_.GetRandomVariable("G")->GetParameter("Mean");
    auto tau = rand_var_.GetRandomVariable("Tau");
    auto ksi = rand_var_.GetRandomVariable("Ksi");

    // This part should be tuned by a xml file
    log_file << "Iteration n: " << iter_num;
    log_file << " - noise: " << noise_->GetVariance();
    log_file << " - G: " << g;
    log_file << " - T0: " << tau->GetParameter("Mean") << " - Var(Tau): "
             << tau->GetParameter("Variance");
    log_file << " - Ksi: " << ksi->GetParameter("Mean") << " - Var(Ksi): "
             << ksi->GetParameter("Variance") << std::endl;
    log_file.close();
  }
  else {
    log_file.open(GV::TEST_DIR + "log_univariate_file.txt", std::ofstream::out | std::ofstream::app);
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

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void UnivariateModel::InitializeFakeRandomVariables()
{
  /// It initialize the model with particular random variables, mainly needed to simulate data
  /// Pitfall : This is not generic for need. Both with the SimulateData function, is has to be redefined / refactored
  /// according to the needs

  /// Population variables
  noise_ = std::make_shared<GaussianRandomVariable>(0.0, 0.0001);
  rand_var_.AddRandomVariable("P", "Gaussian", {0.02, 0.00001* 0.00001});

  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
  rand_var_.AddRandomVariable("Tau", "Gaussian", {70, 0.25});

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void UnivariateModel::ComputeSubjectTimePoint(const Realizations &reals, const int subject_num)
{

  /// The model introduces a time-warp for each subject, namely a time reparametrization
  /// For a given subject i, t_ij becomes alpha_i * ( t_ij - tau_i)
  /// alpha_i is the pace of disease propagation of patient i (Faster/Slower than average)
  /// tau_i is the disease onset of patient i (Beginning of the disease before/after the average time of conversion)
  
  /// If the subject number is -1, then the function recalculates the reparametrization for all the subjects
  if(subject_num != -1) {
    double acc_factor = exp(reals.at("Ksi", subject_num));
    double time_shift = reals.at("Tau", subject_num);
    subj_time_points_[subject_num] = acc_factor * (individual_obs_date_[subject_num] - time_shift);
  }
  else {
    for(size_t i = 0; i < subjects_tot_num_; ++i) {
      double acc_factor = exp(reals.at("Ksi")(i));
      double time_shift = reals.at("Tau")(i);
      
      subj_time_points_[i] = acc_factor * (individual_obs_date_[i] - time_shift);
    }
  }


}

AbstractModel::VectorType UnivariateModel::ComputeParallelCurve(int subject_num, int obs_num)
{
  ScalarType time_point = subj_time_points_[subject_num](obs_num);

  ScalarType parallel_curve = exp( - time_point );
  parallel_curve *= (1.0 / position_ - 1);
  
  
  return VectorType(1, 1.0/parallel_curve);
}
