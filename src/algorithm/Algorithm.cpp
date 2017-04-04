#include "Algorithm.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


Algorithm::Algorithm(io::AlgorithmSettings& settings) {
  /// Initialize the algorithm attributes
  // TODO : check if it is enough, based on future needs
  max_iter_num_     = settings.GetMaximumNumberOfIterations();
  burnin_iter_num_  = settings.GetNumberOfBurnInIterations();
  output_iter_      = settings.GetOutputDisplayIteration();
  data_save_iter_   = settings.GetDataSaveIteration();
}

Algorithm::Algorithm(io::AlgorithmSettings& settings, std::shared_ptr<AbstractModel> model, std::shared_ptr<AbstractSampler> sampler) {
  /// Initialize the algorithm attributes
  // TODO : check if it is enough, based on future needs
  max_iter_num_     = settings.GetMaximumNumberOfIterations();
  burnin_iter_num_  = settings.GetNumberOfBurnInIterations();
  output_iter_      = settings.GetOutputDisplayIteration();
  data_save_iter_   = settings.GetDataSaveIteration();

  SetModel(model);
  SetSampler(sampler);
}


Algorithm::~Algorithm() {
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Core method :
////////////////////////////////////////////////////////////////////////////////////////////////////

void Algorithm::ComputeMCMCSAEM(const Observations& obs) {
  /// This function is core to the software. It initialize parts of the model and sampler
  /// and runs the MCMC-SAEM algorithm. The class attributes define the properties of the MCMC-SAEM
  InitializeModel(obs);
  InitializeSampler();
  InitializeStochasticSufficientStatistics(obs);

  for(int iter = 0; iter < max_iter_num_; iter ++)
  {
    IterationMCMCSAEM(obs, iter);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) : Initialization
////////////////////////////////////////////////////////////////////////////////////////////////////

void Algorithm::InitializeStochasticSufficientStatistics(const Observations& obs)
{
  /// It initialize the stochastic sufficient statistics by copying the one from the model.
  /// Pitfall : it computes the suff stat of the model where only the length is needed

  stochastic_sufficient_stats_ = model_->GetSufficientStatistics(*realizations_, obs);

  for(auto&& it : stochastic_sufficient_stats_){
    std::fill(it.begin(), it.end(), 0.0);
  }

}


void Algorithm::InitializeModel(const Observations& obs)
{
  /// It initialize the model, draw its respective realizations and initialize the acceptance ratios
  /// which are key to observe the algorithm convergence

  model_->Initialize(obs);
  Realizations real = model_->SimulateRealizations();

  realizations_ = std::make_shared<Realizations>(real);

  auto init = {std::make_tuple<int, std::string, int>(-1, "All", 0)};
  model_->UpdateModel(real, init);

  for(auto it = realizations_->begin(); it != realizations_->end(); ++it)
  {

    VectorType v(it->second.size(), 0);
    acceptance_ratio_[it->first] = v;
  }

}

void Algorithm::InitializeSampler()
{
  sampler_->InitializeSampler(*realizations_, *model_);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) : Iteration
////////////////////////////////////////////////////////////////////////////////////////////////////


void Algorithm::IterationMCMCSAEM(const Observations& obs, int iter){
  if( IsOutputIteration(iter) ) {
    DisplayIterations(iter);
  }

  ComputeSimulationStep(obs, iter);
  SufficientStatisticsVector sufficient_stats = model_->GetSufficientStatistics(*realizations_, obs);
  ComputeStochasticApproximation(sufficient_stats, iter);
  model_->UpdateRandomVariables(stochastic_sufficient_stats_);

  if( IsOutputIteration(iter) ) {
    DisplayOutputs();
  }
  if( IsDataSaveIteration(iter)) {
    model_->SaveData(iter, *realizations_);
  }
}

void Algorithm::ComputeSimulationStep(const Observations& obs, int iter)
{
  /// It compute the simulate step to draw new realizations based on the previous one.
  /// The previous realizations are kept to compute the acceptance ratio

  Realizations prev_reals = *realizations_;
  sampler_->Sample(*realizations_, *model_, obs);
  ComputeAcceptanceRatio(prev_reals , iter);
}


void Algorithm::ComputeAcceptanceRatio(Realizations& prev_real_, int iter)
{
  for(auto it = realizations_->begin(); it != acceptance_ratio_.end(); ++it)
  {
      int key_var = it->first;

      VectorType new_real = it->second;
      VectorType prev_real = prev_real_.at(key_var);

      auto iter_prev_real_ = prev_real.begin();
      auto iter_new_real = new_real.begin();
      auto iter_accept_ratio = acceptance_ratio_.at(key_var).begin();

      for(    ; iter_prev_real_ != prev_real.end() && iter_new_real != new_real.end() && iter_accept_ratio != acceptance_ratio_.at(key_var).end()
              ; ++iter_prev_real_, ++iter_new_real, ++iter_accept_ratio)
      {
          bool Change = (*iter_new_real != *iter_prev_real_);
          *iter_accept_ratio = (*iter_accept_ratio * iter + Change ) / (iter + 1);
      }
  }
  if(IsOutputIteration(iter)) {
    DisplayAcceptanceRatio();
  }
}

void Algorithm::ComputeStochasticApproximation(SufficientStatisticsVector& stat_vector, int iter)
{
  /// It comptue the stochastic approximation step S_(k+1) = stat_vector(k) + stochastic variation of the previous state
  assert(stat_vector.size() == stochastic_sufficient_stats_.size());

  double step_size = DecreasingStepSize(iter);
  auto it_stoch_s = stochastic_sufficient_stats_.begin();

  for(auto it_s = stat_vector.begin(); it_s != stat_vector.end(); ++it_s, ++it_stoch_s)
      *it_stoch_s += step_size * (*it_s - *it_stoch_s);

}


double Algorithm::DecreasingStepSize(int iter)
{
    double Q = (double)iter - (double)burnin_iter_num_;
    double epsilon = std::max(1.0, Q);
    return 1.0 / pow(epsilon, 0.6); // TODO : TO CHECK
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) : Conditions
////////////////////////////////////////////////////////////////////////////////////////////////////

bool Algorithm::IsOutputIteration(int iter){
  if(output_iter_ == 0) {
    return false;
  }
  return (iter%output_iter_ == 0);
}

bool Algorithm::IsDataSaveIteration(int iter){
  if(data_save_iter_ == 0) {
    return false;
  }
  return (iter%data_save_iter_ == 0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s)
////////////////////////////////////////////////////////////////////////////////////////////////////

void Algorithm::DisplayIterations(int iter){
    std::cout  << std::endl << "--------------------- Iteration " <<
    iter << " -------------------------------" << std::endl;
}

void Algorithm::DisplayOutputs()
{
    model_->DisplayOutputs(*realizations_);
}

void Algorithm::DisplayAcceptanceRatio() {
    std::cout << "AcceptRatio: ";

    auto names_to_show = {"Tau", "Ksi", "Beta#1", "Delta#3"};

    for(auto it = names_to_show.begin(); it != names_to_show.end(); ++it)
    {
        std::string name = *it;
        int key = realizations_->ReverseNameToKey(name);
        VectorType ratios = acceptance_ratio_.at(key);

        std::cout << name << ": " << ratios.mean_value();
        if(ratios.size() != 1)
            std::cout << " & Min: " << ratios.min_value() << " & Max: " << ratios.max_value();
        std::cout << ". ";
    }
    std::cout << std::endl;

    /// Useless for now because all delta or beta are the same
    /*
    auto InName = {"Delta", "Beta"};
    for(auto it = InName.begin(); it != InName.end(); ++it)
    {
        double Min = 1;
        double Max = 0;
        double Mean = 0;
        int Count = 0;
        for(auto it2 = acceptance_ratio_.begin(); it2 != acceptance_ratio_.end(); ++it2)
        {
            std::string name = realizations_->ReverseKeyToName(it2->first);
            name = name.substr(0, name.find_first_of("#"));
            if(name == *it)
            {
                ++Count;
                double AccepVal = it2->second(0);
                Min = std::min(Min, AccepVal);
                Max = std::max(Max, AccepVal);
                Mean += AccepVal;
            }
        }

        std::cout << *it <<"(Min/Mean/Max): " << Min << "/" << Mean/Count << "/" << Max << ". ";
    }
    std::cout << std::endl;
    */

}
