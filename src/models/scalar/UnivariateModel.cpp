#include "UnivariateModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

UnivariateModel
::UnivariateModel(io::ModelSettings &MS) 
{

}

UnivariateModel
::~UnivariateModel() 
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
UnivariateModel
::Initialize(const Observations &Obs) 
{
  /// Data-related attributes
  m_NumberOfSubjects          = Obs.GetNumberOfSubjects();
  m_IndividualObservationDate = Obs.GetObservations();
  m_SubjectTimePoints         = Obs.GetObservations();
  m_NbTotalOfObservations     = Obs.GetTotalNumberOfObservations();
  m_SumObservations           = Obs.GetTotalSumOfCognitiveScores();
  
  /// Population variables 
  /// (Initialization of the random variable m_Random Variable 
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.00001);
  
  m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.02, 0.0001*0.0001});
  m_RealizationsPerRandomVariable["P"] = 1; 
  
  /// Individual realizations 
  /// (Initialization of the random variable m_Random Variable 
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0.13, 0.0004});
  m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
  
  m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
  m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
}


ScalarType  
UnivariateModel
::InitializePropositionDistributionVariance(std::string Name) 
const 
{
  Name = Name.substr(0, Name.find_first_of("#"));
  
  /// It returns the variance of the proposition distribution 
  /// of a given random variable of the model
  if("P" == Name)
    return 0.00001;
  if("Ksi" == Name)
    return 0.0001;
  if("Tau" == Name)
    return 0.5;
}

void 
UnivariateModel
::UpdateModel(const Realizations &R, int Type, const std::vector<std::string> Names) 
{
  /// Given a list of names (which in fact corresponds to the variables that have potentially changed),
  /// the function updates the parameters associated to these names
  
  
  /// Possible parameters to update, depending on the name being in "vect<> Names"
  bool ComputeP;
  bool ComputeIndividual = (Type > -1);
  
  /// Parameters to update, depending on the names called
  for(auto it = Names.begin(); it != Names.end(); ++it)
  {
    std::string Name = it->substr(0, it->find_first_of("#"));
    if(Name == "None")
    {
      continue;
    }
    else if(Name == "Ksi" or Name == "Tau") 
    {
      continue;
    }
    else if(Name == "P")
    {
      ComputeIndividual = false;
      ComputeP = true;
    }
    else if(Name == "All")
    {
      ComputeSubjectTimePoint(R, -1);
      ComputeIndividual = false;
      ComputeP = true;
    } 
    else
    {
      std::cerr << "The realization does not exist in the multivariate model > update model" << std::endl;
    }
  }
  
  if(ComputeIndividual) ComputeSubjectTimePoint(R, Type);
  if(ComputeP)          m_P = 1.0 / (1 + exp(-R.at("P", 0))); 
}

AbstractModel::SufficientStatisticsVector
UnivariateModel
::GetSufficientStatistics(const Realizations& R, const Observations &Obs) 
{
  /// Computation of the suffisient statistics of the model
  /// Basically, it is a vector (named S1 to S9) of vectors
  
  /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
  VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
  auto itS1 = S1.begin(), itS2 = S2.begin();
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      for(size_t j = 0; j < Obs.GetNumberOfTimePoints(i); ++j)
      {
          VectorType PC = ComputeParallelCurve(i, j);
          *itS1 = dot_product(PC, Obs.GetSubjectCognitiveScore(i, j));
          *itS2 = PC.squared_magnitude();
          ++itS1, ++itS2;
      }
  }
  
  /// S3 <- Ksi_i   &    S4 <- Ksi_i * Ksi_i
  VectorType S3 = R.at("Ksi");
  VectorType S4 = R.at("Ksi") % R.at("Ksi");
  
  /// S5 <- Tau_i   &    S6 <- Tau_i * Tau_i
  VectorType S5 = R.at("Tau");
  VectorType S6 = R.at("Tau") % R.at("Tau");
  
  /// S7 <- G 
  VectorType S7(1, R.at("P", 0));
  
  return {S1, S2, S3, S4, S5, S6, S7};
}

void
UnivariateModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS) 
{
    /// This function updates the random variables of the model m_RandomVariables
  /// According to the sufficient statistic vector SS
  /// See documentation for further information about the update process
  
  /// Update the noise variance, sigma
  ScalarType NoiseVariance = m_SumObservations;
  const ScalarType * itS1 = SS[0].memptr();
  const ScalarType * itS2 = SS[1].memptr();
  for(size_t i = 0; i < SS[0].size(); ++i)
      NoiseVariance += - 2 * itS1[i] + itS2[i];

  NoiseVariance /= m_NbTotalOfObservations * m_ManifoldDimension;
  m_Noise->SetVariance(NoiseVariance);
  
  /// Update Ksi : Mean & Variance
  ScalarType KsiMean = 0.0, KsiVariance = 0.0;
  const ScalarType * itS3 = SS[2].memptr();
  const ScalarType * itS4 = SS[3].memptr();
  
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
    KsiMean     += itS3[i];
    KsiVariance += itS4[i]; 
  }
  
  KsiMean     /= m_NumberOfSubjects;
  KsiVariance -= m_NumberOfSubjects * KsiMean * KsiMean;
  KsiVariance /= m_NumberOfSubjects;
  
  m_RandomVariables.UpdateRandomVariable("Ksi", {{"Mean", KsiMean}, { "Variance", KsiVariance}});
  
  /// Update Tau : Mean & Variance
  ScalarType TauMean = 0.0, TauVariance = 0.0;
  const ScalarType * itS5 = SS[4].memptr();
  const ScalarType * itS6 = SS[5].memptr();
  
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
    TauMean     += itS5[i];
    TauVariance += itS6[i];
  }
  
  TauMean     /= m_NumberOfSubjects;
  TauVariance -= m_NumberOfSubjects * TauMean * TauMean;
  TauVariance /= m_NumberOfSubjects;
  
  m_RandomVariables.UpdateRandomVariable("Tau", {{"Mean", TauMean}, {"Variance", TauVariance}});
  
  /// Update G : Mean
  m_RandomVariables.UpdateRandomVariable("P", {{"Mean", SS[6](0)}});
}

ScalarType 
UnivariateModel
::ComputeLogLikelihood(const Observations &Obs) 
{
  /// It computes the likelihood of the model. For each subject i, it sums its likelihood, namely the distance,
  /// for each time t_ij, between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve
  
  ScalarType LogLikelihood = 0;

  /// For each subject
  for(size_t i = 0; i < m_NumberOfSubjects; ++i) 
  {
    ScalarType N = Obs.GetNumberOfTimePoints(i);
    
    /// For each timepoint
    for(size_t j = 0; j < N; ++j)
    {
      VectorType it = Obs.GetSubjectCognitiveScore(i, j);
      VectorType PC = ComputeParallelCurve(i, j);
      LogLikelihood += (it - PC).squared_magnitude();
    }
  }
  
  LogLikelihood /= -2*m_Noise->GetVariance();
  LogLikelihood -= m_NbTotalOfObservations * log(sqrt(2 * m_Noise->GetVariance() * M_PI));
  
  return LogLikelihood;
}

ScalarType 
UnivariateModel
::ComputeIndividualLogLikelihood(const IndividualObservations &Obs, const int SubjectNumber) 
{
  /// Given a particular subject i, it computes its likelihood, namely the distance, for each observation t_ij,
  /// between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve
  
  ScalarType LogLikelihood = 0;
  auto N = Obs.GetNumberOfTimePoints();
  
  /// For each timepoints of the particular subject
  for(size_t i = 0; i < N; ++i)
  {
    VectorType it = Obs.GetCognitiveScore(i);
    VectorType PC = ComputeParallelCurve(SubjectNumber, i);
    LogLikelihood += (it - PC).squared_magnitude();
  }
  
  LogLikelihood /= - 2 * m_Noise->GetVariance();
  LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
  
  return LogLikelihood;
}

Observations
UnivariateModel
::SimulateData(io::DataSettings &DS) 
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
  m_ManifoldDimension = DS.GetCognitiveScoresDimension();
  m_NumberOfSubjects  = DS.GetNumberOfSimulatedSubjects();
  
  m_RealizationsPerRandomVariable["P"] = 1;
  m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
  m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
  
  auto R = SimulateRealizations();
  
  /// Update the model
  m_P = 1.0 / (1 + exp(-R.at("P", 0)));
  
      
  /// Simulate the data
  std::random_device RD;
  std::mt19937 RNG(RD());
  std::uniform_int_distribution<int> Uni(DS.GetMinimumNumberOfObservations(), DS.GetMaximumNumberOfObservations());
  UniformRandomVariable NumberOfTimePoints(60, 95);
  GaussianRandomVariable Noise(0, m_Noise->GetVariance());
  
  /// Simulate the data
  Observations Obs;
  for(int i = 0; i < m_NumberOfSubjects; ++i)
  { 
    /// Get a random number of timepoints and sort them
    VectorType T = NumberOfTimePoints.Samples(Uni(RNG));
    T.sort();
    m_SubjectTimePoints.push_back(T);
    
    /// Simulate the data base on the time-points
    IndividualObservations IO(T);   
    std::vector<VectorType> Landmarks;
    for(size_t j = 0; j < T.size(); ++j)
    {
      VectorType PC = ComputeParallelCurve(i, j);
      Landmarks.push_back(PC + Noise.Samples(m_ManifoldDimension));
    }
    
    IO.AddLandmarks(Landmarks);
    Obs.AddIndividualData(IO);
  }
  
  /// Initialize the observation and model attributes
  Obs.InitializeGlobalAttributes();
  m_IndividualObservationDate = Obs.GetObservations();
  m_SumObservations           = Obs.GetTotalSumOfLandmarks();
  m_NbTotalOfObservations     = Obs.GetTotalNumberOfObservations();
  
  return Obs;
  
}

std::vector<AbstractModel::SamplerBlock>
UnivariateModel
::GetSamplerBlocks() 
const 
{
  std::vector<SamplerBlock> Blocks;
  
  /// Block P
  MiniBlock P;
  P.push_back(std::make_pair("P", 0));
  Blocks.push_back(std::make_pair(-1, P));
  
  /// Individual blocks
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
    MiniBlock IndividualBlock;
    IndividualBlock.push_back(std::make_pair("Ksi",i));
    IndividualBlock.push_back(std::make_pair("Tau",i));
    
    Blocks.push_back(std::make_pair(i, IndividualBlock));
  }

  return Blocks;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void
UnivariateModel
::DisplayOutputs(const Realizations &AR) 
{
  auto P = m_RandomVariables.GetRandomVariable("P")->GetParameter("Mean");
  auto Tau = m_RandomVariables.GetRandomVariable("Tau");
  auto Ksi = m_RandomVariables.GetRandomVariable("Ksi");
  
  
  std::cout << "Noise: " << m_Noise->GetVariance();
  std::cout << " - P: " << P;
  std::cout << " - T0: " << Tau->GetParameter("Mean") << " - Var(Tau): " << Tau->GetParameter("Variance");
  std::cout << " - Ksi: " << Ksi->GetParameter("Mean") << " - Var(Ksi): " << Ksi->GetParameter("Variance") << std::endl;

}

void 
UnivariateModel
::SaveData(unsigned int IterationNumber, const Realizations &R) 
{
  /// It saves the random variables / realizations / whatever model parameters
  /// Mainly needed for post processing
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
UnivariateModel
::InitializeFakeRandomVariables() 
{
  /// It initialize the model with particular random variables, mainly needed to simulate data
  /// Pitfall : This is not generic for need. Both with the SimulateData function, is has to be redefined / refactored
  /// according to the needs
  
  /// Population variables
  m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.0001);
  m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.02, 0.00001* 0.00001});

  /// Individual variables
  m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
  m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {70, 0.25});

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
UnivariateModel
::ComputeSubjectTimePoint(const Realizations &R, const int SubjectNumber) 
{  

  /// The model introduces a time-warp for each subject, namely a time reparametrization
  /// For a given subject i, t_ij becomes alpha_i * ( t_ij - tau_i)
  /// alpha_i is the pace of disease propagation of patient i (Faster/Slower than average)
  /// tau_i is the disease onset of patient i (Beginning of the disease before/after the average time of conversion)
  
  /// If the subject number is -1, then the function recalculates the reparametrization for all the subjects
  if(SubjectNumber != -1) 
  {
      double AccFactor = exp(R.at("Ksi", SubjectNumber));
      double TimeShift = R.at("Tau", SubjectNumber);
      m_SubjectTimePoints[SubjectNumber] = AccFactor * (m_IndividualObservationDate[SubjectNumber] - TimeShift);
  }
  else
  {

      for(size_t i = 0; i < m_NumberOfSubjects; ++i) 
      {
          double AccFactor = exp(R.at("Ksi")(i));
          double TimeShift = R.at("Tau")(i);

          m_SubjectTimePoints[i] = AccFactor * (m_IndividualObservationDate[i] - TimeShift);
      }
  }
  
  
}

AbstractModel::VectorType 
UnivariateModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
  ScalarType TimePoint = m_SubjectTimePoints[SubjectNumber](ObservationNumber);
  
  ScalarType ParallelCurve = exp( - TimePoint / (m_P * (1 - m_P)));
  ParallelCurve *= (1.0 / m_P - 1);
  
  return VectorType(1, 1.0/ParallelCurve);
}
