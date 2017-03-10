#include "MultivariateModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


MultivariateModel
::MultivariateModel(io::ModelSettings &MS) 
{
  /// Initialize the data dimension and the number of sources
  m_NbIndependentSources = MS.GetIndependentSourcesNumber();
}

MultivariateModel
::~MultivariateModel() 
{
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultivariateModel
::Initialize(const Observations& Obs) 
{
  /// The function initialize different attributes of the model
  /// As well as the specific random variables and their realizations used by the model
  // TODO : Pitfall : How to have a smart initialization, that may come out of a .txt / .csv file
  // TODO : instead of using the same initialization or modifiying it in the code
  // TODO : Some cases may fall in between the two cases (default values needed)
  
  
  /// Data-related attributes
  m_ManifoldDimension         = Obs.GetSubjectObservations(0).GetCognitiveScore(0).size();
  m_NumberOfSubjects          = Obs.GetNumberOfSubjects();
  m_IndividualObservationDate = Obs.GetObservations();
  m_SubjectTimePoints         = Obs.GetObservations();
  m_NbTotalOfObservations     = Obs.GetTotalNumberOfObservations();
  m_SumObservations           = Obs.GetTotalSumOfCognitiveScores();
  
  /// Initialize the size of some parameters
  m_Deltas.set_size(m_ManifoldDimension);
  m_Block.set_size(m_ManifoldDimension);
  
  /// Population variables 
  /// (Initialization of the random variable m_Random Variable 
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.00001);
  
  m_RandomVariables.AddRandomVariable("G", "Gaussian", {0.12, 0.0001 * 0.0001});
  m_RealizationsPerRandomVariable["G"] = 1;
  
  for(size_t i = 1; i < m_ManifoldDimension; ++i)
  {
    std::string NameDelta = "Delta#" + std::to_string(i);
    m_RandomVariables.AddRandomVariable(NameDelta, "Gaussian", {0, 0.003 * 0.003});
    m_RealizationsPerRandomVariable[NameDelta] = 1;
  }
  
  for(int i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
  {
    std::string Name = "Beta#" + std::to_string(i);
    m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.001*0.001});
    m_RealizationsPerRandomVariable[Name] = 1;
  }
  
  /// Individual realizations 
  /// (Initialization of the random variable m_Random Variable 
  /// and the associated number of realizations m_RealizationPerRandomVariable)
  m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0.13, 0.0004});
  m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
  
  m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
  m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
  
  for(int i = 0; i < m_NbIndependentSources; ++i)
  {
    std::string Name = "S#" + std::to_string(i);
    m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0.0, 1});
    m_RealizationsPerRandomVariable[Name] = m_NumberOfSubjects;
  }
}

ScalarType
MultivariateModel
::InitializePropositionDistributionVariance(std::string Name) 
const 
{
  Name = Name.substr(0, Name.find_first_of("#"));
  
  /// It returns the variance of the proposition distribution 
  /// of a given random variable of the model
  if("G" == Name)
    return 0.00001;
  if("Delta" == Name)
    return 0.000001;
  if("Beta" == Name)
    return 0.00001;
  if("Ksi" == Name)
    return 0.0001;
  if("Tau" == Name)
    return 0.5;
  if("S" == Name)
    return 0.7;
}

void
MultivariateModel
::UpdateModel(const Realizations &R, int Type,
            const std::vector<std::string, std::allocator<std::string>> Names)
{
  /// Given a list of names (which in fact corresponds to the variables that have potentially changed),
  /// the function updates the parameters associated to these names
  
  
  /// Possible parameters to update, depending on the name being in "vect<> Names"
  bool ComputeG = false;
  bool ComputeDelta = false;
  bool ComputeBasis = false;
  bool ComputeA = false;
  bool ComputeSpaceShift = false;
  bool ComputeBlock_ = false;
  bool IndividualOnly = (Type > -1);
  
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
    else if("G" == Name)
    {
      IndividualOnly = false;
      ComputeBasis = true;
      ComputeA = true;
      ComputeSpaceShift = true;
      ComputeBlock_ = true;
    }
    else if("Delta" == Name)
    {
      IndividualOnly = false;
      ComputeBasis = true;
      ComputeA = true;
      ComputeSpaceShift = true;
      ComputeBlock_ = true;
    }
    else if("Beta" == Name) 
    {
      IndividualOnly = false;
      ComputeA = true;
      ComputeSpaceShift = true;
    }
    else if("S" == Name)
    {
      ComputeSpaceShift = true;
    }
    else if("All" == Name)
    {
      ComputeSubjectTimePoint(R, -1);
      IndividualOnly = false;
      ComputeG = true;
      ComputeDelta = true;
      ComputeBasis = true;
      ComputeA = true;
      ComputeSpaceShift = true;
      ComputeBlock_ = true;
    } 
    else
    {
      std::cerr << "The realization does not exist in the multivariate model > update model" << std::endl;
    }
  }
  
  
  // TODO : To parse it even faster, update just the coordinates within the names
  if(IndividualOnly) ComputeSubjectTimePoint(R, Type);
  
  if(ComputeG)          m_G = exp(R.at("G", 0));
  if(ComputeDelta)      ComputeDeltas(R);
  if(ComputeBasis)      ComputeOrthonormalBasis();
  if(ComputeA)          ComputeAMatrix(R);
  if(ComputeSpaceShift) ComputeSpaceShifts(R);
  if(ComputeBlock_)     ComputeBlock(R);
}

AbstractModel::SufficientStatisticsVector
MultivariateModel
::GetSufficientStatistics(const Realizations &R, const Observations& Obs) 
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
  VectorType S7(1, R.at("G", 0));
  
  /// S8 <- beta_k
  VectorType S8((m_ManifoldDimension-1) * m_NbIndependentSources);
  ScalarType * itS8 = S8.memptr();
  for(size_t i = 0; i < S8.size(); ++i)
      itS8[i] = R.at("Beta#" + std::to_string(i), 0);
  
  /// S8 <- delta_k
  VectorType S9(m_ManifoldDimension - 1);
  ScalarType * itS9 = S9.memptr();
  for(size_t i = 1; i < S9.size(); ++i)
      itS9[i] = R.at("Delta#" + std::to_string(i), 0);

  return {S1, S2, S3, S4, S5, S6, S7, S8, S9};
  
}


void 
MultivariateModel
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
  m_RandomVariables.UpdateRandomVariable("G", {{"Mean", SS[6](0)}});
  
  /// Update Beta_k : Mean
  const ScalarType * itS8 = SS[7].memptr();
  for(size_t i = 0; i < SS[7].size(); ++i)
    m_RandomVariables.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", itS8[i]}});
  
  /// Update Delta_k : Mean
  const ScalarType * itS9 = SS[8].memptr();
  for(size_t i = 1; i < SS[8].size(); ++i)
    m_RandomVariables.UpdateRandomVariable("Delta#" + std::to_string(i), {{"Mean", itS9[i]}});
  
}


double
MultivariateModel
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
      auto& it = Obs.GetSubjectCognitiveScore(i, j);
      VectorType PC = ComputeParallelCurve(i, j);
      LogLikelihood += (it - PC).squared_magnitude();
    }
  }
  
  LogLikelihood /= -2*m_Noise->GetVariance();
  LogLikelihood -= m_NbTotalOfObservations * log(sqrt(2 * m_Noise->GetVariance() * M_PI));
  
  return LogLikelihood;
}


double
MultivariateModel
::ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber) 
{
  /// Given a particular subject i, it computes its likelihood, namely the distance, for each observation t_ij,
  /// between the observation y_ij and the prediction f(t_ij) = ComputeParallelCurve
  
  ScalarType LogLikelihood = 0;
  auto N = Obs.GetNumberOfTimePoints();
  
  /// For each timepoints of the particular subject
  for(size_t i = 0; i < N; ++i)
  {
    auto& it = Obs.GetCognitiveScore(i);
    VectorType PC = ComputeParallelCurve(SubjectNumber, i);
    LogLikelihood += (it - PC).squared_magnitude();
  }
  
  LogLikelihood /= - 2 * m_Noise->GetVariance();
  LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
  
  return LogLikelihood;
}

Observations
MultivariateModel
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
  
  m_RealizationsPerRandomVariable["G"] = 1;
  
  for(size_t i = 1; i < m_ManifoldDimension; ++i)
    m_RealizationsPerRandomVariable["Delta#" + std::to_string(i)] = 1;

  for(size_t i = 0; i <  m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    m_RealizationsPerRandomVariable["Beta#" + std::to_string(i)] = 1;

  m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
  m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;

  for(int i = 0; i < m_NbIndependentSources; ++i)
    m_RealizationsPerRandomVariable["S#" + std::to_string(i)] = m_NumberOfSubjects;

  auto R = SimulateRealizations();
  
  /// Update the model
  m_G = exp(R.at("G", 0));
  ComputeDeltas(R);
  ComputeOrthonormalBasis();
  ComputeAMatrix(R);
  ComputeSpaceShifts(R);
  ComputeBlock(R);
  
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
      Landmarks.push_back(ComputeParallelCurve(i, j) + Noise.Samples(m_ManifoldDimension));
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
MultivariateModel
::GetSamplerBlocks() 
const 
{
  /// It defines the blocks used in the sampler class. A block is defined by a type, and a vector of pairs
  /// Each pair is composed of <Name of the random variable, Realization number>
  /// TO IMPROVE : These blocks may change during the iterations or they can be random, ...
  
  int PopulationType = -1;
  std::vector<SamplerBlock> Blocks;
  
  /// Block G
  MiniBlock G;
  G.push_back(std::make_pair("G", 0));
  Blocks.push_back(std::make_pair(PopulationType, G));
  
  /// Block Delta
  MiniBlock Delta;
  for(size_t i = 1; i < m_ManifoldDimension; ++i)
    Delta.push_back(std::make_pair("Delta#" + std::to_string(i), 0));
  Blocks.push_back(std::make_pair(PopulationType, Delta));
  
  /// Block Beta
  MiniBlock Beta;
  for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
  Blocks.push_back(std::make_pair(PopulationType, Beta));
  
  /// Individual variables --> Each block corresponds to one individual with all his/her associated random variables
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
    MiniBlock IndividualBlock;
    IndividualBlock.push_back(std::make_pair("Ksi",i));
    IndividualBlock.push_back(std::make_pair("Tau",i));
    for(size_t j = 0; j < m_NbIndependentSources; ++j)
      IndividualBlock.push_back(std::make_pair("S#" + std::to_string(j), i));
    
    Blocks.push_back(std::make_pair(i, IndividualBlock));
  }
  
  return Blocks;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultivariateModel
::DisplayOutputs(const Realizations &AR) 
{
  /// It defines the outputs displayed on the terminal
  
  auto G = m_RandomVariables.GetRandomVariable("G")->GetParameter("Mean");
  auto Tau = m_RandomVariables.GetRandomVariable("Tau");
  auto Ksi = m_RandomVariables.GetRandomVariable("Ksi");
  
  
  std::cout << "Noise: " << m_Noise->GetVariance();
  std::cout << " - G: " << G;
  std::cout << " - T0: " << Tau->GetParameter("Mean") << " - Var(Tau): " << Tau->GetParameter("Variance");
  std::cout << " - Ksi: " << Ksi->GetParameter("Mean") << " - Var(Ksi): " << Ksi->GetParameter("Variance") << std::endl;
}

void 
MultivariateModel
::SaveData(unsigned int IterationNumber, const Realizations &R) 
{
  /// It saves the random variables / realizations / whatever model parameters
  /// Mainly needed for post processing
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////


void
MultivariateModel
::InitializeFakeRandomVariables() 
{
  /// It initialize the model with particular random variables, mainly needed to simulate data
  /// Pitfall : This is not generic for need. Both with the SimulateData function, is has to be redefined / refactored
  /// according to the needs
  
  /// Population variables
  m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.0001);
  
  m_RandomVariables.AddRandomVariable("G", "Gaussian", {0.08, 0.00001* 0.00001});
  
  for(size_t i = 0; i < m_ManifoldDimension; ++i)
    m_RandomVariables.AddRandomVariable("Delta#" + std::to_string(i), "Gaussian", {0.036, 0.004 * 0.004});
  
  for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    m_RandomVariables.AddRandomVariable("Beta#" + std::to_string(i), "Gaussian", {0, 0.0001 * 0.0001});
  
  /// Individual variables
  m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
  m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {70, 0.25});
  for(size_t i = 0; i < m_NbIndependentSources; ++i)
    m_RandomVariables.AddRandomVariable("S#" + std::to_string(i), "Gaussian", {0.0, 0.5});
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultivariateModel
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

void 
MultivariateModel
::ComputeDeltas(const Realizations &R) 
{
  /// Create the vector of time shifts delta accross all coordinates
  /// delta = (delta_1, delta_2, ..., delta_N) where delta_1 = 0
  
  m_Deltas(0) = 0.0;
  ScalarType * d = m_Deltas.memptr();
  
  for(size_t i = 1; i < m_ManifoldDimension; ++i)
    d[i] = R.at("Delta#" + std::to_string(i), 0);
}

void 
MultivariateModel
::ComputeOrthonormalBasis() 
{
  /// It computes a basis of vector orthogonal to the geodesic gamma_derivative at t0, with respect to
  /// the scalar product defined on the riemannian manifold
  /// Further information about the mathematical operation in the documentation
  
  ScalarType V0 = - m_G / (m_G + 1) * exp(m_RandomVariables.GetRandomVariable("Ksi")->GetParameter("Mean"));
  VectorType U(m_ManifoldDimension);
  ScalarType * u = U.memptr();
  ScalarType * d = m_Deltas.memptr();
  
  for(size_t i = 0; i < m_ManifoldDimension; ++i)
    u[i] = V0 * (1.0/m_G + exp(d[i]));
  
  /// Compute the initial pivot vector U
  double Norm = U.magnitude();
  U(0) += copysign(1, -U(0)) * Norm;
      
  // Check the vectorise_row_wise function of the ArmadilloMatrixWrapper
  MatrixType U2(U);
  double NormU2 = U.squared_magnitude();
  MatrixType FinalMatrix2 = (-2.0/NormU2) * U2*U2.transpose();
  for(size_t i = 0; i < m_ManifoldDimension; ++i)
    FinalMatrix2(i, i ) += 1;
  
  m_OrthogonalBasis = FinalMatrix2;
}

void 
MultivariateModel
::ComputeAMatrix(const Realizations &R) 
{
  /// It computes the A matrix (ref documentation) based on the orthonormal basis B and the beta coefficients
  
  MatrixType NewA(m_ManifoldDimension, m_NbIndependentSources);
  
  for(int i = 0; i < m_NbIndependentSources; ++i)
  {
    VectorType Beta(m_ManifoldDimension, 0.0);
    for(size_t j = 0; j < m_ManifoldDimension - 1; ++j)
    {
      std::string Number = std::to_string(int(j + i*(m_ManifoldDimension - 1)));
      Beta(j) = R.at( "Beta#" + Number, 0);
    }
    
    NewA.set_column(i, m_OrthogonalBasis * Beta);
  }
  
  m_AMatrix = NewA;
}

void 
MultivariateModel
::ComputeSpaceShifts(const Realizations &R) 
{
  /// It computes the individual space shifts w_i (ref documentation or NIPS paper)
  MatrixType SS(m_NbIndependentSources, m_NumberOfSubjects);
  
  for(int i = 0; i < m_NbIndependentSources; ++i) 
    SS.set_row(i, R.at("S#" + std::to_string(i)));
  
  m_SpaceShifts = m_AMatrix * SS;
}

void 
MultivariateModel
::ComputeBlock(const Realizations &R) 
{
  /// It compute a block used to optimise the number of computation within the likelihood function
  ScalarType * b = m_Block.memptr();
  ScalarType * d = m_Deltas.memptr();
  
  for(size_t i = 0; i < m_ManifoldDimension; ++i) 
  {
    ScalarType Q = m_G * exp(-d[i]);
    b[i] = (Q + 1) * (Q + 1) / Q;
  }
}

AbstractModel::VectorType
MultivariateModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
  /// Given a subject i = SubjectNumber and a timepoint t_ij where j = Observation number,
  /// it computes the f(t_i) value corresponding to the current model
  
  VectorType ParallelCurve(m_ManifoldDimension);
  ScalarType * p = ParallelCurve.memptr();
  
  double Time = m_SubjectTimePoints[SubjectNumber](ObservationNumber);
  ScalarType * d = m_Deltas.memptr();
  ScalarType * b = m_Block.memptr();
  ScalarType * w = m_SpaceShifts.get_column(SubjectNumber).memptr();
  
  for(size_t i = 0; i < m_ManifoldDimension; ++i)
    p[i] = 1.0/(1.0 + m_G * exp(-d[i] - w[i]*b[i] - Time));
  
  return ParallelCurve;
}