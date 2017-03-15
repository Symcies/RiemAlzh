#include "MeshworkModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

MeshworkModel
::MeshworkModel(io::ModelSettings &MS) 
{
  m_NbIndependentSources= MS.GetIndependentSourcesNumber();
  
  std::string KernelMatrixPath = MS.GetInvertKernelPath();
  std::string InterpolationMatrixPath = MS.GetInterpolationKernelPath();

  m_InvertKernelMatrix = io::ReadData::OpenKernel(KernelMatrixPath).transpose();
  m_InterpolationMatrix = io::ReadData::OpenKernel(InterpolationMatrixPath);
  
  m_NbControlPoints = m_InvertKernelMatrix.columns();
  manifold_dim_ = m_InterpolationMatrix.rows();
  m_Thicknesses.set_size(manifold_dim_);
  m_Deltas.set_size(manifold_dim_);
  block1_.set_size(manifold_dim_);
  
}

MeshworkModel
::~MeshworkModel() 
{
  
}

  
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
MeshworkModel
::Initialize(const Observations& Obs) 
{
  /// Data-related attributes
  subjects_tot_num_          = Obs.GetNumberOfSubjects();
  indiv_obs_date_ = Obs.GetObservations();
  indiv_time_points_         = Obs.GetObservations();
  obs_tot_num_     = Obs.GetTotalNumberOfObservations();
  sum_obs_           = Obs.GetTotalSumOfLandmarks();
  
  /// Population variables
  m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.01 );

  rand_var_.AddRandomVariable("P", "Gaussian", {0.13, 0.00005 * 0.00005});
  asso_num_real_per_rand_var_["P"] = m_NbControlPoints;
  
  for(size_t i = 1; i < m_NbControlPoints; ++i)
  {
      std::string Name = "Delta#" + std::to_string(i);
      rand_var_.AddRandomVariable(Name, "Gaussian", {0, 0.003 * 0.003});
      asso_num_real_per_rand_var_[Name] = 1;
  }
  
  for(size_t i = 0; i < m_NbIndependentSources*(manifold_dim_ - 1); ++i)
  {
      std::string Name = "Beta#" + std::to_string(i);
      rand_var_.AddRandomVariable(Name, "Gaussian", {0, 0.001 * 0.001});
      asso_num_real_per_rand_var_[Name] = 1;
  }
  
  /// Individual variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", {-3.1971, 0.000000004});
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  
  rand_var_.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
      
  for(int i = 0; i < m_NbIndependentSources; ++i)
  {
      std::string Name = "S#" + std::to_string(i);
      rand_var_.AddRandomVariable(Name, "Gaussian", {0.0, 1});
      asso_num_real_per_rand_var_[Name] =  subjects_tot_num_;
  }
  
}

ScalarType 
MeshworkModel
::InitializePropositionDistributionVariance(std::string Name) 
const 
{
  Name = Name.substr(0, Name.find_first_of("#"));
  
  if("P" == Name)
      return 0.0000001;
  if("Delta" == Name)
      return 0.0000001;
  if("Beta" == Name)
      return 0.000007*0.000007;
  if("Ksi" == Name)
      return 0.00003;
  if("Tau" == Name)
      return 0.04 * 0.04;
  if("S" == Name)
      return 0.4;
     
}


void 
MeshworkModel
::UpdateModel(const Realizations &R, int Type,
            const std::vector<std::string, std::allocator<std::string>> Names) 
{
  bool ComputeThickness = false;
  bool ComputeDelta = false;
  bool ComputeBasis = false;
  bool ComputeA = false;
  bool ComputeSpaceShift = false;
  bool ComputeBlock_ = false;
  
  bool IndividualOnly = (Type > -1);
  
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
          IndividualOnly = false;
          ComputeThickness = true;
          ComputeBasis = true;
          ComputeA = true;
          ComputeSpaceShift = true;
          ComputeBlock_ = true;
          continue;
      }
      else if(Name == "Delta")
      { 
          IndividualOnly = false;
          ComputeDelta = true;
          ComputeBasis = true;
          ComputeA = true;
          ComputeSpaceShift = true;
          ComputeBlock_ = true;
          continue;
      }
      else if(Name == "Beta")
      {
          IndividualOnly = false;
          ComputeA = true;
          ComputeSpaceShift = true;
          continue;
      }
      else if(Name == "S")
      {
          ComputeSpaceShift = true;
      }
      else if(Name == "All")
      {
          IndividualOnly = false;
          ComputeSubjectTimePoint(R, -1);
          
          ComputeThickness = true;
          ComputeDelta = true;
          ComputeBasis = true;
          ComputeA = true;
          ComputeSpaceShift = true;
          ComputeBlock_ = true;
          break;
      } 
      else
      {
          std::cerr << "The random variable name " << Name << "is unknown to the meshwork model" << std::endl; 
      }
  }
  
  if(IndividualOnly)    ComputeSubjectTimePoint(R, Type);
  if(ComputeThickness)  ComputeThicknesses(R);
  if(ComputeDelta)      ComputeDeltas(R);
  if(ComputeBasis)      ComputeOrthonormalBasis();
  if(ComputeA)          ComputeAMatrix(R);
  if(ComputeSpaceShift) ComputeSpaceShifts(R);
  if(ComputeBlock_)     ComputeBlock();
}


MeshworkModel::SufficientStatisticsVector
MeshworkModel
::GetSufficientStatistics(const Realizations &R, const Observations& Obs) 
{
  /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
  VectorType S1(obs_tot_num_), S2(obs_tot_num_);
  auto itS1 = S1.begin(), itS2 = S2.begin();
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      for(size_t j = 0; j < Obs.GetNumberOfTimePoints(i); ++j)
      {
          VectorType PC = ComputeParallelCurve(i, j);
          *itS1 = dot_product(PC, Obs.GetSubjectLandmark(i, j));
          *itS2 = PC.squared_magnitude();
          ++itS1, ++itS2;
      }
  }
  
  /// Sufficient Statistic Ksi and Ksi*Ksi
  VectorType S3 = R.at("Ksi");
  VectorType S4 = R.at("Ksi") % R.at("Ksi");
  
  /// Sufficient statistic Tau and Tau*Tau
  VectorType S5 = R.at("Tau");
  VectorType S6 = R.at("Tau") % R.at("Tau");
  
  /// Sufficient statistic beta_k
  VectorType S7((manifold_dim_-1) * m_NbIndependentSources);
  ScalarType * itS7 = S7.memptr();
  for(size_t i = 0; i < ((manifold_dim_-1) * m_NbIndependentSources); ++i)
      itS7[i] = R.at("Beta#" + std::to_string(i), 0);
  
  /// Sufficient statistic p_k and p_k*p_k
  VectorType S8 = R.at("P");
  VectorType S9 = R.at("P") % R.at("P");
  
  /// Sufficient statistic delta_k
  VectorType S10(m_NbControlPoints - 1);
  ScalarType * itS10 = S10.memptr();
  for(size_t i = 1; i < m_NbControlPoints; ++i)
      itS10[i-1] = R.at("Delta#" + std::to_string(i), 0);
  
  
  /// Return the sufficient statistic vector
  SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10};
  return S;
}

void 
MeshworkModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS) 
{
  /// Update the noise variance, sigma
  ScalarType NoiseVariance = sum_obs_;
  const ScalarType * itS1 = SS[0].memptr();
  const ScalarType * itS2 = SS[1].memptr();
  for(size_t i = 0; i < SS[0].size(); ++i)
      NoiseVariance += - 2 * itS1[i] + itS2[i];
  
  NoiseVariance /= obs_tot_num_ * manifold_dim_;
  m_Noise->SetVariance(NoiseVariance);
  
  /// Update Ksi : Mean and Variance
  ScalarType KsiMean = 0.0, KsiVariance = 0.0;
  const ScalarType * itS3 = SS[2].memptr();
  const ScalarType * itS4 = SS[3].memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i) 
  {
      KsiMean     += itS3[i];
      KsiVariance += itS4[i];
  }
      
  KsiMean     /= subjects_tot_num_;
  KsiVariance -= subjects_tot_num_ * KsiMean * KsiMean;
  KsiVariance /= subjects_tot_num_;
  
  rand_var_.UpdateRandomVariable("Ksi", {{"Mean", KsiMean}, {"Variance", KsiVariance}});
  
  
  /// Update Tau : Mean and Variance
  ScalarType TauMean = 0.0, TauVariance = 0.0;
  const ScalarType * itS5 = SS[4].memptr();
  const ScalarType * itS6 = SS[5].memptr();
  
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      TauMean     += itS5[i];
      TauVariance += itS6[i];
  }
  
  TauMean     /= subjects_tot_num_;
  TauVariance -= subjects_tot_num_ * TauMean * TauMean;
  TauVariance /= subjects_tot_num_;
  
  rand_var_.UpdateRandomVariable("Tau", {{"Mean", TauMean}, {"Variance", TauVariance}});
  
  
  /// Update Beta_k : Mean
  const ScalarType * itS7 = SS[6].memptr();
  for(size_t i = 0; i < SS[6].size(); ++i)
      rand_var_.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", itS7[i]}});
   
  /// Update P_k : Mean and Var
  ScalarType PMean = 0.0, PVariance = 0.0;
  const ScalarType * itS8 = SS[7].memptr();
  const ScalarType * itS9 = SS[8].memptr();
  
  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
      PMean     += itS8[i];
      PVariance += itS9[i];
  }
  
  PMean     /= m_NbControlPoints;
  PVariance -= m_NbControlPoints * PMean * PMean;
  PVariance /= m_NbControlPoints;
  
  rand_var_.UpdateRandomVariable("P", {{"Mean", PMean}, {"Variance", PVariance}});
  
  /// Update Delta_k : Mean
  const ScalarType * itS10 = SS[9].memptr();
  for(size_t i = 0; i < m_NbControlPoints - 1; ++i)
      rand_var_.UpdateRandomVariable("Delta#" + std::to_string(i+1), {{"Mean", itS10[i]}});
}

ScalarType 
MeshworkModel
::ComputeLogLikelihood(const Observations &Obs) 
{
  double LogLikelihood = 0;

  for(size_t i = 0; i < subjects_tot_num_; ++i) 
  {
     ScalarType N = Obs.GetNumberOfTimePoints(i);
      
      for(size_t j = 0; j < N; ++j)
      {
          auto& it = Obs.GetSubjectLandmark(i, j);
          VectorType ParallelCurve = ComputeParallelCurve(i, j);
          LogLikelihood += (it - ParallelCurve).squared_magnitude();
      }
  }
  
  LogLikelihood  /= -2*m_Noise->GetVariance();
  LogLikelihood -= obs_tot_num_*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
  
  return LogLikelihood;
}

ScalarType 
MeshworkModel
::ComputeIndividualLogLikelihood(const IndividualObservations& Obs, const int SubjectNumber) 
{
  /// Get the data
  double LogLikelihood = 0;
  auto N = Obs.GetNumberOfTimePoints();
  
#pragma omp parallel for reduction(+:LogLikelihood)   
  for(size_t i = 0; i < N; ++i)
  {
      auto& it = Obs.GetLandmark(i);
      VectorType P2 = ComputeParallelCurve(SubjectNumber, i);
      LogLikelihood += (it - P2).squared_magnitude();
  }
  
  LogLikelihood /= -2*m_Noise->GetVariance();
  LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
  
  return LogLikelihood;
}

Observations
MeshworkModel
::SimulateData(io::DataSettings &DS) 
{
  typedef std::vector< std::pair< VectorType, double> > IndividualData;
  
  subjects_tot_num_ = DS.GetNumberOfSimulatedSubjects();
  
  /// Initialize the realizations and simulate them
  asso_num_real_per_rand_var_["P"] = m_NbControlPoints;
  
  for(int i = 1; i < m_NbControlPoints; ++i)
      asso_num_real_per_rand_var_["Delta#" + std::to_string(i)] = 1;
  
  
  for(size_t i = 0; i <  m_NbIndependentSources*(manifold_dim_ - 1); ++i)
      asso_num_real_per_rand_var_["Beta#" + std::to_string(i)] = 1;
      
  asso_num_real_per_rand_var_["Ksi"] = subjects_tot_num_;
  asso_num_real_per_rand_var_["Tau"] = subjects_tot_num_;
  
  for(int i = 0; i < m_NbIndependentSources; ++i)
      asso_num_real_per_rand_var_["S#" + std::to_string(i)] = subjects_tot_num_;
  
  auto R = SimulateRealizations();
  
  /// Update the model
  ComputeDeltas(R);
  ComputeThicknesses(R);
  ComputeOrthonormalBasis();
  ComputeAMatrix(R);
  ComputeSpaceShifts(R);
  
  ComputeBlock();
  
  /// Simulate the data
  std::random_device RD;
  std::mt19937 RNG(RD());
  std::uniform_int_distribution<int> Uni(DS.GetMinimumNumberOfObservations(), DS.GetMaximumNumberOfObservations());
  UniformRandomVariable NumberOfTimePoints(60, 95);
  GaussianRandomVariable Noise(0, m_Noise->GetVariance());
  
  /// Simulate the data
  Observations Obs;
  for(int i = 0; i < subjects_tot_num_; ++i)
  { 
    /// Get a random number of timepoints and sort them
    VectorType T = NumberOfTimePoints.Samples(Uni(RNG));
    T.sort();
    indiv_time_points_.push_back(T);
    
    /// Simulate the data base on the time-points
    IndividualObservations IO(T);   
    std::vector<VectorType> Landmarks;
    for(size_t j = 0; j < T.size(); ++j)
    {
      Landmarks.push_back(ComputeParallelCurve(i, j) + Noise.Samples(manifold_dim_));
    }
    
    IO.AddLandmarks(Landmarks);
    Obs.AddIndividualData(IO);
  }
  
  /// Initialize the observation and model attributes
  Obs.InitializeGlobalAttributes();
  indiv_obs_date_ = Obs.GetObservations();
  sum_obs_           = Obs.GetTotalSumOfLandmarks();
  obs_tot_num_     = Obs.GetTotalNumberOfObservations();
  
  return Obs;
}


std::vector<AbstractModel::SamplerBlock>
MeshworkModel
::GetSamplerBlocks() 
const
{
  int PopulationType = -1;
  int NbBeta = 3;
  int NbDelta = 0;
  int NbP = 5;
  
  
  std::vector<SamplerBlock> Blocks;
  
  
  /// Insert P;
  MiniBlock P;
  
  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
      P.push_back(std::make_pair("P", i));
      //P.clear();
  }
  Blocks.push_back(std::make_pair(PopulationType, P));
  
  
  /// Insert Beta_k
  /*
  MiniBlock Beta;
  int BetaModulo = (int)m_NbIndependentSources*(manifold_dim_ - 1) /NbBeta;
  for(size_t i = 0; i < m_NbIndependentSources*(manifold_dim_ - 1); ++i)
  {
      Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
      bool HingeCondition = (i%BetaModulo == 0 && i != 0);
      bool FinalCondition = (i == m_NbIndependentSources*(manifold_dim_ - 1) - 1);
      if(FinalCondition || HingeCondition)
      {
          Blocks.push_back(std::make_pair(PopulationType, Beta));
          Beta.clear();
      }
  }
   */
  MiniBlock Beta;
  for(size_t i = 0; i < m_NbIndependentSources*(manifold_dim_ - 1); ++i)
  {
      Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));   
  }
  Blocks.push_back(std::make_pair(PopulationType, Beta));
  
  /// Insert Delta_k
  MiniBlock Delta;
  for(size_t i = 1; i < m_NbControlPoints; ++i)
  {
      Delta.push_back(std::make_pair("Delta#" + std::to_string(i), 0));
  }
  Blocks.push_back(std::make_pair(PopulationType, Delta));
  
  /// Individual variables
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      MiniBlock IndividualBlock;
      IndividualBlock.push_back(std::make_pair("Ksi", i));
      IndividualBlock.push_back(std::make_pair("Tau", i));
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
MeshworkModel
::DisplayOutputs(const Realizations &R) 
{
  auto P = rand_var_.GetRandomVariable("P");
  auto PReal = R.at("P");
  
  auto Tau = rand_var_.GetRandomVariable("Tau");
  auto Ksi = rand_var_.GetRandomVariable("Ksi");
  
  double DeltaMin = R.at("Delta#1", 0);
  double DeltaMax = DeltaMin;
  for(size_t i = 1; i < 258; ++i)
  {
      double DeltaK = R.at("Delta#" + std::to_string(i), 0);
      DeltaMax = std::max(DeltaMax, DeltaK);
      DeltaMin = std::min(DeltaMin, DeltaK);
  }
  
  
  std::cout << "Noise: " << m_Noise->GetVariance() << " - PMean: " << exp(P->GetParameter("Mean"));
  std::cout << " - PVar: " << P->GetParameter("Variance") << " - PMin: " << exp(PReal.min_value()) << " - PMax: " << exp(PReal.max_value());
  std::cout << " - T0: " << Tau->GetParameter("Mean") << " - TauVar: " << Tau->GetParameter("Variance");
  std::cout << " - V0: " << exp(Ksi->GetParameter("Mean")) << " - KsiVar: " << Ksi->GetParameter("Variance");
  std::cout << " - MinDelta: " << DeltaMin << " - MaxDelta: " << DeltaMax << std::endl;
}


void 
MeshworkModel
::SaveData(unsigned int IterationNumber, const Realizations &R) 
{
  
  std::ofstream Outputs;    
  std::string FileName = "/Users/igor.koval/Documents/Work/RiemAlzh/src/io/outputs/Meshwork/Parameters" + std::to_string(IterationNumber) + ".txt";
  Outputs.open(FileName, std::ofstream::out | std::ofstream::trunc);
  
  /// Save the final noise variance
  Outputs << m_Noise->GetVariance() << std::endl;
  
  /// Save the number of subjects, the manifold dimension, the number of sources, and, the number of control points
  Outputs << subjects_tot_num_ << ", " << manifold_dim_ << ", " << m_NbIndependentSources << ", " << m_NbControlPoints << std::endl;
  
  /// Save the delta_mean -> First one being equal to 0 as the reference
  Outputs << 0 << ", ";
  for(size_t i = 1; i < m_NbControlPoints; ++i)
  {
      Outputs << rand_var_.GetRandomVariable("Delta#" + std::to_string(i))->GetParameter("Mean");
      if(i != m_NbControlPoints - 1) { Outputs << ", "; }
  }
  Outputs << std::endl;
  
  /// Save the thicknesses
  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
      Outputs << R.at("P", i);
      if(i != m_NbControlPoints - 1) { Outputs << ", ";}
  }
  Outputs << std::endl;
  
  /// Save the tau
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      Outputs << R.at("Tau", i) ;
      if(i != subjects_tot_num_ - 1) { Outputs << ", ";}
  }
  Outputs << std::endl;
  
  /// Save the ksi
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      Outputs << R.at("Ksi", i) ;
      if(i != subjects_tot_num_ - 1) { Outputs << ", ";}
  }
  
      /// Save (S_i)
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      for(size_t j = 0; j < m_NbIndependentSources; ++j)
      {
          Outputs << R.at("S#" + std::to_string(j))(i);
          if(i != m_NbIndependentSources - 1) { Outputs << ", "; }
      }
      Outputs << std::endl;
  }
  
  /// Save (W_i)
  auto SizeW = subjects_tot_num_;
  for(size_t i = 0; i < subjects_tot_num_; ++i)
  {
      VectorType W = space_shifts_.get_column(i);
      for(auto it = W.begin(); it != W.end(); ++it)
      {
          Outputs << *it;
          if(i != SizeW - 1) { Outputs << ", "; }
      }
      Outputs << std::endl;
  }
  
}


void
MeshworkModel
::InitializeFakeRandomVariables() 
{
  /// Noise 
  m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.00001 );
  
  /// Population random variables
  rand_var_.AddRandomVariable("P", "Gaussian", {0.1, 0.001 * 0.001});
  
  for(size_t i = 0; i < m_NbControlPoints; ++i)
      rand_var_.AddRandomVariable("Delta#" + std::to_string(i), "Gaussian", {0, 0.0001 * 0.0001});
  
  for(size_t i = 0; i < m_NbIndependentSources*(manifold_dim_ - 1); ++i)
      rand_var_.AddRandomVariable("Beta#" + std::to_string(i), "Gaussian", {0, 0.0001 * 0.0001});
  
  
  /// Individual random variables
  rand_var_.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
  rand_var_.AddRandomVariable("Tau", "Gaussian", {62, 0.25});
  
  for(int i = 0; i < m_NbIndependentSources; ++i)
      rand_var_.AddRandomVariable("S#" + std::to_string(i), "Gaussian", {0.0, 0.5});

  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MeshworkModel
::ComputeSubjectTimePoint(const Realizations &R, const int SubjectNumber) 
{
  if(SubjectNumber != -1) 
  {
      double AccFactor = exp(R.at("Ksi", SubjectNumber));
      double TimeShift = R.at("Tau", SubjectNumber);
      indiv_time_points_[SubjectNumber] = AccFactor * (indiv_obs_date_[SubjectNumber] - TimeShift);
  }
  else
  {

      for(size_t i = 0; i < subjects_tot_num_; ++i) 
      {
          double AccFactor = exp(R.at("Ksi")(i));
          double TimeShift = R.at("Tau")(i);

          indiv_time_points_[i] = AccFactor * (indiv_obs_date_[i] - TimeShift);
      }
  }
}

void 
MeshworkModel
::ComputeDeltas(const Realizations &R) 
{
  VectorType Delta(m_NbControlPoints);
  Delta(0) = 0.0;
  ScalarType * d = Delta.memptr();

#pragma omp parallel for
  for(size_t i = 1; i < m_NbControlPoints; ++i)
  {
      d[i] = R.at("Delta#" + std::to_string(i), 0);
  }
  
  m_Deltas = m_InterpolationMatrix * m_InvertKernelMatrix * Delta;
}


void 
MeshworkModel
::ComputeThicknesses(const Realizations &R) 
{
  m_Thicknesses = m_InterpolationMatrix * m_InvertKernelMatrix * R.at("P").exp();
}

void 
MeshworkModel
::ComputeOrthonormalBasis() 
{
      /// Get the data
  auto V0 = rand_var_.GetRandomVariable("Ksi")->GetParameter("Mean");
  V0 = exp(V0);
  
  
  VectorType U(manifold_dim_);
  ScalarType * t = m_Thicknesses.memptr();
  ScalarType * d = m_Deltas.memptr();
  ScalarType * u = U.memptr();
  
#pragma omp simd
  for(int i = 0; i < manifold_dim_; ++i)
  {
      u[i] = V0 / (t[i] * t[i])* exp(-d[i]); 
  }
  
  

  
  /// Compute the initial pivot vector U
  double Norm = U.magnitude();
  U(0) += copysign(1, -U(0)) * Norm;
      
  // Check the vectorise_row_wise function of the ArmadilloMatrixWrapper
  MatrixType U2(U);
  double NormU2 = U.squared_magnitude();
  MatrixType FinalMatrix2 = (-2.0/NormU2) * U2*U2.transpose();
  for(size_t i = 0; i < manifold_dim_; ++i)
      FinalMatrix2(i, i ) += 1;
  

  orthog_basis_ = FinalMatrix2;
}

void 
MeshworkModel
::ComputeAMatrix(const Realizations &R) 
{
      
  MatrixType NewA(manifold_dim_, m_NbIndependentSources);
  
  for(int i = 0; i < m_NbIndependentSources; ++i)
  {
      VectorType Beta(manifold_dim_, 0.0);
      for(size_t j = 0; j < manifold_dim_ - 1; ++j)
      {
          std::string Number = std::to_string(int(j + i*(manifold_dim_ - 1)));
          Beta(j) = R.at( "Beta#" + Number, 0);
      }
      
      NewA.set_column(i, orthog_basis_ * Beta);
  }
  
  a_matrix_ = NewA;
}

void 
MeshworkModel
::ComputeSpaceShifts(const Realizations &R) 
{
  
  MatrixType SS(m_NbIndependentSources, subjects_tot_num_);
  for(int i = 0; i < m_NbIndependentSources; ++i) 
  {
      SS.set_row(i, R.at("S#" + std::to_string(i)));
  }
  space_shifts_ = a_matrix_ * SS;
}

void
MeshworkModel
::ComputeBlock() 
{
  ScalarType * t = m_Thicknesses.memptr();
  ScalarType * d = m_Deltas.memptr();
  ScalarType * b = block1_.memptr();

#pragma omp simd
  for(size_t i = 0; i < manifold_dim_; ++i)
      b[i] = 1 / (t[i] * exp(d[i]));
}


AbstractModel::VectorType
MeshworkModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
  double TimePoint = indiv_time_points_[SubjectNumber](ObservationNumber);
  
  VectorType ParallelCurve(manifold_dim_);
  
  ScalarType * p = ParallelCurve.memptr();
  ScalarType * t = m_Thicknesses.memptr();
  ScalarType * d = m_Deltas.memptr();
  ScalarType * b = block1_.memptr();
  ScalarType * w = space_shifts_.get_column(SubjectNumber).memptr();
  
  for(size_t i = 0; i < manifold_dim_; ++i)
      p[i] = t[i] * exp( w[i] / (t[i]*exp(d[i])) + d[i] - TimePoint/t[i]);
      
  
  return ParallelCurve;
}





