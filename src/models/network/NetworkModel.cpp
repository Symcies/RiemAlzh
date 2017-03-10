#include "NetworkModel.h"


NetworkModel
::NetworkModel(io::ModelSettings &MS)
{
  m_NbIndependentSources= MS.GetIndependentSourcesNumber();

  std::string KernelMatrixPath = MS.GetInvertKernelPath();
  std::string InterpolationMatrixPath = MS.GetInterpolationKernelPath();

  m_InvertKernelMatrix = io::ReadData::OpenKernel(KernelMatrixPath).transpose();
  m_InterpolationMatrix = io::ReadData::OpenKernel(InterpolationMatrixPath);

  m_NbControlPoints = m_InvertKernelMatrix.columns();
  m_ManifoldDimension = m_InterpolationMatrix.rows();
  m_Thicknesses.set_size(m_ManifoldDimension);
  m_Nus.set_size(m_ManifoldDimension);
  m_Block.set_size(m_ManifoldDimension);
}

NetworkModel
::~NetworkModel()
{

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
NetworkModel
::Initialize(const Observations& Obs)
{
  /// Data-related attributes
  m_NumberOfSubjects          = Obs.GetNumberOfSubjects();
  m_IndividualObservationDate = Obs.GetObservations();
  m_SubjectTimePoints         = Obs.GetObservations();
  m_NbTotalOfObservations     = Obs.GetTotalNumberOfObservations();
  m_SumObservations           = Obs.GetTotalSumOfLandmarks();


  /// Noise
  m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.00001 );

  /// Population variables
  m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.13, 0.0001 * 0.0001});
  m_RealizationsPerRandomVariable["P"] = m_NbControlPoints;

  m_RandomVariables.AddRandomVariable("Nu", "Gaussian", {0.04088, 0.001*0.001});
  m_RealizationsPerRandomVariable["Nu"] = m_NbControlPoints;

  for(int i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
  {
      std::string Name = "Beta#" + std::to_string(i);
      m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.001*0.001});
      m_RealizationsPerRandomVariable.insert({Name, 1});
  }


  /// Individual variables
  m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
  m_RealizationsPerRandomVariable.insert({"Ksi", m_NumberOfSubjects});

  m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {75, 0.5});
  m_RealizationsPerRandomVariable.insert({"Tau", m_NumberOfSubjects});

  for(int i = 0; i < m_NbIndependentSources; ++i)
  {
      std::string Name = "S#" + std::to_string(i);
      m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0.0, 1});
      m_RealizationsPerRandomVariable.insert({Name, m_NumberOfSubjects});
  }
}

ScalarType
NetworkModel
::InitializePropositionDistributionVariance(std::string Name)
const
{
  Name = Name.substr(0, Name.find_first_of("#"));
  if("P" == Name)
      return 0.000001;
  if("Nu" == Name)
      return 0.000000006;
  if("Beta" == Name)
      return 0.000008*0.000008;
  if("Ksi" == Name)
      return 0.001;
  if("Tau" == Name)
      return 0.02 * 0.02;
  if("S" == Name)
      return 1;
}

void
NetworkModel
::UpdateModel(const Realizations &R, int Type,
            const std::vector<std::string, std::allocator<std::string>> Names)
{
  bool ComputeThickness = false;
  bool ComputeNu = false;
  bool ComputeBlock1 = false;
  bool ComputeBasis = false;
  bool ComputeA = false;
  bool ComputeSpace = false;

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
      else if(Name == "Nu")
      {
          IndividualOnly = false;
          ComputeNu = true;
          ComputeBlock1 = true;
          ComputeBasis = true;
          ComputeA = true;
          ComputeSpace = true;
      }
      else if(Name == "P")
      {
          IndividualOnly = false;
          ComputeThickness = true;
          ComputeBlock1 = true;
          ComputeBasis = true;
          ComputeA = true;
          ComputeSpace = true;
      }
      else if(Name == "Beta")
      {
          IndividualOnly = false;
          ComputeA = true;
          ComputeSpace = true;
          continue;
      }
      else if(Name == "S")
      {
          ComputeSpace= true;
      }
      else if(Name == "All")
      {
          IndividualOnly = false;
          ComputeSubjectTimePoint(R, -1);

          ComputeThickness = true;
          ComputeNu = true;
          ComputeBasis = true;
          ComputeA = true;
          ComputeSpace = true;
          ComputeBlock1 = true;
          break;
      }
      else
      {
          std::cerr << "The random variable name " << Name << "is unknown to the meshwork model" << std::endl;
      }
  }


  if(IndividualOnly)   ComputeSubjectTimePoint(R, Type);
  if(ComputeThickness) ComputeThicknesses(R);
  if(ComputeNu)        ComputeNus(R);
  if(ComputeBlock1)    ComputeBlock();
  if(ComputeBasis)     ComputeOrthonormalBasis();
  if(ComputeA)         ComputeAMatrix(R);
  if(ComputeSpace)     ComputeSpaceShifts(R);
}

AbstractModel::SufficientStatisticsVector
NetworkModel
::GetSufficientStatistics(const Realizations &R, const Observations& Obs)
{
  /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
  VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
  auto itS1 = S1.begin(), itS2 = S2.begin();
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      for(size_t j = 0; j < Obs.GetNumberOfTimePoints(i); ++j)
      {
          VectorType PC = ComputeParallelCurve(i, j);
          *itS1 = dot_product(PC, Obs.GetSubjectLandmark(i, j));
          *itS2 = PC.squared_magnitude();
          ++itS1, ++itS2;
      }
  }

  /// Sufficient Statistic P_k and P_k * P_k
  VectorType S3 = R.at("P");
  VectorType S4 = R.at("P") % R.at("P");
  
  /// Sufficient Statistic Nu_k and Nu_k*Nu_k
  VectorType S5 = R.at("Nu");
  VectorType S6 = R.at("Nu") % R.at("Nu");

  /// Sufficient statistic beta_k
  VectorType S7((m_ManifoldDimension-1) * m_NbIndependentSources);
  ScalarType * itS7 = S7.memptr();
  for(size_t i = 0; i < ((m_ManifoldDimension-1) * m_NbIndependentSources); ++i)
      itS7[i] = R.at("Beta#" + std::to_string(i), 0);

  /// Sufficient statistic Ksi_i * Ksi_i
  VectorType S8 = R.at("Ksi") % R.at("Ksi");


  /// Sufficient statistic Tau_i and Tau_i * Tau_i
  VectorType S9 = R.at("Tau");
  VectorType S10 = R.at("Tau") % R.at("Tau");

  return {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10};
}

void
NetworkModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS)
{
  /// Update the noise variance, sigma
  ScalarType NoiseVariance = m_SumObservations;
  const ScalarType * itS1 = SS[0].memptr();
  const ScalarType * itS2 = SS[1].memptr();
  for(size_t i = 0; i < SS[0].size(); ++i)
      NoiseVariance += - 2 * itS1[i] + itS2[i];

  NoiseVariance /= m_NbTotalOfObservations * m_ManifoldDimension;
  m_Noise->SetVariance(NoiseVariance);


  /// Update P_k : Mean and Var
  ScalarType PMean = 0.0, PVariance = 0.0;
  const ScalarType * itS3 = SS[2].memptr();
  const ScalarType * itS4 = SS[3].memptr();

  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
    PMean     += itS3[i];
    PVariance += itS4[i];
  }

  PMean     /= m_NbControlPoints;
  PVariance -= m_NbControlPoints * PMean * PMean;
  PVariance /= m_NbControlPoints;

  m_RandomVariables.UpdateRandomVariable("P", {{"Mean", PMean}, {"Variance", PVariance}});

  /// Update Nu_k : Mean and Var
  ScalarType NuMean = 0.0, NuVariance = 0.0;
  const ScalarType * itS5 = SS[4].memptr();
  const ScalarType * itS6 = SS[5].memptr();

  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
    NuMean     += itS5[i];
    NuVariance += itS6[i];
  }

  NuMean     /= m_NbControlPoints;
  NuVariance -= m_NbControlPoints * NuMean * NuMean;
  NuVariance /= m_NbControlPoints;

  m_RandomVariables.UpdateRandomVariable("Nu", {{"Mean", NuMean}, {"Variance", NuVariance}});

  /// Update Beta_k : Mean
  const ScalarType * itS7 = SS[6].memptr();
  for(size_t i = 0; i < SS[6].size(); ++i)
    m_RandomVariables.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", itS7[i]}});

  /// Update Ksi : Mean and Variance
  ScalarType KsiVariance = 0.0;
  const ScalarType * itS8 = SS[7].memptr();

  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
    KsiVariance += itS4[i];

  KsiVariance /= m_NumberOfSubjects;

  m_RandomVariables.UpdateRandomVariable("Ksi", {{"Variance", KsiVariance}});


  /// Update Tau : Mean and Variance
  ScalarType TauMean = 0.0, TauVariance = 0.0;
  const ScalarType * itS9  = SS[8].memptr();
  const ScalarType * itS10 = SS[9].memptr();

  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      TauMean     += itS9[i];
      TauVariance += itS10[i];
  }

  TauMean     /= m_NumberOfSubjects;
  TauVariance -= m_NumberOfSubjects * TauMean * TauMean;
  TauVariance /= m_NumberOfSubjects;

  m_RandomVariables.UpdateRandomVariable("Tau", {{"Mean", TauMean}, {"Variance", TauVariance}});



}

ScalarType
NetworkModel
::ComputeLogLikelihood(const Observations &Obs)
{
  double LogLikelihood = 0;

  for(size_t i = 0; i < m_NumberOfSubjects; ++i) 
  {
   ScalarType N = Obs.GetNumberOfTimePoints(i);
    
    for(size_t j = 0; j < N; ++j)
    {
      auto& it = Obs.GetSubjectLandmark(i, j);
      VectorType ParallelCurve = ComputeParallelCurve(i, j);
      LogLikelihood += (it - ParallelCurve).squared_magnitude();
    }
  }
  
  LogLikelihood /= -2*m_Noise->GetVariance();
  LogLikelihood -= m_NbTotalOfObservations*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
  
  return LogLikelihood;
}

ScalarType
NetworkModel
::ComputeIndividualLogLikelihood( const IndividualObservations& Obs, const int SubjectNumber)
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
NetworkModel
::SimulateData(io::DataSettings &DS)
{
  m_NumberOfSubjects = DS.GetNumberOfSimulatedSubjects();

  /// Initialize the realizations and simulate them
  m_RealizationsPerRandomVariable["P"] = m_NbControlPoints;

  m_RealizationsPerRandomVariable["Nu"] = m_NbControlPoints;

  for(size_t i = 0; i <  m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
      m_RealizationsPerRandomVariable["Beta#" + std::to_string(i)] = 1;

  m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
  m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;

  for(int i = 0; i < m_NbIndependentSources; ++i)
      m_RealizationsPerRandomVariable["S#" + std::to_string(i)] = m_NumberOfSubjects;

  auto R = SimulateRealizations();

  /// Update the model
  ComputeThicknesses(R);
  ComputeNus(R);
  ComputeBlock();
  ComputeOrthonormalBasis();
  ComputeAMatrix(R);
  ComputeSpaceShifts(R);

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
NetworkModel
::GetSamplerBlocks()
const
{
  int PopulationType = -1;

  std::vector<SamplerBlock> Blocks;

  /// Insert P
  MiniBlock P;
  for(size_t i = 0; i < m_NbControlPoints; ++i)
    P.push_back(std::make_pair("P", i));
  Blocks.push_back(std::make_pair(PopulationType, P));

  // Insert Nu
  MiniBlock Nu;
  for(size_t i = 0; i < m_NbControlPoints; ++i)
    Nu.push_back(std::make_pair("Nu", i));
  Blocks.push_back(std::make_pair(PopulationType, Nu));

  /// Insert Beta
  MiniBlock Beta;
  for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
  Blocks.push_back(std::make_pair(PopulationType, Beta));

  /// Individual variables
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
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
// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void
NetworkModel
::DisplayOutputs(const Realizations &R)
{
  auto P = m_RandomVariables.GetRandomVariable("P");
  auto PReal = R.at("P");

  auto Nu = m_RandomVariables.GetRandomVariable("Nu");
  auto NuReal = R.at("Nu");

  auto Tau = m_RandomVariables.GetRandomVariable("Tau");
  auto Ksi = m_RandomVariables.GetRandomVariable("Ksi");

  std::cout << "Noise: " << m_Noise->GetVariance() << " - PMean: " << exp(P->GetParameter("Mean"));
  std::cout << " - PVar: " << P->GetParameter("Variance") << " - PMin: " << exp(PReal.min_value()) << " - PMax: " << exp(PReal.max_value());
  std::cout << " - T0: " << Tau->GetParameter("Mean") << " - TauVar: " << Tau->GetParameter("Variance");
  std::cout << " - KsiVar: " << Ksi->GetParameter("Variance");
  std::cout << " - NuMean/V0: " << Nu->GetParameter("Mean") << " - NuVar: " << Nu->GetParameter("Variance");
  std::cout << " - NuMin: " << NuReal.min_value() << " - NuMax: " << NuReal.max_value() << std::endl;

}

void
NetworkModel
::SaveData(unsigned int IterationNumber, const Realizations &R)
{
  std::ofstream Outputs;
  std::string FileName = "/Users/igor.koval/Documents/Work/RiemAlzh/src/io/outputs/Network/Parameters" + std::to_string(IterationNumber) + ".txt";
  Outputs.open(FileName, std::ofstream::out | std::ofstream::trunc);

  /// Save the final noise variance
  Outputs << m_Noise->GetVariance() << std::endl;

  /// Save the number of subjects, the manifold dimension, the number of sources, and, the number of control points
  Outputs << m_NumberOfSubjects << ", " << m_ManifoldDimension << ", " << m_NbIndependentSources << ", " << m_NbControlPoints << std::endl;

  /// Save the thicknesses
  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
      Outputs << R.at("P", i);
      if(i != m_NbControlPoints - 1) { Outputs << ", ";}
  }
  Outputs << std::endl;

  /// Save the velocity nu
  for(size_t i = 0; i < m_NbControlPoints; ++i)
  {
      Outputs << R.at("Nu", i);
      if(i != m_NbControlPoints - 1) { Outputs << ", ";}
  }
  Outputs << std::endl;

      /// Save the tau
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      Outputs << R.at("Tau", i) ;
      if(i != m_NumberOfSubjects - 1) { Outputs << ", ";}
  }
  Outputs << std::endl;

  /// Save the ksi
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      Outputs << R.at("Ksi", i) ;
      if(i != m_NumberOfSubjects - 1) { Outputs << ", ";}
  }

  /// Save (S_i)
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      for(size_t j = 0; j < m_NbIndependentSources; ++j)
      {
          Outputs << R.at("S#" + std::to_string(j))(i);
          if(i != m_NbIndependentSources - 1) { Outputs << ", "; }
      }
      Outputs << std::endl;
  }

  /// Save (W_i)
  auto SizeW = m_NumberOfSubjects;
  for(size_t i = 0; i < m_NumberOfSubjects; ++i)
  {
      VectorType W = m_SpaceShifts.get_column(i);
      for(auto it = W.begin(); it != W.end(); ++it)
      {
          Outputs << *it;
          if(i != SizeW - 1) { Outputs << ", "; }
      }
      Outputs << std::endl;
  }


}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void
NetworkModel
::InitializeFakeRandomVariables()
{
  /// Noise
  m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.00001 );

  /// Population variables
  m_RandomVariables.AddRandomVariable("P",  "Gaussian", {0.13, 0.00005 * 0.00005});
  m_RandomVariables.AddRandomVariable("Nu", "Gaussian", {0.04088, 0.001*0.001});

  for(int i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
      m_RandomVariables.AddRandomVariable("Beta#" + std::to_string(i), "Gaussian", {0, 0.001*0.001});

  /// Individual variables
  m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
  m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {75, 0.025});

  for(int i = 0; i < m_NbIndependentSources; ++i)
      m_RandomVariables.AddRandomVariable("S#" + std::to_string(i), "Gaussian", {0.0, 1});

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
NetworkModel
::ComputeSubjectTimePoint(const Realizations &R, const int SubjectNumber)
{
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
NetworkModel
::ComputeNus(const Realizations &R)
{
  m_Nus = m_InterpolationMatrix * m_InvertKernelMatrix * R.at("Nu");
}

void
NetworkModel
::ComputeThicknesses(const Realizations &R)
{
  m_Thicknesses = m_InterpolationMatrix * m_InvertKernelMatrix * R.at("P").exp();
}

void
NetworkModel
::ComputeOrthonormalBasis()
{

  VectorType U(m_ManifoldDimension);
  ScalarType * u = U.memptr();
  ScalarType * t = m_Thicknesses.memptr();
  ScalarType * n = m_Nus.memptr();

  for(size_t i = 0; i < m_ManifoldDimension; ++i)
      u[i] = n[i] / (t[i] * t[i]);

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
NetworkModel
::ComputeAMatrix(const Realizations &R)
{
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
NetworkModel
::ComputeSpaceShifts(const Realizations &R)
{
  MatrixType SS(m_NbIndependentSources, m_NumberOfSubjects);
  for(int i = 0; i < m_NbIndependentSources; ++i)
  {
      SS.set_row(i, R.at("S#" + std::to_string(i)));
  }
  m_SpaceShifts = m_AMatrix * SS;
}

void
NetworkModel
::ComputeBlock()
{
  ScalarType * b = m_Block.memptr();
  ScalarType * n = m_Nus.memptr();
  ScalarType * t = m_Thicknesses.memptr();

  for(size_t i = 0; i < m_ManifoldDimension; ++i)
      b[i] = n[i] / (t[i] * t[i]);
}

AbstractModel::VectorType
NetworkModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber)
{
  VectorType ParallelCurve(m_ManifoldDimension);
  ScalarType * p = ParallelCurve.memptr();

  double TimePoint = m_SubjectTimePoints[SubjectNumber](ObservationNumber);
  ScalarType * t = m_Thicknesses.memptr();
  ScalarType * w = m_SpaceShifts.get_column(SubjectNumber).memptr();
  ScalarType * b = m_Block.memptr();

  for(size_t i = 0; i < m_ManifoldDimension; ++i)
      p[i] = t[i] * exp(w[i] + b[i]*TimePoint);


  return ParallelCurve;
}
