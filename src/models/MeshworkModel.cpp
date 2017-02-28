#include "MeshworkModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

MeshworkModel
::MeshworkModel(io::ModelSettings &MS) 
{
    m_ManifoldDimension = MS.GetManifoldDimension();
    m_NbIndependentSources= MS.GetNumberOfIndependentSources();
    
    std::string KernelMatrixPath = MS.GetInvertKernelPath();
    std::string InterpolationMatrixPath = MS.GetInterpolationKernelPath();

    m_InvertKernelMatrix = io::ReadData::OpenKernel(KernelMatrixPath).transpose();
    m_InterpolationMatrix = io::ReadData::OpenKernel(InterpolationMatrixPath);
    
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_Thicknesses.set_size(m_ManifoldDimension);
    m_Deltas.set_size(m_ManifoldDimension);
    m_Block1.set_size(m_ManifoldDimension);
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
    m_NumberOfSubjects          = Obs.GetNumberOfSubjects();
    m_IndividualObservationDate = Obs.GetObservations();
    m_SubjectTimePoints         = Obs.GetObservations();
    m_NbTotalOfObservations     = Obs.GetTotalNumberOfObservations();
    m_SumObservations           = Obs.GetTotalSumOfLandmarks();
    
    /// Population variables
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.01 );

    m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.13, 0.00005 * 0.00005});
    m_RealizationsPerRandomVariable["P"] = m_NbControlPoints;
    
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        std::string Name = "Delta#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.003 * 0.003});
        m_RealizationsPerRandomVariable[Name] = 1;
    }
    
    for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    {
        std::string Name = "Beta#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.001 * 0.001});
        m_RealizationsPerRandomVariable[Name] = 1;
    }
    
    /// Individual variables
    m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {-3.1971, 0.000000004});
    m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
    
    m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
    m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
        
    for(int i = 0; i < m_NbIndependentSources; ++i)
    {
        std::string Name = "S#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0.0, 1});
        m_RealizationsPerRandomVariable[Name] =  m_NumberOfSubjects;
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
    
    /// Sufficient Statistic Ksi and Ksi*Ksi
    VectorType S3 = R.at("Ksi");
    VectorType S4 = R.at("Ksi") % R.at("Ksi");
    
    /// Sufficient statistic Tau and Tau*Tau
    VectorType S5 = R.at("Tau");
    VectorType S6 = R.at("Tau") % R.at("Tau");
    
    /// Sufficient statistic beta_k
    VectorType S7((m_ManifoldDimension-1) * m_NbIndependentSources);
    ScalarType * itS7 = S7.memptr();
    for(size_t i = 0; i < ((m_ManifoldDimension-1) * m_NbIndependentSources); ++i)
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
    ScalarType NoiseVariance = m_SumObservations;
    const ScalarType * itS1 = SS[0].memptr();
    const ScalarType * itS2 = SS[1].memptr();
    for(size_t i = 0; i < SS[0].size(); ++i)
        NoiseVariance += - 2 * itS1[i] + itS2[i];
    
    NoiseVariance /= m_NbTotalOfObservations * m_ManifoldDimension;
    m_Noise->SetVariance(NoiseVariance);
    
    /// Update Ksi : Mean and Variance
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
    
    m_RandomVariables.UpdateRandomVariable("Ksi", {{"Mean", KsiMean}, {"Variance", KsiVariance}});
    
    
    /// Update Tau : Mean and Variance
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
    
    
    /// Update Beta_k : Mean
    const ScalarType * itS7 = SS[6].memptr();
    for(size_t i = 0; i < SS[6].size(); ++i)
        m_RandomVariables.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", itS7[i]}});
     
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
    
    m_RandomVariables.UpdateRandomVariable("P", {{"Mean", PMean}, {"Variance", PVariance}});
    
    /// Update Delta_k : Mean
    const ScalarType * itS10 = SS[9].memptr();
    for(size_t i = 0; i < m_NbControlPoints - 1; ++i)
        m_RandomVariables.UpdateRandomVariable("Delta#" + std::to_string(i+1), {{"Mean", itS10[i]}});
}

ScalarType 
MeshworkModel
::ComputeLogLikelihood(const OldData &D) 
{
    double LogLikelihood = 0;
//#pragma omp parallel for reduction(+:LogLikelihood)   
    for(size_t j = 0; j < m_NumberOfSubjects; ++j) 
    {
       double N = D.at(j).size();

        for(size_t i = 0; i < N; ++i)
        {
            auto& it = D.at(j).at(i);
            VectorType ParallelCurve = ComputeParallelCurve(j, i);
            LogLikelihood += (it.first - ParallelCurve).squared_magnitude();
        }
        
    }
    
    LogLikelihood  /= -2*m_Noise->GetVariance();
    LogLikelihood -= m_NbTotalOfObservations*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
    
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

MeshworkModel::OldData
MeshworkModel
::SimulateData(io::DataSettings &DS) 
{
    typedef std::vector< std::pair< VectorType, double> > IndividualData;
    
    m_NumberOfSubjects = DS.GetNumberOfSimulatedSubjects();
    
    /// Initialize the realizations and simulate them
    m_RealizationsPerRandomVariable["P"] = m_NbControlPoints;
    
    for(int i = 1; i < m_NbControlPoints; ++i)
        m_RealizationsPerRandomVariable["Delta#" + std::to_string(i)] = 1;
    
    
    for(size_t i = 0; i <  m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
        m_RealizationsPerRandomVariable["Beta#" + std::to_string(i)] = 1;
        
    m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
    m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
    
    for(int i = 0; i < m_NbIndependentSources; ++i)
        m_RealizationsPerRandomVariable["S#" + std::to_string(i)] = m_NumberOfSubjects;
    
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
    std::uniform_real_distribution<double> ObsDistrib(60, 95);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    OldData D;
    double RealNoise = 0.0;
    
    for(int i = 0; i < m_NumberOfSubjects; ++i)
    {
        IndividualData ID;
        std::vector<double> TimePoints;
        for(int j = 0; j < Uni(RNG); ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));
        }
        std::sort(TimePoints.begin(), TimePoints.end());
                
        VectorType T(TimePoints.size());
        int y = 0;
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it, ++y) 
        {
            T(y) = *it;
        }
        
        m_SubjectTimePoints.push_back(T);
        
        int ObsNb = 0;
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it, ++ObsNb)
        {
            VectorType Scores(m_ManifoldDimension);
            VectorType ParallelCurve = ComputeParallelCurve(i, ObsNb);
            for(    auto IterCurve = ParallelCurve.begin(), IterScore = Scores.begin() 
                    ; IterCurve != ParallelCurve.end() && IterScore!= Scores.end()
                    ; ++IterCurve, ++IterScore)
            {
                double Noise = NoiseDistrib(RNG);
                RealNoise += Noise;
                *IterScore = *IterCurve + Noise;
            }
            ID.push_back(std::pair<VectorType, double>(Scores, *it));
        }
        D.push_back(ID);
    }
    
    std::vector<VectorType> IndividualObservationDate;
    double SumObservations = 0.0, K = 0.0;
    unsigned int Indiv = 0;
    for(auto it = D.begin(); it != D.end(); ++it, ++Indiv)
    {
        VectorType IndividualObs(it->size());
        K += it->size();
        int i = 0;
        for(auto it2 = it->begin(); it2 != it->end(); ++it2, ++i)
        {
            SumObservations += it2->first.squared_magnitude();
            IndividualObs(i) = it2->second;
        }
        IndividualObservationDate.push_back(IndividualObs);
    }
    m_IndividualObservationDate = IndividualObservationDate;
    m_SumObservations = SumObservations;
    m_NbTotalOfObservations = K;
    
    std::cout << "Real Likelihood = " << ComputeLogLikelihood(D) << std::endl;
    
    return D;
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
    int BetaModulo = (int)m_NbIndependentSources*(m_ManifoldDimension - 1) /NbBeta;
    for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    {
        Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
        bool HingeCondition = (i%BetaModulo == 0 && i != 0);
        bool FinalCondition = (i == m_NbIndependentSources*(m_ManifoldDimension - 1) - 1);
        if(FinalCondition || HingeCondition)
        {
            Blocks.push_back(std::make_pair(PopulationType, Beta));
            Beta.clear();
        }
    }
     */
    MiniBlock Beta;
    for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
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
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////
    

void 
MeshworkModel
::DisplayOutputs(const Realizations &R) 
{
    auto P = m_RandomVariables.GetRandomVariable("P");
    auto PReal = R.at("P");
    
    auto Tau = m_RandomVariables.GetRandomVariable("Tau");
    auto Ksi = m_RandomVariables.GetRandomVariable("Ksi");
    
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
    Outputs << m_NumberOfSubjects << ", " << m_ManifoldDimension << ", " << m_NbIndependentSources << ", " << m_NbControlPoints << std::endl;
    
    /// Save the delta_mean -> First one being equal to 0 as the reference
    Outputs << 0 << ", ";
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        Outputs << m_RandomVariables.GetRandomVariable("Delta#" + std::to_string(i))->GetParameter("Mean");
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


void
MeshworkModel
::InitializeFakeRandomVariables() 
{
    /// Noise 
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.00001 );
    
    /// Population random variables
    m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.1, 0.001 * 0.001});
    
    for(size_t i = 0; i < m_NbControlPoints; ++i)
        m_RandomVariables.AddRandomVariable("Delta#" + std::to_string(i), "Gaussian", {0, 0.0001 * 0.0001});
    
    for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
        m_RandomVariables.AddRandomVariable("Beta#" + std::to_string(i), "Gaussian", {0, 0.0001 * 0.0001});
    
    
    /// Individual random variables
    m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
    m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {62, 0.25});
    
    for(int i = 0; i < m_NbIndependentSources; ++i)
        m_RandomVariables.AddRandomVariable("S#" + std::to_string(i), "Gaussian", {0.0, 0.5});

    
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
    auto V0 = m_RandomVariables.GetRandomVariable("Ksi")->GetParameter("Mean");
    V0 = exp(V0);
    
    
    VectorType U(m_ManifoldDimension);
    ScalarType * t = m_Thicknesses.memptr();
    ScalarType * d = m_Deltas.memptr();
    ScalarType * u = U.memptr();
    
#pragma omp simd
    for(int i = 0; i < m_ManifoldDimension; ++i)
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
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        FinalMatrix2(i, i ) += 1;
    

    m_OrthogonalBasis = FinalMatrix2;
}

void 
MeshworkModel
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
MeshworkModel
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
MeshworkModel
::ComputeBlock() 
{
    ScalarType * t = m_Thicknesses.memptr();
    ScalarType * d = m_Deltas.memptr();
    ScalarType * b = m_Block1.memptr();

#pragma omp simd
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        b[i] = 1 / (t[i] * exp(d[i]));
}


AbstractModel::VectorType
MeshworkModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
    double TimePoint = m_SubjectTimePoints[SubjectNumber](ObservationNumber);
    
    VectorType ParallelCurve(m_ManifoldDimension);
    
    ScalarType * p = ParallelCurve.memptr();
    ScalarType * t = m_Thicknesses.memptr();
    ScalarType * d = m_Deltas.memptr();
    ScalarType * b = m_Block1.memptr();
    ScalarType * w = m_SpaceShifts.get_column(SubjectNumber).memptr();
    
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        p[i] = t[i] * exp( w[i] / (t[i]*exp(d[i])) + d[i] - TimePoint/t[i]);
        
    
    return ParallelCurve;
}





