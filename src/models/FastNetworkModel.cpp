#include "FastNetworkModel.h"
#include <armadillo>


FastNetworkModel
::FastNetworkModel(io::ModelSettings &MS) 
{
    m_ManifoldDimension = MS.GetManifoldDimension();
    m_NbIndependentComponents = MS.GetNumberOfIndependentSources();
    
    std::string KernelMatrixPath = MS.GetInvertKernelPath();
    std::string InterpolationMatrixPath = MS.GetInterpolationKernelPath();

    m_InvertKernelMatrix = io::ReadData::OpenKernel(KernelMatrixPath).transpose();
    m_InterpolationMatrix = io::ReadData::OpenKernel(InterpolationMatrixPath);
    
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_Nus.set_size(m_ManifoldDimension);
    m_Deltas.set_size(m_ManifoldDimension);
    m_Block1.set_size(m_ManifoldDimension);
    m_Block2.set_size(m_ManifoldDimension);
}


FastNetworkModel
::FastNetworkModel(const unsigned int NbIndependentComponents,
                          std::shared_ptr<MatrixType> KernelMatrix,
                          std::shared_ptr<MatrixType> InterpolationMatrix) 
{
    m_NbIndependentComponents = NbIndependentComponents;
    m_InvertKernelMatrix = KernelMatrix->transpose();
    m_InterpolationMatrix = *InterpolationMatrix;
    
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_Nus.set_size(m_ManifoldDimension);
    m_Deltas.set_size(m_ManifoldDimension);
    m_Block1.set_size(m_ManifoldDimension);
    m_Block2.set_size(m_ManifoldDimension);

}


FastNetworkModel
::~FastNetworkModel() 
{
    
   
    
}

void
FastNetworkModel
::Initialize(const OldData& D) 
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    /// Data-related attributes
    m_NumberOfSubjects = D.size();
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
    m_SubjectTimePoints = IndividualObservationDate;
    m_SumObservations = SumObservations;
    m_NbTotalOfObservations = K;
    
    
    
     /// Population variables
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.000001 );

    m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.12, 0.0001 * 0.0001});
    m_RealizationsPerRandomVariable.insert({"P", 1});
        
    for(int i = 1; i < m_NbControlPoints; ++i)
    {
        std::string NameDelta = "Delta#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(NameDelta, "Gaussian", {0, 0.003 * 0.003});
        m_RealizationsPerRandomVariable.insert({NameDelta, 1});
    }
    
  
    m_RandomVariables.AddRandomVariable("Nu", "Gaussian", {0.04088, 0.001*0.001});
    m_RealizationsPerRandomVariable["Nu"] = m_NbControlPoints;
    
    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
    {
        std::string Name = "Beta#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.001*0.001});
        m_RealizationsPerRandomVariable.insert({Name, 1});
    }
    
    
    /// Individual variables
    m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
    m_RealizationsPerRandomVariable.insert({"Ksi", m_NumberOfSubjects});
    
    m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {75, 0.025});
    m_RealizationsPerRandomVariable.insert({"Tau", m_NumberOfSubjects});
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        std::string Name = "S#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0.0, 1});
        m_RealizationsPerRandomVariable.insert({Name, m_NumberOfSubjects});
    }
    

    
}



ScalarType 
FastNetworkModel
::InitializePropositionDistributionVariance(std::string Name) 
const
{
    Name = Name.substr(0, Name.find_first_of("#"));
    if("P" == Name)
        return 0.000001;
    if("Delta" == Name)
        return 0.00000001;
    if("Nu" == Name)
        return 0.00001;
    if("Beta" == Name)
        return 0.00001*0.00001;
    if("Ksi" == Name)
        return 0.0008;
    if("Tau" == Name)
        return 0.02 * 0.02;
    if("S" == Name)
        return 0.7;
}

void 
FastNetworkModel
::UpdateModel(const Realizations& AR, int Type,
              const std::vector<std::string> Names) 
{
    
    bool ComputePosition = false;
    bool ComputeDelta = false;
    bool ComputeNu = false;
    bool ComputeBasis = false;
    bool ComputeA = false;
    bool ComputeSpaceShift = false;
    bool ComputeBlock_1 = false;
    bool ComputeBlock_2 = false;
    
    bool IndividualOnly = true;
    if(Type == -1) IndividualOnly = false;
        
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
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_2 = true;
            continue;
        }
        else if(Name == "Delta")
        {
            IndividualOnly = false;
            ComputeDelta = true;
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_1 = true;
            continue;
        }
        else if(Name == "P")
        {
            IndividualOnly = false;
            ComputePosition = true;
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_1 = true;
            ComputeBlock_2 = true;
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
            continue;
        }
        else if(Name == "All")
        {
            ComputeSubjectTimePoint(AR, -1);
            IndividualOnly = false;
            ComputePosition = true;
            ComputeDelta = true;
            ComputeNu = true;
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_1 = true;
            ComputeBlock_2 = true;
            break;
        }
        else
        {
            std::cerr << "PROBLEM WITH FAST NETWORK MODEL" << std::endl;
        }
    }
    
    // TODO : To parse it even faster, update just the coordinates within the names
    if(IndividualOnly) ComputeSubjectTimePoint(AR, Type);
    
    if(ComputePosition) { 
        m_P0 = exp(AR.at("P", 0)); 
        //std::cout << "m_P0 " << m_P0 << " & " << AR.at("P", 0) <<  std::endl;
        //std::cout << "P0 : " << m_RandomVariables.GetRandomVariable("P")->GetParameter(0) << " & " << m_RandomVariables.GetRandomVariable("P")->GetParameter(1) << std::endl;
        //int a = 0;
    }
    if(ComputeDelta) ComputeDeltas(AR);
    if(ComputeNu) ComputeNus(AR);
    if(ComputeBasis) ComputeOrthonormalBasis();
    if(ComputeA) ComputeAMatrix(AR);
    if(ComputeSpaceShift) ComputeSpaceShifts(AR);
    if(ComputeBlock_1) ComputeBlock1();
    if(ComputeBlock_2) ComputeBlock2();
    
}

FastNetworkModel::OldData
FastNetworkModel
::SimulateData(io::DataSettings& DS) 
{
typedef std::vector< std::pair< VectorType, double> > IndividualData;
    
    m_NumberOfSubjects = DS.GetNumberOfSimulatedSubjects();
    
    
    /// Initialize the realizations and simulate them
    m_RealizationsPerRandomVariable["P"] = 1;
    
    for(int i = 1; i < m_NbControlPoints; ++i)
        m_RealizationsPerRandomVariable["Delta#" + std::to_string(i)] = 1;
    
    m_RealizationsPerRandomVariable["Nu"] = m_NbControlPoints;

    
    for(size_t i = 0; i <  m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
        m_RealizationsPerRandomVariable["Beta#" + std::to_string(i)] = 1;
        
    m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
    m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
        m_RealizationsPerRandomVariable["S#" + std::to_string(i)] = m_NumberOfSubjects;
    
    auto R = SimulateRealizations();
    
    
    /// Update the model
    m_P0 = exp(R.at("P")(0));
    ComputeDeltas(R);
    ComputeNus(R);
    ComputeOrthonormalBasis();
    ComputeAMatrix(R);
    ComputeSpaceShifts(R);
    ComputeBlock1();
    ComputeBlock2();
    
    
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


double 
FastNetworkModel
::ComputeLogLikelihood(const OldData& D) 
{
    double LogLikelihood = 0;
#pragma omp parallel for reduction(+:LogLikelihood)   
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

double
FastNetworkModel
::ComputeIndividualLogLikelihood(const OldData& D, const int SubjectNumber) 
{
    /// Get the data
    double LogLikelihood = 0;
    auto N = D.at(SubjectNumber).size();
    auto Did = D.at(SubjectNumber);
    
#pragma omp parallel for reduction(+:LogLikelihood)   
    for(size_t i = 0; i < N; ++i)
    {
        auto& it = Did.at(i);
        VectorType P2 = ComputeParallelCurve(SubjectNumber, i);
        LogLikelihood += (it.first - P2).squared_magnitude();
    }
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
    
    return LogLikelihood;
}

FastNetworkModel::SufficientStatisticsVector
FastNetworkModel
::GetSufficientStatistics(const Realizations& R, const OldData& D) 
{
    
    /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
    VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
    auto itS1 = S1.begin(), itS2 = S2.begin();
    int i = 0;
    for(auto itD = D.begin(); itD != D.end(); ++itD, ++i)
    {        
        int j = 0;
        for(auto itD2 = itD->begin(); itD2 != itD->end(); ++itD2, ++j)
        {
            VectorType P2 = ComputeParallelCurve(i, j);
            *itS1 = dot_product(P2, itD2->first);
            *itS2 = P2.squared_magnitude();
            ++itS1, ++itS2;
        }
    }
    
    /// S3 <- Ksi_i * Ksi_i
    VectorType S3 = R.at("Ksi") % R.at("Ksi");

    
    /// S4 <- Tau_i   &    S5 <- Tau_i * Tau_i
    VectorType S4 = R.at("Tau");
    VectorType S5 = R.at("Tau") % R.at("Tau");
    
    /// S6 <- P0 
    VectorType S6(1, R.at("P", 0));
    
    /// S7 <- beta_k
    VectorType S7((m_ManifoldDimension-1) * m_NbIndependentComponents);
    i = 0;
    for(auto it = S7.begin(); it != S7.end(); ++it, ++i)
    {
        *it = R.at("Beta#" + std::to_string(i), 0);
    }
    
    /// S8 <- delta_k
    VectorType S8(m_NbControlPoints - 1);
    i = 1;
    for( auto itS8 = S8.begin(); itS8 != S8.end(); ++itS8, ++i)
    {
        *itS8 = R.at("Delta#" + std::to_string(i), 0);
        
    }
    
    /// S9 <- nu_k, S10 = nu_k * nu_k
    VectorType S9 = R.at("Nu");
    VectorType S10 = R.at("Nu") % R.at("Nu");
    
    
    SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10};
    return S;
}   

void
FastNetworkModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS, const OldData& D) 
{
    /// Update sigma
    double NoiseVariance = m_SumObservations;
    for(auto itS1 = SS[0].begin(), itS2 = SS[1].begin(); itS1 != SS[0].end() && itS2!= SS[1].end(); ++itS1, ++itS2)
    {
        NoiseVariance += -2* *itS1 + *itS2;
    }
    NoiseVariance /= m_NbTotalOfObservations * m_ManifoldDimension;
    m_Noise->SetVariance(NoiseVariance);
    
    /// Update ksi and sigma_ksi
    double KsiVariance = 0.0;
    for(auto it = SS[2].begin(); it != SS[2].end(); ++it)
    {
        KsiVariance += *it;
    }
    KsiVariance /= m_NumberOfSubjects;
    
    m_RandomVariables.UpdateRandomVariable("Ksi", {{"Variance", KsiVariance}});
    
    /// Update tau and sigma_tau
    double TauMean = 0.0, TauVariance = 0.0;
    for(auto it = SS[3].begin(); it != SS[3].end(); ++it)
    {
        TauMean += *it;
    }
    TauMean /= m_NumberOfSubjects;
    for(auto it = SS[4].begin(); it != SS[4].end(); ++it)
    {
        TauVariance += *it;
    }
    TauVariance -= m_NumberOfSubjects * TauMean * TauMean;
    TauVariance /= m_NumberOfSubjects;
    
    m_RandomVariables.UpdateRandomVariable("Tau", {{"Mean", TauMean}, {"Variance", TauVariance}});
    
    /// Update P0
    m_RandomVariables.UpdateRandomVariable("P", {{"Mean", SS[5](0)}});
    
    /// Update beta_k
    int i = 0;
    for(auto it = SS[6].begin(); it != SS[6].end(); ++it, ++i)
    {
        m_RandomVariables.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", *it}});
    }
        
    /// Update delta_k
    i = 1;
    for(auto it = SS[7].begin(); it != SS[7].end(); ++it, ++i)
    {
        m_RandomVariables.UpdateRandomVariable("Delta#" + std::to_string(i), {{"Mean", *it}});
    }
    
    /// Update nu_k = v0 and sigma_nu
    double V0 = 0.0, NuVariance = 0.0;
    for(auto it = SS[8].begin(); it != SS[8].end(); ++it)
    {
        V0 += *it;
    }
    V0 /= m_NbControlPoints;
    for(auto it = SS[9].begin(); it != SS[9].end(); ++it)
    {
        NuVariance += *it;
    }
    NuVariance -= m_NbControlPoints * V0 * V0;
    NuVariance /= m_NbControlPoints;

    m_RandomVariables.UpdateRandomVariable("Nu", {{"Mean", V0}, {"Variance", NuVariance}});
}


std::vector<AbstractModel::SamplerBlock>
FastNetworkModel
::GetSamplerBlocks() 
const 
{
    int PopulationType = -1;
    int NbBeta = 1;
    int NbDelta = 0;
    int NbNu = 0;
    
    
    std::vector<SamplerBlock> Blocks;
    
    /// Insert P0;
    MiniBlock P = {std::make_pair("P", 0)};
    Blocks.push_back(std::make_pair(PopulationType, P));
    
    /// Insert Beta_k
    /*
    MiniBlock Beta;
    int BetaModulo = (int)m_NbIndependentComponents*(m_ManifoldDimension - 1) /NbBeta;
    for(size_t i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
    {
        Beta.push_back(std::make_pair("Beta#" + std::to_string(i), 0));
        bool HingeCondition = (i%BetaModulo == 0 && i != 0);
        bool FinalCondition = (i == m_NbIndependentComponents*(m_ManifoldDimension - 1) - 1);
        if(FinalCondition || HingeCondition)
        {
            Blocks.push_back(std::make_pair(PopulationType, Beta));
            Beta.clear();
        }
    }
    */
    

    MiniBlock Beta;
    for(size_t i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
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
    
    /// Insert Nu_k
    MiniBlock Nu;
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        Nu.push_back(std::make_pair("Nu", i));
    }
    Blocks.push_back(std::make_pair(PopulationType, Nu));
    
    /// Individual variables
    
    for(size_t i = 0; i < m_NumberOfSubjects; ++i)
    {
        MiniBlock IndividualBlock;
        IndividualBlock.push_back(std::make_pair("Ksi", i));
        IndividualBlock.push_back(std::make_pair("Tau", i));
        for(size_t j = 0; j < m_NbIndependentComponents; ++j)
            IndividualBlock.push_back(std::make_pair("S#" + std::to_string(j), i));
        
        Blocks.push_back(std::make_pair(i, IndividualBlock));
    }
    
    /*
    /// Individual variables bis
    MiniBlock TauBlock;
    for(size_t i = 0; i < m_NumberOfSubjects; ++i) {
        TauBlock.push_back(std::make_pair("Tau", i));
        Blocks.push_back(std::make_pair(i, TauBlock));
        TauBlock.clear();
    }
    
    MiniBlock KsiBlock;
    for(size_t i = 0; i < m_NumberOfSubjects; ++i) {
        KsiBlock.push_back(std::make_pair("Ksi", i));
        Blocks.push_back(std::make_pair(i, KsiBlock));
        KsiBlock.clear();
    }
    
    for(size_t j = 0; j < m_NbIndependentComponents; ++j)
    {
        MiniBlock SBlock;
        for(size_t i = 0; i < m_NumberOfSubjects; ++i) {
            SBlock.push_back(std::make_pair("S#" + std::to_string(j), i));
            Blocks.push_back(std::make_pair(i, SBlock));
            SBlock.clear();
        }
    }
     */
    
    return Blocks;
}

void 
FastNetworkModel
::DisplayOutputs(const Realizations& R) 
{
    auto P0 = m_RandomVariables.GetRandomVariable("P")->GetParameter("Mean");
    auto Tau = m_RandomVariables.GetRandomVariable("Tau");
    
    auto Nu = m_RandomVariables.GetRandomVariable("Nu"); 

    double NuMax = R.at("Nu").max_value();
    double NuMin = R.at("Nu").min_value();

    
    double DeltaMin = R.at("Delta#1", 0);
    double DeltaMax = DeltaMin;
    for(size_t i = 1; i < 258; ++i)
    {
        double DeltaK = R.at("Delta#" + std::to_string(i), 0);
        DeltaMax = std::max(DeltaMax, DeltaK);
        DeltaMin = std::min(DeltaMin, DeltaK);
    }
    
    std::cout << "Noise: " << m_Noise->GetVariance();
    std::cout << " - p0: " << exp(P0) << " - t0: " << Tau->GetParameter("Mean") << " - S_tau:" << Tau->GetParameter("Variance");
    std::cout << " - v0: " << Nu->GetParameter("Mean") << " - S_nu:" << Nu->GetParameter("Variance");
    std::cout << " - MaxNu: " << NuMax << " - MinNu: " << NuMin;
    std::cout << " - MaxDelta: " << DeltaMax << " - MinDelta: " << DeltaMin << std::endl;
    
}


void 
FastNetworkModel
::SaveData(unsigned int IterationNumber, const Realizations& R) 
{
    
    unsigned int NumberOfSubjects = R.at("Tau").size();
    std::ofstream Outputs;    
    std::string FileName = "/Users/igor.koval/Documents/Work/RiemAlzh/src/io/outputs/Simulated/Parameters" + std::to_string(IterationNumber) + ".txt";
    Outputs.open(FileName, std::ofstream::out | std::ofstream::trunc);
    
    /// Save the final noise variance
    Outputs << m_Noise->GetVariance() << std::endl;
    
   /// Save Number of subject, Dimensions, Number of Sources, Number of control points
    Outputs << NumberOfSubjects << ", " << m_ManifoldDimension << ", " << m_NbIndependentComponents  << ", " << m_NbControlPoints << std::endl;
    
    /// Save P0, P0_mean and P0_var
    auto P0 = m_RandomVariables.GetRandomVariable("P");
    Outputs << R.at("P")(0) << ", " << P0->GetParameter("Mean") << ", " << P0->GetParameter("Variance") << std::endl;
    // std::cout  << R.at("P")(0) << ", " << P0->GetParameter("Mean") << ", " << P0->GetParameter("Variance") << std::endl;
    
    /// Save (Ksi_i)
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        Outputs << R.at("Ksi")(i);
        if(i != NumberOfSubjects - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save Ksi_mean and Ksi_Var
    auto Ksi = m_RandomVariables.GetRandomVariable("Ksi");
    Outputs << Ksi->GetParameter("Mean") << ", " << Ksi->GetParameter("Variance") << std::endl;
    
    /// Save V0 
    auto Nu = m_RandomVariables.GetRandomVariable("Nu");
    Outputs << Nu->GetParameter("Mean") << std::endl;
    
    /// Save (Tau_i)
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        Outputs << R.at("Tau")(i);
        if(i != NumberOfSubjects - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save Tau_Mean and Tau_Var
    auto Tau = m_RandomVariables.GetRandomVariable("Tau");
    Outputs << Tau->GetParameter("Mean") << ", " << Tau->GetParameter("Variance") << std::endl;
    
    /// Save (Delta_tilde_k)
    Outputs << 0 << ", ";
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        Outputs << R.at("Delta#" + std::to_string(i))(0);
        if(i != m_NbControlPoints - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Delta_k)
    Outputs << 0 << ", ";
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        std::string Name = "Delta#" + std::to_string(i);
        Outputs << m_RandomVariables.GetRandomVariable(Name)->GetParameter("Mean");
        if(i != m_NbControlPoints - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Nu_k)
    for(size_t i = 0; i < m_NbControlPoints; ++i)  
    {
        Outputs << R.at("Nu", i);  
        if(i != m_NbControlPoints - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Nu_k_mean)
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        Outputs << m_RandomVariables.GetRandomVariable("Nu")->GetParameter("Mean");
        if(i != m_NbControlPoints - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (S_i)
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        for(size_t j = 0; j < m_NbIndependentComponents; ++j)
        {
            Outputs << R.at("S#" + std::to_string(j))(i);
            if(i != m_NbIndependentComponents - 1) { Outputs << ", "; }
        }
        Outputs << std::endl;
    }
    
    /// Save (W_i)
    auto SizeW = m_NumberOfSubjects;
    for(size_t i = 0; i < NumberOfSubjects; ++i)
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
FastNetworkModel
::InitializeFakeRandomVariables() 
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
       
    
     /// Population variables
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.1 );

    m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.08, 0.00001 * 0.00001});
        
    for(int i = 1; i < m_NbControlPoints; ++i)
        m_RandomVariables.AddRandomVariable("Delta#" + std::to_string(i), "Gaussian", {0.5, 0.0001 * 0.0001});

    m_RandomVariables.AddRandomVariable("Nu", "Gaussian", {0.036, 0.004*0.004});
    
    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
        m_RandomVariables.AddRandomVariable("Beta#" + std::to_string(i), "Gaussian", {0, 0.0001 * 0.0001});
    
    
    /// Individual variables    
    m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
    
    m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {70, 0.25});
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
        m_RandomVariables.AddRandomVariable("S#" + std::to_string(i), "Gaussian", {0.0, 0.5});
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
FastNetworkModel
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
    /*
    if(SubjectNumber != -1) {
        double AccFactor = exp(R.at("Ksi", SubjectNumber));
        double TimeShift = R.at("Tau", SubjectNumber);
        
        auto N = m_IndividualObservationDate[SubjectNumber].size();
        
        ScalarType * real = m_IndividualObservationDate[SubjectNumber].memptr();
        ScalarType * reparam = m_SubjectTimePoints[SubjectNumber].memptr();
    
#pragma omp simd
        for(size_t i = 0; i < N; ++i)
            reparam[i] = AccFactor * (real[i] - TimeShift);
        
    }
    else
    {
#pragma parallel for
        for(size_t i = 0; i < m_NumberOfSubjects; ++i)
        {
            double AccFactor = exp(R.at("Ksi")(i));
            double TimeShift = R.at("Tau")(i);
            
            auto N = m_IndividualObservationDate[i].size();
            
            ScalarType * real = m_IndividualObservationDate[i].memptr();
            ScalarType * reparam = m_SubjectTimePoints[i].memptr();
            
            for(size_t j = 0; j < N; ++j)
                reparam[j] = AccFactor * (real[j] - TimeShift);
            
            }
    }
     */
}


void
FastNetworkModel
::ComputeDeltas(const Realizations& R) 
{    
    VectorType Delta(m_NbControlPoints);
    Delta(0) = 0.0;
    ScalarType * d = Delta.memptr();

#pragma omp parallel for
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        d[i] = R.at("Delta#" + std::to_string(i), 0);
    }
    
    auto InterpolationCoeff = m_InvertKernelMatrix * Delta;
    m_Deltas = m_InterpolationMatrix * InterpolationCoeff;
}

void
FastNetworkModel
::ComputeNus(const Realizations& R) 
{
    m_Nus = m_InterpolationMatrix * m_InvertKernelMatrix * R.at("Nu");
    /*
    VectorType Nu(m_NbControlPoints);
    ScalarType * n = Nu.memptr();

#pragma omp parallel for
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        n[i] = R.at("Nu", i);
    }
    
    auto InterpolationCoeff = m_InvertKernelMatrix * Nu;
    m_Nus = m_InterpolationMatrix * InterpolationCoeff;
     */
}

void
FastNetworkModel
::ComputeOrthonormalBasis() 
{
    /// Get the data
    auto V0 = m_RandomVariables.GetRandomVariable("Nu")->GetParameter("Mean");
    V0 = exp(V0);
    
    /*
    VectorType U(m_ManifoldDimension);
    ScalarType * n = m_Nus.memptr();
    ScalarType * d = m_Deltas.memptr();
    ScalarType * u = U.memptr();
    
#pragma omp simd
    for(int i = 0; i < m_ManifoldDimension; ++i)
    {
        u[i] = n[i] * V0 / (m_P0 * m_P0)* exp(-d[i]); 
    }
    */
    
    VectorType U = (V0/ (m_P0*m_P0)) * m_Nus % m_Deltas.exp();
    
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
FastNetworkModel
::ComputeAMatrix(const Realizations& R) 
{   
    
    MatrixType NewA(m_ManifoldDimension, m_NbIndependentComponents);
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
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
FastNetworkModel
::ComputeSpaceShifts(const Realizations& R) 
{
    
    MatrixType SS(m_NbIndependentComponents, m_NumberOfSubjects);
    for(int i = 0; i < m_NbIndependentComponents; ++i) 
    {
        
        SS.set_row(i, R.at("S#" + std::to_string(i)));
    }
    m_SpaceShifts = m_AMatrix * SS;
    
}


void
FastNetworkModel
::ComputeBlock1() 
{
    ScalarType * d = m_Deltas.memptr();
    ScalarType *b = m_Block1.memptr();

#pragma omp simd
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        b[i] = 1 / (m_P0 * exp(d[i]));
}


void 
FastNetworkModel
::ComputeBlock2() 
{
    ScalarType * n = m_Nus.memptr();
    ScalarType * b = m_Block2.memptr();

#pragma omp simd
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        b[i] = n[i] / m_P0;
}

FastNetworkModel::VectorType
FastNetworkModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
    double TimePoint = m_SubjectTimePoints[SubjectNumber](ObservationNumber);
    
    VectorType ParallelCurve(m_ManifoldDimension);
    
    
    ScalarType * d = m_Deltas.memptr();
    ScalarType * b1 = m_Block1.memptr();
    ScalarType * n = m_Nus.memptr();
    ScalarType * b2 = m_Block2.memptr();
    ScalarType * w = m_SpaceShifts.get_column(SubjectNumber).memptr();
    ScalarType * p = ParallelCurve.memptr();

//#pragma omp simd
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
    {
        p[i] = m_P0 * exp(d[i] + w[i] * b1[i]  - b2[i] * TimePoint);
        //std::cout << "Observation at : " << m_IndividualObservationDate[SubjectNumber](ObservationNumber) << std::endl;
        //std::cout << p[i] << " = " << m_P0 << " , " << d[i] << " , " << w[i] << " , " << n[i] << " , " << TimePoint << std::endl;
        //int a = 0;
    }
        
    
    //ParallelCurve = m_P0 * (m_Deltas + m_SpaceShifts.get_column(SubjectNumber) % m_Block1 - TimePoint * m_Block2).exp();
    
    return ParallelCurve;
}