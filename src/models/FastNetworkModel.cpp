#include "FastNetworkModel.h"
#include <armadillo>


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
::Initialize(const Data& D) 
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
     /// Population variables
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.000001 );
    auto P0 = std::make_shared<GaussianRandomVariable>(0.1, 0.0001 * 0.0001);
    m_PopulationRandomVariables.insert(RandomVariable("P0", P0));
        
    for(int i = 1; i < m_NbControlPoints; ++i)
    {
        auto Delta = std::make_shared<GaussianRandomVariable>(0, 0.003*0.003);
        std::string Name1 = "Delta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name1, Delta));
    }
    
    // TODO : VERY IMPORTANT : This can be refactored such that there is only one random variable
    //                         with Multiple Realizations. It would be m_NbControlPoints realizations here
    double V0 = 0.04088;
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        auto Nu = std::make_shared<GaussianRandomVariable>(V0, 0.0004*0.0004);
        std::string Name2 = "Nu#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name2, Nu));
    }
    
    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
    {
        auto Beta = std::make_shared<GaussianRandomVariable>(0, 0.001);
        std::string Name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Beta));
    }
    
    
    /// Individual variables
    auto Ksi = std::make_shared<GaussianRandomVariable>(0, 0.000000004);
    auto Tau = std::make_shared<GaussianRandomVariable>(62, 0.25);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        auto S = std::make_shared<GaussianRandomVariable>(0.0, 1);
        std::string Name = "S#" + std::to_string(i);
        m_IndividualRandomVariables.insert(RandomVariable(Name, S));
    }
    
    /// Other
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
    
}


void 
FastNetworkModel
::UpdateModel(const Realizations &R, int Type,
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
        else if(Name == "P0")
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
            ComputeSubjectTimePoint(R, -1);
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
    if(IndividualOnly) ComputeSubjectTimePoint(R, Type);
    
    if(ComputePosition) { m_P0 = exp(R.at("P0")(0)); }
    if(ComputeDelta) ComputeDeltas(R);
    if(ComputeNu) ComputeNus(R);
    if(ComputeBasis) ComputeOrthonormalBasis();
    if(ComputeA) ComputeAMatrix(R);
    if(ComputeSpaceShift) ComputeSpaceShifts(R);
    if(ComputeBlock_1) ComputeBlock1();
    if(ComputeBlock_2) ComputeBlock2();
    
}

FastNetworkModel::Data
FastNetworkModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    m_NumberOfSubjects = NumberOfSubjects;
    auto R = SimulateRealizations(NumberOfSubjects);
    m_P0 = exp(R.at("P0")(0));
    ComputeDeltas(R);
    ComputeNus(R);
    ComputeOrthonormalBasis();
    ComputeAMatrix(R);
    ComputeSpaceShifts(R);
    ComputeBlock1();
    ComputeBlock2();
    
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::uniform_real_distribution<double> ObsDistrib(60, 95);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    Data D;
    double RealNoise = 0.0;
    
    for(int i = 0; i < NumberOfSubjects; ++i)
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
    
    std::cout << "Real Likelihood = " << ComputeLogLikelihood(R, D) << std::endl;
    
    return D;
}


double 
FastNetworkModel
::ComputeLogLikelihood(const Realizations& R, const Data& D) 
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
::ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, const int SubjectNumber) 
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
::GetSufficientStatistics(const Realizations& R, const Data& D) 
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
    VectorType S3(m_NumberOfSubjects);
    auto itKsi = R.at("Ksi").begin();
    for(auto itS3 = S3.begin(); itKsi != R.at("Ksi").end() ; ++itKsi, ++itS3)
    {
        *itS3 = *itKsi * *itKsi;
    }
    
    /// S4 <- Tau_i   &    S5 <- Tau_i * Tau_i
    VectorType S4(m_NumberOfSubjects), S5(m_NumberOfSubjects);
    auto itTau = R.at("Tau").begin();
    auto itS4 = S4.begin(), itS5 = S5.begin();
    for(    ; itTau != R.at("Tau").end(); ++itTau, ++itS4, ++itS5)
    {
        *itS4 = *itTau;
        *itS5 = *itTau * *itTau;
    }
    
    /// S6 <- P0 
    VectorType S6(1, R.at("P0")(0));
    
    /// S7 <- beta_k
    VectorType S7((m_ManifoldDimension-1) * m_NbIndependentComponents);
    i = 0;
    for(auto it = S7.begin(); it != S7.end(); ++it, ++i)
    {
        *it = R.at("Beta#" + std::to_string(i))(0);
    }
    
    /// S8 <- delta_k
    VectorType S8(m_NbControlPoints - 1);
    i = 1;
    
    for( auto itS8 = S8.begin(); itS8 != S8.end(); ++itS8, ++i)
    {
        *itS8 = R.at("Delta#" + std::to_string(i))(0);
        
    }
    
    // S9 <- nu_k, S10 = nu_k * nu_k
    VectorType S9(m_NbControlPoints), S10(m_NbControlPoints);
    i = 0;
    auto itS9 = S9.begin(), itS10 = S10.begin();
    for( ; itS9 != S9.end(); ++itS9, ++itS10, ++i)
    {
        double Nu_ = R.at("Nu#" + std::to_string(i))(0);
        *itS9 = Nu_;
        *itS10 = Nu_*Nu_;   
    }
    
    
    SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10};
    return S;
}   

void
FastNetworkModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS, const Data& D) 
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
    
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    Ksi->SetVariance(KsiVariance);
    
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
    
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    Tau->SetMean(TauMean);
    Tau->SetVariance(TauVariance);
    
    /// Update P0
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    P0->SetMean(SS[5](0));
    
    /// Update beta_k
    int i = 0;
    for(auto it = SS[6].begin(); it != SS[6].end(); ++it, ++i)
    {
        auto AbstractBeta = m_PopulationRandomVariables.at("Beta#" + std::to_string(i));
        auto Beta = std::static_pointer_cast<GaussianRandomVariable>(AbstractBeta);
        Beta->SetMean(*it);
    }
        
    /// Update delta_k
    i = 1;
    for(auto it = SS[7].begin(); it != SS[7].end(); ++it, ++i)
    {
        auto AbstractDelta = m_PopulationRandomVariables.at("Delta#" + std::to_string(i));
        auto Delta = std::static_pointer_cast<GaussianRandomVariable>(AbstractDelta);
        Delta->SetMean(*it);
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
    
    for(i = 0; i < m_NbControlPoints; ++i)
    {
        auto AbstractNu = m_PopulationRandomVariables.at("Nu#" + std::to_string(i));
        auto Nu = std::static_pointer_cast<GaussianRandomVariable>(AbstractNu);
        Nu->SetMean(V0);
        Nu->SetVariance(NuVariance);
    }
}

void 
FastNetworkModel
::DisplayOutputs(const Realizations& R) 
{
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));   
    
    auto Nu = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Nu#0"));
    double NuMax = R.at("Nu#0")(0);
    double NuMin = NuMax;
    for(size_t i = 1; i < 258; ++i)
    {
        double NuK = R.at("Nu#" + std::to_string(i))(0);
        NuMax = std::max(NuMax, NuK);
        NuMin = std::min(NuMin, NuK);
    }
    
    double DeltaMin = R.at("Delta#1")(0);
    double DeltaMax = DeltaMin;
    for(size_t i = 1; i < 258; ++i)
    {
        double DeltaK = R.at("Delta#" + std::to_string(i))(0);
        DeltaMax = std::max(DeltaMax, DeltaK);
        DeltaMin = std::min(DeltaMin, DeltaK);
    }
    
    
    std::cout << "Noise: " << m_Noise->GetVariance();
    std::cout << " - p0: " << exp(P0->GetMean()) << " - t0: " << Tau->GetMean() << " - S_tau:" << Tau->GetVariance();
    std::cout << " - v0: " << Nu->GetMean() << " - S_nu:" << Nu->GetVariance();
    std::cout << " - MaxNu: " << NuMax << " - MinNu: " << NuMin;
    std::cout << " - MaxDelta: " << DeltaMax << " - MinDelta: " << DeltaMin << std::endl;
}


void 
FastNetworkModel
::SaveData(unsigned int IterationNumber, const Realizations& R) 
{
    unsigned int NumberOfSubjects = R.at("Tau").size();
    std::ofstream Outputs;    
    std::string FileName = "/Users/igor.koval/Documents/Work/RiemAlzh/src/io/outputs/Normalized/Parameters" + std::to_string(IterationNumber) + ".txt";
    //std::string FileName = "Outputs/FastNetwork/Parameters.txt";
    Outputs.open(FileName, std::ofstream::out | std::ofstream::trunc);
    
    
    /// Save the final noise variance
    Outputs << m_Noise->GetVariance() << std::endl;
    
   /// Save Number of subject, Dimensions, Number of Sources, Number of control points
    Outputs << NumberOfSubjects << ", " << m_ManifoldDimension << ", " << m_NbIndependentComponents  << ", " << m_NbControlPoints << std::endl;
    
    /// Save P0, P0_mean and P0_var
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    Outputs << R.at("P0")(0) << ", " << P0->GetMean() << ", " << P0->GetVariance() << std::endl;
    
    /// Save (Ksi_i)
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        Outputs << R.at("Ksi")(i);
        if(i != NumberOfSubjects - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save Ksi_mean and Ksi_Var
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    Outputs << Ksi->GetMean() << ", " << Ksi->GetVariance() << std::endl;
    
    /// Save V0 
    auto Nu = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Nu#0"));
    Outputs << Nu->GetMean() << std::endl;
    
    /// Save (Tau_i)
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        Outputs << R.at("Tau")(i);
        if(i != NumberOfSubjects - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save Tau_Mean and Tau_Var
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    Outputs << Tau->GetMean() << ", " << Tau->GetVariance() << std::endl;
    
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
        auto Delta = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at(Name));
        Outputs << Delta->GetMean();
        if(i != m_NbControlPoints - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Nu_k)
    for(size_t i = 0; i < m_NbControlPoints; ++i)  
    {
        Outputs << R.at("Nu#" + std::to_string(i))(0);  
        if(i != m_NbControlPoints - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Nu_k_mean)
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        std::string Name = "Nu#" + std::to_string(i);
        auto Nu = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at(Name));
        Outputs << Nu->GetMean();
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
        //VectorType W = m_SpaceShifts[i];
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
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.0000000001 );
    auto P0 = std::make_shared<GaussianRandomVariable>(0.932, 0.0001 * 0.0001);
    m_PopulationRandomVariables.insert(RandomVariable("P0", P0));
        
    for(int i = 1; i < m_NbControlPoints; ++i)
    {
        auto Delta = std::make_shared<GaussianRandomVariable>(0, 0.0008*0.0008);
        std::string Name1 = "Delta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name1, Delta));
    }
    
    // TODO : VERY IMPORTANT : This can be refactored such that there is only one random variable
    //                         with Multiple Realizations. It would be m_NbControlPoints realizations here
    double V0 = 0.04688;
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        auto Nu = std::make_shared<GaussianRandomVariable>(V0, 0.0006*0.0006);
        std::string Name2 = "Nu#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name2, Nu));
    }
    
    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
    {
        auto Beta = std::make_shared<GaussianRandomVariable>(0, 0.0005);
        std::string Name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Beta));
    }
    
    
    /// Individual variables
    auto Ksi = std::make_shared<GaussianRandomVariable>(0, 0.000000004);
    auto Tau = std::make_shared<GaussianRandomVariable>(65, 0.25);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        auto S = std::make_shared<GaussianRandomVariable>(0.0, 0.5);
        std::string Name = "S#" + std::to_string(i);
        m_IndividualRandomVariables.insert(RandomVariable(Name, S));
    }
    
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
FastNetworkModel
::ComputeSubjectTimePoint(const Realizations &R, const int SubjectNumber) 
{
    if(SubjectNumber != -1) {
        double AccFactor = exp(R.at("Ksi")(SubjectNumber));
        double TimeShift = R.at("Tau")(SubjectNumber);
        
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
        d[i] = R.at("Delta#" + std::to_string(i))(0);
    }
    
    auto InterpolationCoeff = m_InvertKernelMatrix * Delta;
    m_Deltas = m_InterpolationMatrix * InterpolationCoeff;
}

void
FastNetworkModel
::ComputeNus(const Realizations& R) 
{
    VectorType Nu(m_NbControlPoints);
    ScalarType * n = Nu.memptr();

#pragma omp parallel for
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        n[i] = R.at("Nu#" + std::to_string(i))(0);
    }
    
    auto InterpolationCoeff = m_InvertKernelMatrix * Nu;
    m_Nus = m_InterpolationMatrix * InterpolationCoeff;
}

void
FastNetworkModel
::ComputeOrthonormalBasis() 
{
    /// Get the data
    auto V0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    V0 = exp(V0);
    
    VectorType U(m_ManifoldDimension);
    ScalarType * n = m_Nus.memptr();
    ScalarType * d = m_Deltas.memptr();
    ScalarType * u = U.memptr();
    
#pragma omp simd
    for(int i = 0; i < m_ManifoldDimension; ++i)
    {
        u[i] = n[i] * V0 / (m_P0 * m_P0)* exp(-d[i]); 
    }
    
    /// Compute the initial pivot vector U
    double Norm = U.magnitude();
    U(0) += copysign(1, -U(0)) * Norm;
        
    // Check the vectorise_row_wise function of the ArmadilloMatrixWrapper
    //MatrixType Id = diagonal_matrix(m_ManifoldDimension, 1.0);
    //MatrixType Id = identity_matrix(m_ManifoldDimension);
    MatrixType U2(U);
    double NormU2 = U.squared_magnitude();
    //MatrixType FinalMatrix2 = Id - 2.0/NormU2 * U2*U2.transpose();
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
            Beta(j) = R.at( "Beta#" + Number)(0);
        }
        
        NewA.set_column(i, m_OrthogonalBasis * Beta);
    }
    
    m_AMatrix = NewA;
     
    /*
    MatrixType AMatrix(m_ManifoldDimension, m_NbIndependentComponents);
    
    for(int i = 0; i < m_NbIndependentComponents ; ++i)
    {
        VectorType Beta(m_ManifoldDimension - 1);
        for(int j = 0; j < m_ManifoldDimension - 1 ; ++j)
        {
            std::string Number = std::to_string(int(j + i*(m_ManifoldDimension - 1)));
            Beta(j) = R.at( "Beta#" + Number)(0);
        }
        
        VectorType V = LinearCombination(Beta, m_OrthogonalBasis)   ;     
        AMatrix.set_column(i, V);
    }
    
    m_AMatrix = std::move(AMatrix);
    */
}

void 
FastNetworkModel
::ComputeSpaceShifts(const Realizations& R) 
{
    
    
    MatrixType SS(m_NbIndependentComponents, m_NumberOfSubjects);
    for(int i = 0; i < m_NumberOfSubjects; ++i) 
    {
        for (int j = 0; j < m_NbIndependentComponents; ++j) 
        {
            SS(j, i) = R.at("S#" + std::to_string(j))(i);
        }
    }
    m_SpaceShifts = m_AMatrix * SS;
    
    /*
    std::vector<VectorType> SpaceShifts(m_NumberOfSubjects);
    for(int i = 0; i < m_NumberOfSubjects; ++i)
    {
        VectorType Si(m_NbIndependentComponents);
        for(int j = 0; j < m_NbIndependentComponents; ++j)
        {
            Si(j) = R.at("S#" + std::to_string(j))(i);
        }
        
        SpaceShifts[i] = m_AMatrix * Si;
    }
    
    m_SpaceShifts = std::move(SpaceShifts);
    */
}


void
FastNetworkModel
::ComputeBlock1() 
{
    ScalarType * d = m_Deltas.memptr();
    ScalarType *b = m_Block1.memptr();

#pragma omp simd
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        b[i] = m_P0 * exp(d[i]);
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
    ScalarType * b2 = m_Block2.memptr();
    ScalarType * w = m_SpaceShifts.get_column(SubjectNumber).memptr();
    //ScalarType * w = m_SpaceShifts[SubjectNumber].memptr();
    ScalarType * p = ParallelCurve.memptr();

//#pragma omp simd
    for(size_t i = 0; i < m_ManifoldDimension; ++i)
        p[i] = m_P0 * exp(d[i] + w[i] / b1[i]  - b2[i] * TimePoint);
    
    return ParallelCurve;
}