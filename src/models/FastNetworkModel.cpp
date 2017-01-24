#include "FastNetworkModel.h"


FastNetworkModel
::FastNetworkModel(const unsigned int NbIndependentComponents,
                          std::shared_ptr<MatrixType> KernelMatrix,
                          std::shared_ptr<MatrixType> InterpolationMatrix) 
{
    m_NbIndependentComponents = NbIndependentComponents;
    m_InvertKernelMatrix = KernelMatrix->transpose();
    m_InterpolationMatrix = *InterpolationMatrix;
    
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_InterpolationCoeffNu.set_size(m_NbControlPoints);
    m_InterpolationCoeffDelta.set_size(m_NbControlPoints);
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
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.0000001 );
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
        auto Nu = std::make_shared<GaussianRandomVariable>(V0, 0.0003*0.0003);
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
    auto Tau = std::make_shared<GaussianRandomVariable>(65, 1.0);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        auto S = std::make_shared<GaussianRandomVariable>(0.0, 1);
        std::string Name = "S#" + std::to_string(i);
        m_IndividualRandomVariables.insert(RandomVariable(Name, S));
    }
    
    /// Other
    double SumObservations = 0.0, K = 0.0;
    for(auto it = D.begin(); it != D.end(); ++it)
    {
        K += it->size();
        for(auto it2 = it->begin(); it2 != it->end(); ++it2)
        {
            SumObservations += it2->first.squared_magnitude();
        }
    }
    m_SumObservations = SumObservations;
    m_NbTotalOfObservations = K;
}


void 
FastNetworkModel
::UpdateModel(const Realizations &R,
              const std::vector<std::string> Names) 
{
    // TODO : it is not the best case : separate nu and delta
    int UpdateCase = 1;
    
    for(auto it = Names.begin(); it != Names.end(); ++it) 
    {
        std::string Name = it->substr(0, it->find_first_of("#"));

        if (Name == "None" or Name == "Ksi" or Name == "Tau") {
            continue;
        } else if (Name == "S") {
            UpdateCase = std::max(UpdateCase, 2);
        } else if (Name == "Beta") {
            UpdateCase = std::max(UpdateCase, 3);
        } else if(Name == "P0") {
            UpdateCase = std::max(UpdateCase, 4);
        }
        else if (Name == "Nu" or Name == "Delta" or Name == "All") {
            UpdateCase = 5;
            break;
        } else {
            UpdateCase = 5;
            std::cout << "Should be" << Name << "be in FastNetworkModel > Update Parameters?" << std::endl;
            break;
        }
    }        
    
    switch(UpdateCase)
    {
        case 1:
            break;
        case 2:
            ComputeSpaceShifts(R);
            break;
        case 3:
            ComputeAMatrix(R);
            ComputeSpaceShifts(R);
            break;
        case 4:
            ComputeOrthonormalBasis(R);
            ComputeAMatrix(R);
            ComputeSpaceShifts(R);
            break;
        case 5:
            ComputeInterpoCoeffDelta(R);
            ComputeInterpoCoeffNu(R);
            ComputeOrthonormalBasis(R);
            ComputeAMatrix(R);
            ComputeSpaceShifts(R);
            break;
        default:
            std::cout << "Error? FastNetworkModel > UpdateModel";
            break;
    }
}

FastNetworkModel::Data
FastNetworkModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    Data D;
    return D;
}


double 
FastNetworkModel
::ComputeLogLikelihood(const Realizations& R, 
                       const Data& D) 
{
    /// Get the data
    return 0;
}

double
FastNetworkModel
::ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, const int SubjectNumber) 
{
    /// Get the data
    double  P0 = exp(R.at("P0")(0));
    VectorType Delta = GetDelta(R);
    VectorType Nu = GetNu(R);
    
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(SubjectNumber));

    auto Block1 = ComputeBlock1(P0, SpaceShift, Delta);
    auto Block2 = ComputeBlock2(P0, Nu);
    
    double LogLikelihood = 0;
    double N = D.at(SubjectNumber).size();
    
#pragma omp parallel for reduction(+:LogLikelihood)   
    for(size_t i = 0; i < N; ++i)
    {
        auto& it = D.at(SubjectNumber).at(i);
        double TimePoint = SubjectTimePoint(it.second);
        VectorType ParallelCurve = ComputeParallelCurve(P0, Block1, Block2, TimePoint);
        LogLikelihood += (it.first - ParallelCurve).squared_magnitude();
    }
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
    
    return LogLikelihood;
}

FastNetworkModel::SufficientStatisticsVector
FastNetworkModel
::GetSufficientStatistics(const Realizations& R, const Data& D) 
{
    
    auto Delta = GetDelta(R);
    auto Nu = GetNu(R);
    double NumberOfSubjects = R.at("Ksi").size();
    
    /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
    VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
    auto itS1 = S1.begin(), itS2 = S2.begin();
    int i = 0;
    for(auto itD = D.begin(); itD != D.end(); ++itD, ++i)
    {
        
        auto TimePoint = GetSubjectTimePoint(i, R);
        auto SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        auto Block1 = ComputeBlock1(exp(R.at("P0")(0)), SpaceShift, Delta);
        auto Block2 = ComputeBlock2(exp(R.at("P0")(0)), Nu);
        
        for(auto itD2 = itD->begin(); itD2 != itD->end(); ++itD2)
        {
            double Time = TimePoint(itD2->second);
            VectorType ParallelVCurve = ComputeParallelCurve(exp(R.at("P0")(0)), Block1, Block2, Time);
            *itS1 = dot_product(ParallelVCurve, itD2->first);
            *itS2 = ParallelVCurve.squared_magnitude();
            ++itS1, ++itS2;
        }
    }
    
    /// S3 <- Ksi_i * Ksi_i
    VectorType S3(NumberOfSubjects);
    auto itKsi = R.at("Ksi").begin();
    for(auto itS3 = S3.begin(); itKsi != R.at("Ksi").end() ; ++itKsi, ++itS3)
    {
        *itS3 = *itKsi * *itKsi;
    }
    
    /// S4 <- Tau_i   &    S5 <- Tau_i * Tau_i
    VectorType S4(NumberOfSubjects), S5(NumberOfSubjects);
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
    double NumberOfSubjects = D.size();
    
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
    KsiVariance /= NumberOfSubjects;
    
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    Ksi->SetVariance(KsiVariance);
    
    /// Update tau and sigma_tau
    double TauMean = 0.0, TauVariance = 0.0;
    for(auto it = SS[3].begin(); it != SS[3].end(); ++it)
    {
        TauMean += *it;
    }
    TauMean /= NumberOfSubjects;
    for(auto it = SS[4].begin(); it != SS[4].end(); ++it)
    {
        TauVariance += *it;
    }
    TauVariance -= NumberOfSubjects * TauMean * TauMean;
    TauVariance /= NumberOfSubjects;
    
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
    std::string FileName = "Outputs/FastNetwork/Parameters" + std::to_string(IterationNumber) + ".txt";
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
    auto SizeW = m_SpaceShifts.at("W0").size();
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        VectorType W = m_SpaceShifts.at("W" + std::to_string(i));
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
    // Todo 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



FastNetworkModel::VectorType
FastNetworkModel
::GetDelta(const Realizations& R) 
{
    return m_InterpolationMatrix * m_InterpolationCoeffDelta;
}


FastNetworkModel::VectorType
FastNetworkModel
::GetNu(const Realizations& R) 
{
    return m_InterpolationMatrix * m_InterpolationCoeffNu;
}


std::function<double(double)> 
FastNetworkModel
::GetSubjectTimePoint(const int SubjectNumber, const Realizations& R) 
{
    double AccFactor = exp(R.at("Ksi")(SubjectNumber));
    double TimeShift = R.at("Tau")(SubjectNumber);
    
    return [AccFactor, TimeShift](double t) { return AccFactor * (t - TimeShift); };
}



void
FastNetworkModel
::ComputeInterpoCoeffDelta(const Realizations& R) 
{
    /*
    VectorType Delta(m_NbControlPoints, 0.0);
    int i = 1;
    for(auto it = Delta.begin() + 1; it != Delta.end(); ++it, ++i)
    {
        *it = R.at("Delta#" + std::to_string(i))(0);
    }
    m_InterpolationCoeffDelta = m_InvertKernelMatrix * Delta;
    */
    
    VectorType Delta(m_NbControlPoints);
    Delta(0) = 0.0;
    ScalarType * d = Delta.memptr();

#pragma omp parallel for
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        d[i] = R.at("Delta#" + std::to_string(i))(0);
    }
    
    m_InterpolationCoeffDelta = m_InvertKernelMatrix * Delta;
}

void
FastNetworkModel
::ComputeInterpoCoeffNu(const Realizations& R) 
{
    /*
    VectorType Nu(m_NbControlPoints, 1.0);
    int i = 1;
    for(auto it = Nu.begin() + 1; it != Nu.end(); ++it, ++i)
    {
        *it = R.at("Nu#" + std::to_string(i))(0);
    }
    m_InterpolationCoeffNu = m_InvertKernelMatrix * Nu;
    */
    
    VectorType Nu(m_NbControlPoints);
    ScalarType * n = Nu.memptr();

#pragma omp parallel for
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        n[i] = R.at("Nu#" + std::to_string(i))(0);
    }
    
    m_InterpolationCoeffNu = m_InvertKernelMatrix * Nu;
}


void
FastNetworkModel
::ComputeOrthonormalBasis(const Realizations& R) 
{
    /// Get the data
    auto P0 = exp(R.at("P0")(0));
    auto Nu = GetNu(R);
    auto Delta = GetDelta(R);
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    auto V0 = exp(Ksi->GetMean());
    
    VectorType U(m_ManifoldDimension);
    
    ScalarType * n = Nu.memptr();
    ScalarType * d = Delta.memptr();
    ScalarType * u = U.memptr();
    
#pragma omp simd
    for(int i = 0; i < m_ManifoldDimension; ++i)
    {
        u[i] = n[i] * V0 / (P0 * P0)* exp(-d[i]); 
    }
    
    /// Compute the initial pivot vector U
    double Norm = U.magnitude();
    U(0) += copysign(1, -U(0)) * Norm;
    double NormU = U.squared_magnitude();
    
        /// Compute Q = I(N) - 2U . Ut / (Ut . U)
    std::vector<VectorType> Q;
    for (auto it = U.begin(); it != U.end(); ++it) {
        VectorType Coordinate(U.size());
        int i = 0;
        for (auto it2 = U.begin(); it2 != U.end(); ++it2, ++i) 
        {
            Coordinate(i) = -2 * *it * *it2 / NormU;
        }
        Q.push_back(Coordinate);
    }

    for (int i = 0; i < Q.size(); ++i) {
        Q[i](i) += 1;
    }
    
    /// TESTS NEEDED
    // TODO : Mettre des tests
    /// END OF TESTS
    
    Q.erase(Q.begin());
    m_OrthogonalBasis = std::move(Q);
    
}

void 
FastNetworkModel
::ComputeAMatrix(const Realizations& R) 
{
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
    
    /// TESTS NEEDED
    // TODO : Mettre des tests
    /// END TESTS
}

void 
FastNetworkModel
::ComputeSpaceShifts(const Realizations& R) 
{
    std::map< std::string, VectorType> SpaceShifts;
    int NumberOfSubjects = (int)R.at("Tau").size();
    TestAssert::WarningEquality_Object(NumberOfSubjects, (int)R.at("Ksi").size(), "FastNetwork > ComputeSpaceShifts");

    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        VectorType Si(m_NbIndependentComponents);
        for(int j = 0; j < m_NbIndependentComponents; ++j)
        {
            Si(j) = R.at("S#" + std::to_string(j))(i);
        }

        VectorType V = m_AMatrix * Si;
        std::pair< std::string, VectorType > SpaceShift("W"+std::to_string(i),  V);
        SpaceShifts.insert( SpaceShift);
    }

    m_SpaceShifts = std::move(SpaceShifts);
    
    /// TESTS NEEDED
    // TODO : Mettre des tests
    /// END TESTS
}


FastNetworkModel::VectorType
FastNetworkModel
::ComputeParallelCurve(double P0, VectorType &Block1, VectorType &Block2, double Timepoint) 
{
    auto N = Block1.size();
    VectorType ParallelCurve(N);
    
    ScalarType * curve = ParallelCurve.memptr();
    ScalarType * b1 = Block1.memptr();
    ScalarType * b2 = Block2.memptr();
    
#pragma omp simd
    for(size_t i = 0; i < N; ++i)
        curve[i] = P0 * exp(b1[i] - b2[i] * Timepoint);
    
    return ParallelCurve;
}


FastNetworkModel::VectorType
FastNetworkModel
::ComputeBlock1(double P0, VectorType& SpaceShift, VectorType& Delta)
{
    auto N = SpaceShift.size();
    VectorType Block1(N);
    
    ScalarType * s = SpaceShift.memptr();
    ScalarType * d = Delta.memptr();
    ScalarType * b = Block1.memptr();

#pragma omp simd
    for(size_t i = 0; i < N; ++i)
        b[i] = s[i] / (P0 * exp(d[i])) + d[i];
    
    return Block1;
}

FastNetworkModel::VectorType
FastNetworkModel
::ComputeBlock2(double P0, VectorType& Nu)
{
    auto N = Nu.size();
    VectorType Block2(N);
    
    ScalarType * n = Nu.memptr();
    ScalarType * b = Block2.memptr();

#pragma omp simd
    for(size_t i = 0; i < N; ++i) 
        b[i] = n[i] / P0;

    return Block2;
}