#include "NetworkPropagationModel.h"

NetworkPropagationModel
::NetworkPropagationModel(const unsigned int NbIndependentComponents,
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


NetworkPropagationModel
::~NetworkPropagationModel() 
{
    
   
    
}

void
NetworkPropagationModel
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
        auto Nu = std::make_shared<GaussianRandomVariable>(1, 0.003*0.003);
        
        std::string Name1 = "Delta#" + std::to_string(i);
        std::string Name2 = "Nu#" + std::to_string(i);
        
        m_PopulationRandomVariables.insert(RandomVariable(Name1, Delta));
        m_PopulationRandomVariables.insert(RandomVariable(Name2, Nu));
        
    }
    
    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension - 1); ++i)
    {
        auto Beta = std::make_shared<GaussianRandomVariable>(0, 0.0005);
        std::string Name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Beta));
    }
    
    
    /// Individual variables
    auto Ksi = std::make_shared<GaussianRandomVariable>(-3.06, 0.4);
    auto Tau = std::make_shared<GaussianRandomVariable>(65, 1.0);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        auto S = std::make_shared<GaussianRandomVariable>(0.0, 0.5);
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
NetworkPropagationModel
::UpdateModel(const Realizations &R,
              const std::vector<std::string> Names) 
{
    // TODO : it is not the best case : separate rho, nu and delta
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
            std::cout << "Should be" << Name << "be in NetworkPropagationModel > Update Parameters?" << std::endl;
            break;
        }
    }
    
    auto R1 = std::make_shared<Realizations>(R);
        
    
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
            std::cout << "Error? NetworkPropagationModel > UpdateModel";
            break;
    }
}

NetworkPropagationModel::Data
NetworkPropagationModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    // Osef
    Data D;
    return D;
}


double 
NetworkPropagationModel
::ComputeLogLikelihood(const Realizations& R, 
                       const Data& D) 
{
    /// Get the data
    return 0;
}

double
NetworkPropagationModel
::ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, const int SubjectNumber) 
{
    /// Get the data
    VectorType Rho(1, exp(R.at("P0")(0)));
    VectorType Delta = GetDelta(R);
    VectorType Nu = GetNu(R);
    
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(SubjectNumber));

    double LogLikelihood = 0;
    double N = D.at(SubjectNumber).size();
    
#pragma omp parallel for reduction(+:LogLikelihood)   
    for(size_t i = 0; i < N; ++i)
    {
        auto& it = D.at(SubjectNumber).at(i);
        double TimePoint = SubjectTimePoint(it.second);
        VectorType ParallelCurve = ComputeParallelCurve(Rho, Delta, Nu, SpaceShift, TimePoint);
        LogLikelihood += (it.first - ParallelCurve).squared_magnitude();
    }
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
    
    return LogLikelihood;
}

NetworkPropagationModel::SufficientStatisticsVector
NetworkPropagationModel
::GetSufficientStatistics(const Realizations& R,
                          const Data& D) 
{
    /// 
    VectorType P0(1, exp(R.at("P0")(0)));
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
        for(auto itD2 = itD->begin(); itD2 != itD->end(); ++itD2)
        {
            double Time = TimePoint(itD2->second);
            VectorType ParallelVCurve = ComputeParallelCurve(P0, Delta, Nu, SpaceShift, Time);
            *itS1 = dot_product(ParallelVCurve, itD2->first);
            *itS2 = ParallelVCurve.squared_magnitude();
            ++itS1, ++itS2;
        }
    }
    
    /// S3 <- Ksi_i   &    S4 <- Ksi_i * Ksi_i
    VectorType S3(NumberOfSubjects), S4(NumberOfSubjects);
    auto itKsi = R.at("Ksi").begin();
    auto itS3 = S3.begin(), itS4 = S4.begin();
    for(    ; itKsi != R.at("Ksi").end() ; ++itKsi, ++itS3, ++itS4)
    {
        *itS3 = *itKsi;
        *itS4 = *itKsi * *itKsi;
    }
    
    /// S5 <- Tau_i   &    S6 <- Tau_i * Tau_i
    VectorType S5(NumberOfSubjects), S6(NumberOfSubjects);
    auto itTau = R.at("Tau").begin();
    auto itS5 = S5.begin(), itS6 = S6.begin();
    for(    ; itTau != R.at("Tau").end(); ++itTau, ++itS5, ++itS6)
    {
        *itS5 = *itTau;
        *itS6 = *itTau * *itTau;
    }
    
    /// S7 <- P0 
    VectorType S7(1, R.at("P0")(0));
    
    /// S8 <- beta_k
    VectorType S8((m_ManifoldDimension-1) * m_NbIndependentComponents);
    i = 0;
    for(auto it = S8.begin(); it != S8.end(); ++it, ++i)
    {
        *it = R.at("Beta#" + std::to_string(i))(0);
    }
    
    /// S9 <- rho_k, S10 <- nu_k
    VectorType S9(m_NbControlPoints - 1), S10(m_NbControlPoints - 1);
    i = 1;
    auto itS9 = S9.begin(), itS10 = S10.begin();
    for(    ; itS9 != S9.end(); ++itS9, ++itS10, ++i)
    {
        *itS9 = R.at("Delta#" + std::to_string(i))(0);
        *itS10 = R.at("Nu#" + std::to_string(i))(0);
    }
    
    SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10};
    return S;
}   

void
NetworkPropagationModel
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
    double KsiMean = 0.0, KsiVariance = 0.0;
    for(auto it = SS[2].begin(); it != SS[2].end(); ++it)
    {
        KsiMean += *it;
    }
    KsiMean /= NumberOfSubjects;
    for(auto it = SS[3].begin(); it != SS[3].end(); ++it)
    {
        KsiVariance += *it;
    }
    KsiVariance -= NumberOfSubjects * KsiMean * KsiMean;
    KsiVariance /= NumberOfSubjects;
    
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    Ksi->SetMean(KsiMean);
    Ksi->SetVariance(KsiVariance);
    
    /// Update tau and sigma_tau
    double TauMean = 0.0, TauVariance = 0.0;
    for(auto it = SS[4].begin(); it != SS[4].end(); ++it)
    {
        TauMean += *it;
    }
    TauMean /= NumberOfSubjects;
    for(auto it = SS[5].begin(); it != SS[5].end(); ++it)
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
    P0->SetMean(SS[6](0));
    
    /// Update beta_k
    int i = 0;
    for(auto it = SS[7].begin(); it != SS[7].end(); ++it, ++i)
    {
        auto AbstractBeta = m_PopulationRandomVariables.at("Beta#" + std::to_string(i));
        auto Beta = std::static_pointer_cast<GaussianRandomVariable>(AbstractBeta);
        Beta->SetMean(*it);
    }
        
    /// Update delta_k
    i = 1;
    for(auto it = SS[8].begin(); it != SS[8].end(); ++it, ++i)
    {
        auto AbstractDelta = m_PopulationRandomVariables.at("Delta#" + std::to_string(i));
        auto Delta = std::static_pointer_cast<GaussianRandomVariable>(AbstractDelta);
        Delta->SetMean(*it);
    }
    
    /// Update nu_k
    i = 1;
    for(auto it = SS[9].begin(); it != SS[9].end(); ++it, ++i)
    {
        auto AbstractNu = m_PopulationRandomVariables.at("Nu#" + std::to_string(i));
        auto Nu = std::static_pointer_cast<GaussianRandomVariable>(AbstractNu);
        Nu->SetMean(*it);
    }
}

void 
NetworkPropagationModel
::DisplayOutputs() 
{
    double Sigma = m_Noise->GetVariance();
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    auto Delta1 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#1"));
    auto Delta10 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#10"));
    auto Delta100 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#100"));
    auto Nu1 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Nu#1"));
    auto Nu10 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Nu#10"));
    auto Nu100 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Nu#100"));
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));   
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    
    std::cout << "P1 : " << exp(P0->GetMean()) << ". T0 : " << Tau->GetMean() << ". Var(Tau) : " << Tau->GetVariance();
    std::cout << ". V0 : " << exp(Ksi->GetMean()) << ". Var(Ksi) : " << Ksi->GetVariance();
    std::cout << ". Delta1 : " << Delta1->GetMean() << ". Delta10 : " << Delta10->GetMean() << ". Delta100 : " << Delta100->GetMean();
    std::cout << ". Nu1 : " << Nu1->GetMean() << ". Nu10 : " << Nu10->GetMean() << ". Nu100 : " << Nu100->GetMean();
    std::cout << ". Sigma : " << Sigma << std::endl;
}




void 
NetworkPropagationModel
::SaveData(unsigned int IterationNumber, const Realizations& R) 
{
    std::ofstream Outputs;
    Outputs.open("MultiSignalNetwork_Parameters.txt", std::ofstream::out | std::ofstream::trunc);
    
    Outputs << m_NbControlPoints << std::endl;
    
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    
    Outputs << exp(P0->GetMean()) << std::endl;
    Outputs << Tau->GetMean()  << std::endl;
    Outputs << exp(Ksi->GetMean()) << std::endl;
    
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        auto AbstractDelta = m_PopulationRandomVariables.at("Delta#" + std::to_string(i));
        auto Delta = std::static_pointer_cast<GaussianRandomVariable>(AbstractDelta);
        Outputs << Delta->GetMean() << std::endl;
    }
    
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        auto AbstractNu = m_PopulationRandomVariables.at("Nu#" + std::to_string(i));
        auto Nu = std::static_pointer_cast<GaussianRandomVariable>(AbstractNu);
        Outputs << Nu->GetMean() << std::endl;
    }
}


void 
NetworkPropagationModel
::InitializeFakeRandomVariables() 
{
    // Todo 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



NetworkPropagationModel::VectorType
NetworkPropagationModel
::GetDelta(const Realizations& R) 
{
    return m_InterpolationMatrix * m_InterpolationCoeffDelta;
}


NetworkPropagationModel::VectorType
NetworkPropagationModel
::GetNu(const Realizations& R) 
{
    return m_InterpolationMatrix * m_InterpolationCoeffNu;
}


std::function<double(double)> 
NetworkPropagationModel
::GetSubjectTimePoint(const int SubjectNumber, const Realizations& R) 
{
    double AccFactor = exp(R.at("Ksi")(SubjectNumber));
    double TimeShift = R.at("Tau")(SubjectNumber);
    
    return [AccFactor, TimeShift](double t) { return AccFactor * (t - TimeShift); };
}



void
NetworkPropagationModel
::ComputeInterpoCoeffDelta(const Realizations& R) 
{
    VectorType Delta(m_NbControlPoints, 0.0);
    int i = 1;
    for(auto it = Delta.begin() + 1; it != Delta.end(); ++it)
    {
        *it = R.at("Delta#" + std::to_string(i))(0);
    }
    m_InterpolationCoeffDelta = m_InvertKernelMatrix * Delta;
}

void
NetworkPropagationModel
::ComputeInterpoCoeffNu(const Realizations& R) 
{
    VectorType Nu(m_NbControlPoints, 1.0);
    int i = 1;
    for(auto it = Nu.begin() + 1; it != Nu.end(); ++it)
    {
        *it = R.at("Nu#" + std::to_string(i))(0);
    }
    m_InterpolationCoeffNu = m_InvertKernelMatrix * Nu;
}


void
NetworkPropagationModel
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
        for (auto it2 = U.begin(); it2 != U.end(); ++it2, ++i) {
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
NetworkPropagationModel
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
NetworkPropagationModel
::ComputeSpaceShifts(const Realizations& R) 
{
    std::map< std::string, VectorType> SpaceShifts;
    int NumberOfSubjects = (int)R.at("Tau").size();
    TestAssert::WarningEquality_Object(NumberOfSubjects, (int)R.at("Ksi").size(), "NetworkPropagation > ComputeSpaceShifts");

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


NetworkPropagationModel::VectorType
NetworkPropagationModel
::ComputeParallelCurve(VectorType &P0, VectorType &Delta, VectorType &Nu, VectorType& SpaceShift, double Timepoint) 
{
    auto N = Delta.size();
    double Position = P0(0);
    VectorType ParallelCurve(N);
    
    ScalarType * curve = ParallelCurve.memptr();
    ScalarType * d = Delta.memptr();
    ScalarType * n = Nu.memptr();
    ScalarType * s = SpaceShift.memptr();
    
#pragma omp simd
    for(size_t i = 0; i < N; ++i)
        curve[i] = Position * exp(s[i] / (Position * exp(d[i])) + d[i] - n[i] * Timepoint / Position);
    
    return ParallelCurve;
}