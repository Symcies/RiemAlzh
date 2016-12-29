#include "NetworkPropagationModel2.h"
#include "../Manifolds/ExponentialCurveManifold.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////
    
NetworkPropagationModel2
::NetworkPropagationModel2(const unsigned int NbIndependentComponents,
                           std::shared_ptr<AbstractManifold> M,
                           std::shared_ptr<MatrixType> KernelMatrix,
                           std::shared_ptr<MatrixType> InterpolationMatrix) 
{
    m_NbIndependentComponents = NbIndependentComponents;
    m_Manifold = M;
    m_InvertKernelMatrix = *KernelMatrix;
    m_InterpolationMatrix = *InterpolationMatrix;
    m_OutputParameters.open("ParametersMciConverters_Linear.txt", std::ofstream::out | std::ofstream::trunc);
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_InterpolationCoefficients.set_size(m_NbControlPoints);
    
    //// Tests
    TestAssert::WarningEquality_Object(m_Manifold->GetDimension(), (double)m_InterpolationMatrix.rows(), "Error Network Propagation Constructor > 1");
    TestAssert::WarningEquality_Object(m_NbControlPoints, m_InterpolationMatrix.columns(), "Error Network Propagation Constructor > 2");
    TestAssert::WarningEquality_Object(m_NbControlPoints, m_InvertKernelMatrix.rows(), "Error Network Propagation Constructor > 3");
    TestAssert::WarningEquality_Object(m_NbControlPoints, m_InvertKernelMatrix.columns(), "Network Propagation Constructor > 4");
    
}

NetworkPropagationModel2
::~NetworkPropagationModel2() 
{
    
}


void
NetworkPropagationModel2
::Initialize(const Data& D) 
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    m_IndividualRandomVariables.clear();
    m_PopulationRandomVariables.clear();
    
    
    /////////////////////////////
    /// Population Parameters ///
    /////////////////////////////
    auto P0 = std::make_shared<GaussianRandomVariable>(0.932, 0.00000001);
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0) );
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.0000001 );
    
    //////////////////////////////////
    ///   Read the initial delta   ///
    //////////////////////////////////
    /*
    std::ifstream DeltaFile ("/Users/igor.koval/Documents/Work/RiemAlzh/datatest/delta_init.csv");
    if(DeltaFile.is_open())
    {
        unsigned int i = 0;
        std::string line;
        while(getline(DeltaFile, line))
        {
            if(i != 0) {
                double DeltaMean = std::stod(line);
                auto Delta = std::make_shared<GaussianRandomVariable>(DeltaMean, 0.0001);
                std::string Name = "Delta#" + std::to_string(i);
                m_PopulationRandomVariables.insert(RandomVariable(Name, Delta));
            }
            ++i;
        }
        TestAssert::WarningEquality_Object(i, m_NbControlPoints, "Network Propagation Model > Initialize - deltas");
    }
    else { std::cout << "Unable to open the initial delta"; }
    */
    

     
    
    /////////////////////////////
    /// Individual Parameters ///
    /////////////////////////////
    
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
    
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension() - 1); ++i)
    {
        auto Beta = std::make_shared<GaussianRandomVariable>(0, 0.0005);
        std::string Name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Beta));
    }
    
    for(int i = 1; i < m_NbControlPoints; ++i)
    {
        auto Delta = std::make_shared<GaussianRandomVariable>(0, 0.0008*0.0008);
        std::string Name = "Delta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Delta));
    }
    
    ///////////////////////////////////
    ///   Sufficient Statistic S0   ///
    /// Corresponding to Sum(y_ijk) ///
    ///////////////////////////////////
    
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
    
    
    /////////////////////////////
    ///   Read Initialisation ///
    /// From previous run     ///
    /////////////////////////////
    /*
    std::ifstream Parameters("/Users/igor.koval/Documents/Work/RiemAlzh/ModelParametersMCIConverters.txt");
    if(Parameters.is_open())
    {
        std::string line;
        getline(Parameters, line);
        while(getline(Parameters, line))
        {
            std::stringstream LineStream(line);
            std::string cell;
            bool Pop = false;
            std::string Name;
            double Mean;
            double Var;
            int i = 0;
            while(std::getline(LineStream, cell, ','))
            {
                if(i == 0)  { Pop = ( std::stod(cell) == 0); }
                else if(i == 1) { Name = cell; }
                else if(i == 2) { Mean = std::stod(cell); }
                else if(i == 3) { Var = std::stod(cell); }
                ++i;
            }
            Name = Name.substr(1);
            if(Name.substr(0, Name.find_first_of("#")) != "Beta")
            {
                std::string ok = Name.substr(0, Name.find_first_of("#"));
                double a = 1+1;
            }
            
            if(Name == "Tau") { Var = 5.0; }
            //if(Name == "Ksi") { Var = 0.001; }
            //if(Name == "Delta") {Var = 0.0001; }
            //if(Name == "Beta") { Var = 0.0001; }
            
            auto RV = std::make_shared<GaussianRandomVariable>(Mean, Var);
            if(Pop) { m_PopulationRandomVariables.insert(RandomVariable(Name, RV)); }
            else { m_IndividualRandomVariables.insert(RandomVariable(Name, RV)); }
        }
        m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.0000001 );
    }
    else { std::cout << "Unable to open the parameters"; }
    */
    
    
    
    /// Tests
    TestAssert::WarningInequality_GreaterThan(m_Noise->GetVariance(), 0.0, "NetworkPropagationModel : wrong noise variance");
    
}

void 
NetworkPropagationModel2
::UpdateParameters(const Realizations& R,
                   const std::vector<std::string> Names) 
{
    /// This first part inspects the parameters names to update
    int UpdateCase = 1;
    
    for(auto it = Names.begin(); it != Names.end(); ++it) 
    {
        std::string Name = it->substr(0, it->find_first_of("#"));
    
        if(Name == "None" or Name == "Ksi" or Name == "Tau")
        {
            continue;
        }
        else if(Name == "S")
        {
            UpdateCase = std::max(UpdateCase, 2);
        }
        else if(Name == "Beta")
        {
            UpdateCase = std::max(UpdateCase, 3);
        }
        else if(Name == "P0")
        {
            UpdateCase = std::max(UpdateCase, 4);
        }
        else if(Name == "Delta" or Name == "All")
        {
            UpdateCase = 5;
            break;
        } 
        else
        {
            UpdateCase = 5;
            std::cout << "Should be" << Name << "be in NetworkPropagationModel > Update Parameters?" << std::endl;
            break;
        }
    }
    
    auto R1 = std::make_shared<Realizations>(R);
    
    /// Update Case in fonction of the the vector of names
    switch(UpdateCase)
    {
        case 1:
            break;
        case 2:
            ComputeSpaceShifts(*R1);
            break;
        case 3:
            ComputeAMatrix(*R1);
            ComputeSpaceShifts(*R1);
        case 4:
            ComputeOrthonormalBasis(*R1);
            ComputeAMatrix(*R1);
            ComputeSpaceShifts(*R1);
            break;
        case 5:
            ComputeInterpolationCoefficients(*R1);
            ComputeOrthonormalBasis(*R1);
            ComputeAMatrix(*R1);
            ComputeSpaceShifts(*R1);
            break;
        default:
            std::cout << "Error? NetworkPropagationModel > UpdateParameters";
            break;
    }
    
}

NetworkPropagationModel2::Data
NetworkPropagationModel2
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    auto R = SimulateRealizations(NumberOfSubjects);
    ComputeInterpolationCoefficients(R);
    ComputeOrthonormalBasis(R);
    ComputeAMatrix(R);
    ComputeSpaceShifts(R);
    
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::uniform_real_distribution<double> ObsDistrib(50, 95);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    std::shared_ptr<ExponentialCurveManifold> CastedManifold = std::static_pointer_cast<ExponentialCurveManifold>(m_Manifold);
    double T0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetMean();
    double KsiMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    VectorType V0(1, exp(KsiMean));
    VectorType P0(1, exp(R.at("P0")(0)));
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    
    Data D;
    double RealNoise = 0.0, Q = 0.0;
    
    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        IndividualData ID;
        std::vector<double> TimePoints;
        for(int j = 0; j < Uni(RNG); ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));
        }
        std::sort(TimePoints.begin(), TimePoints.end());
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        auto SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it, ++Q)
        {
            VectorType Scores(SpaceShift.size());
            double TimePoint = SubjectTimePoint(*it);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0,SpaceShift, TimePoint, PropagationCoefficients);
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
    
    std::cout << "Real Noise = " << RealNoise / Q << std::endl;
    std::cout << "Real Likelihood = " << ComputeLogLikelihood(R, D) << std::endl;
    
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
    
    
    return D;
}


double
NetworkPropagationModel2
::ComputeLogLikelihood(const Realizations& R, const Data& D) 
{
    /// Get the data
    std::shared_ptr<ExponentialCurveManifold> CastedManifold = std::static_pointer_cast<ExponentialCurveManifold>(m_Manifold);
    double T0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetMean();
    double KsiMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    VectorType V0(1, exp(KsiMean));
    VectorType P0(1, exp(R.at("P0")(0)));
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    
    /// Compute the likelihood
    double LogLikelihood = 0.0;
    unsigned int i = 0;
    for(auto IterData = D.begin(); IterData != D.end(); ++IterData, ++i)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        for(auto it = IterData->begin(); it != IterData->end(); ++it)
        {
            double TimePoint = SubjectTimePoint(it->second);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, PropagationCoefficients);
            LogLikelihood += (it->first - ParallelCurve).squared_magnitude();
        }
    }
    
    TestAssert::WarningInequality_GreaterThan(LogLikelihood, 0.0, "NetworkPropagationModel>ComputeIndividualLogLikelihood ; wrong Likelihood");
    
    LogLikelihood  /= -2*m_Noise->GetVariance();
    LogLikelihood -= m_NbTotalOfObservations*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
    
    return LogLikelihood ;
    
}

double
NetworkPropagationModel2
::ComputeIndividualLogLikelihood(const Realizations& R,
                                 const Data& D, const int SubjectNumber) 
{
    /// Get the data
    std::shared_ptr<ExponentialCurveManifold> CastedManifold = std::static_pointer_cast<ExponentialCurveManifold>(m_Manifold);
    double T0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetMean();
    double KsiMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    VectorType V0(1, exp(KsiMean));
    VectorType P0(1, exp(R.at("P0")(0)));
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(SubjectNumber));
    
    /// Compute the likelihood
    double LogLikelihood = 0.0;
    auto N = D.at(SubjectNumber).size(); 
    
#pragma omp parallel for reduction(+:LogLikelihood)    
    for(size_t i = 0; i < N; ++i)
    {
        auto& it = D.at(SubjectNumber).at(i);
        double TimePoint = SubjectTimePoint(it.second);
        VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, PropagationCoefficients);
        LogLikelihood += (it.first - ParallelCurve).squared_magnitude();
    }
    
    TestAssert::WarningInequality_GreaterThan(LogLikelihood, 0.0, "NetworkPropagationModel>ComputeIndividualLogLikelihood ; wrong Likelihood");

    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
    
    return LogLikelihood;
}

NetworkPropagationModel2::SufficientStatisticsVector
NetworkPropagationModel2
::GetSufficientStatistics(const Realizations& R,
                          const Data& D) 
{
    //// Initialization
    std::shared_ptr<ExponentialCurveManifold> CastedManifold = std::static_pointer_cast<ExponentialCurveManifold>(m_Manifold);
    double T0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetMean();
    double KsiMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    VectorType V0(1, exp(KsiMean));
    VectorType P0(1, exp(R.at("P0")(0)));
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    double NumberOfSubjects = R.at("Ksi").size();
    
    
    /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
    
    VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
    int i = 0;
    auto IterS1 = S1.begin(), IterS2 = S2.begin();
    for(auto IterD = D.begin(); IterD != D.end() && IterS1 != S1.end() && IterS2 != S2.end(); ++IterD, ++i)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        for(auto it = IterD->begin(); it != IterD->end(); ++it)
        {
            double TimePoint = SubjectTimePoint(it->second);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, PropagationCoefficients);
            *IterS1 = dot_product(ParallelCurve, it->first);
            *IterS2 = ParallelCurve.squared_magnitude();
            ++IterS1, ++IterS2;
        }
    }
     
 /*   
    VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
    ScalarType * it1 = S1.memptr();
    ScalarType * it2 = S2.memptr();
    

#pragma omp simd
    for(size_t i = 0; i < m_NbTotalOfObservations, ++i)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        auto NbObs = D->at(i).size();
        
#pragma omp simd        
        for(size_t j = 0; j < NbObs; ++j)
        {
            auto& it = D->at(i).at(j);
            double TimePoint = SubjectTimePoint(it.second);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, PropagationCoefficients);
            it1[i*+j]
        }
    }
    */
    
    
    /// S3 <- Ksi_i   &    S4 <- Ksi_i * Ksi_i
    VectorType S3(NumberOfSubjects), S4(NumberOfSubjects);
    auto IterKsi = R.at("Ksi").begin();
    auto IterS3 = S3.begin(), IterS4 = S4.begin();
    for(    ; IterKsi != R.at("Ksi").end() && IterS3 != S3.end() && IterS4 != S4.end(); ++IterKsi, ++IterS3, ++IterS4)
    {
        *IterS3 = *IterKsi;
        *IterS4 = *IterKsi * *IterKsi;
    }
    
    
    /// S5 <- Tau_i   &    S6 <- Tau_i * Tau_i
    VectorType S5(NumberOfSubjects), S6(NumberOfSubjects);
    auto IterTau = R.at("Tau").begin();
    auto IterS5 = S5.begin(), IterS6 = S6.begin();
    for(    ; IterTau != R.at("Tau").end() && IterS5 != S5.end() && IterS6 != S6.end(); ++IterTau, ++IterS5, ++IterS6)
    {
        *IterS5 = *IterTau;
        *IterS6 = *IterTau * *IterTau;
    }
    
    /// S6 <- P0(k)
    VectorType S7(1, R.at("P0")(0));
    
    /// S8 <- beta_k
    VectorType S8((m_Manifold->GetDimension() - 1)*m_NbIndependentComponents);
    i = 0;
    for(auto it = S8.begin(); it != S8.end(); ++it, ++i)
    {
        *it = R.at("Beta#" + std::to_string(i))(0);
    }
    
    /// S9 <- delta_k
    VectorType S9(m_NbControlPoints - 1);
    i = 1;
    for(auto it = S9.begin();  it != S9.end(); ++it, ++i)
    {
        *it = R.at("Delta#" + std::to_string(i))(0);
    }
    
    SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9};
    
    /// TESTS NEEDED
    
    /// END TESTS
    
    return S;
    

}

void
NetworkPropagationModel2
::UpdateRandomVariables(const SufficientStatisticsVector &SS, const Data& D) 
{
    double NumberOfSubjects = D.size();
    
    /// Update P0
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    P0->SetMean(SS[6](0));
    
    /// Update Beta_k
    int i = 0;
    for(auto it = SS[7].begin(); it != SS[7].end(); ++it, ++i)
    {
        auto AbstractBeta = m_PopulationRandomVariables.at("Beta#" + std::to_string(i));
        auto Beta = std::static_pointer_cast<GaussianRandomVariable>(AbstractBeta);
        Beta->SetMean(*it);
    }
    
    /// Update Delta_k
    i = 1;
    for(auto it = SS[8].begin(); it != SS[8].end(); ++it, ++i)
    {
        auto AbstractDelta = m_PopulationRandomVariables.at("Delta#" + std::to_string(i));
        auto Delta = std::static_pointer_cast<GaussianRandomVariable>(AbstractDelta);
        Delta->SetMean(*it);
    }
    
    /// Update Ksi_mean and Ksi_var
    double KsiMean = 0.0, KsiVar = 0.0;
    for(auto it = SS[2].begin(); it != SS[2].end(); ++it)
    {
        KsiMean += *it;
    }
    KsiMean /= NumberOfSubjects;
    for(auto it = SS[3].begin(); it != SS[3].end(); ++it)
    {
        KsiVar += *it;
    }
    KsiVar -= NumberOfSubjects * KsiMean * KsiMean;
    KsiVar /= NumberOfSubjects;
    
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    Ksi->SetMean(KsiMean);
    Ksi->SetVariance(KsiVar);
    
    /// Update Tau_Mean = T0 and Tau_Var
    double TauMean = 0.0, TauVar = 0.0;
    for(auto it = SS[4].begin(); it != SS[4].end(); ++it)
    {
        TauMean += *it;
    }
    TauMean /= NumberOfSubjects;
    for(auto it = SS[5].begin(); it != SS[5].end(); ++it)
    {
        TauVar += *it;
    }
    TauVar -= NumberOfSubjects * TauMean * TauMean;
    TauVar /= NumberOfSubjects;
    
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    Tau->SetMean(TauMean);
    Tau->SetVariance(TauVar);
    
    /// Update Noise
    double NoiseVariance = m_SumObservations;
    for(auto IterS1 = SS[0].begin(), IterS2 = SS[1].begin(); IterS1 != SS[0].end() && IterS2 != SS[1].end(); ++IterS1, ++IterS2)
    {
        NoiseVariance += - 2 * *IterS1 + *IterS2;
    }
    NoiseVariance /= m_NbTotalOfObservations * m_Manifold->GetDimension();
    m_Noise->SetVariance(NoiseVariance);
        
    TestAssert::WarningInequality_GreaterThan(NoiseVariance, 0.0, "LongitudinalModel>UpdateRandomVariables ; wrong noise variance");
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void
NetworkPropagationModel2
::ComputeOutputs() 
 {
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));   
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));   
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));   
    auto Beta0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Beta#0"));   
    auto Beta1 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Beta#1"));
    auto Delta1 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#1"));
    auto Delta2 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#2"));
    auto Delta170 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#170"));
    double Sigma = m_Noise->GetVariance();
    
    std::cout << "P0 : " << exp(P0->GetMean()) << " - T0 : " << Tau->GetMean() << " - Var(Tau) : " << Tau->GetVariance();
    std::cout << " - V0 : " << exp(Ksi->GetMean()) << " - Var(Ksi) : " << Ksi->GetVariance() << " - Beta0 : ";
    std::cout << Beta0->GetMean() << " - Beta1 : " << Beta1->GetMean() << " - Delta1 : " << Delta1->GetMean();
    std::cout << " Delta2 : " << Delta2->GetMean() << " - Delta#170 : " << Delta170->GetMean() << " - Sigma : " << Sigma << std::endl;
}   


void 
NetworkPropagationModel2
::SaveData(unsigned int IterationNumber) 
{
    /// Save the data in case it crashes and a re-init is needed
    std::ofstream ModelParameters;    
    ModelParameters.open("ModelParametersMciConverters_Linear.txt", std::ofstream::out | std::ofstream::trunc);
    ModelParameters << "Iteration : " << IterationNumber << std::endl;
    for(auto it = m_PopulationRandomVariables.begin(); it != m_PopulationRandomVariables.end(); ++it)
    {
        auto RV = std::static_pointer_cast<GaussianRandomVariable>(it->second);
        ModelParameters << "0, " << it->first << ", " << RV->GetMean() << ", " << RV->GetVariance() << std::endl;
    }
    for(auto it = m_IndividualRandomVariables.begin(); it != m_IndividualRandomVariables.end(); ++it)
    {
        auto RV = std::static_pointer_cast<GaussianRandomVariable>(it->second);
        ModelParameters << "1, " << it->first << ", " << RV->GetMean() << ", " << RV->GetVariance() << std::endl;
    }
    
    
    
    /// Save the delta_k for visualization
    std::ofstream DeltaViz;    
    DeltaViz.open("DeltaVisualizationMciConverters_Linear.txt", std::ofstream::out | std::ofstream::trunc);    
        
    for(auto it = m_PopulationRandomVariables.begin(); it != m_PopulationRandomVariables.end(); ++it)
    {
        std::string Name = it->first;
        if(Name.substr(0, Name.find_first_of("#")) == "Delta")
        {
            std::string Number = Name.substr(Name.find_first_of("#") + 1);
            auto RV = std::static_pointer_cast<GaussianRandomVariable>(it->second);
            DeltaViz << Number << ", " << RV->GetMean() << std::endl;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void
NetworkPropagationModel2
::InitializeFakeRandomVariables() 
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    /////////////////////////////
    /// Population Parameters ///
    /////////////////////////////
    
    auto P0 = std::make_shared<GaussianRandomVariable>(20, 0.0000001);
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0) );
    
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.01 );
    
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension() - 1); ++i)
    {
        auto Beta = std::make_shared<GaussianRandomVariable>((double)i/20.0, 0.0001);
        std::string Name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Beta));
    }
    
    for(int i = 0; i < m_NbControlPoints - 1; ++i)
    {
        auto Delta = std::make_shared<GaussianRandomVariable>(-0.5 - (double)i/20, 0.0001);
        std::string Name = "Delta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(Name, Delta));
    }    
    
    /////////////////////////////
    /// Individual Parameters ///
    /////////////////////////////
    
    auto Ksi = std::make_shared<GaussianRandomVariable>(0.0, 0.00001);
    auto Tau = std::make_shared<GaussianRandomVariable>(70.0, 3.0);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        auto S = std::make_shared<GaussianRandomVariable>(0.0, 0.05);
        std::string Name = "S#" + std::to_string(i);
        m_IndividualRandomVariables.insert(RandomVariable(Name, S));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

NetworkPropagationModel2::VectorType
NetworkPropagationModel2
::GetPropagationCoefficients(const Realizations& R) 
{
    VectorType PropagationCoefficients = m_InterpolationMatrix * m_InterpolationCoefficients;
    return PropagationCoefficients;    
}

std::function<double(double)>
NetworkPropagationModel2
::GetSubjectTimePoint(const int SubjectNumber, const Realizations& R) 
{
    double AccFactor = exp(R.at("Ksi")(SubjectNumber));
    double TimeShift = R.at("Tau")(SubjectNumber);
    
    return [AccFactor, TimeShift](double t) { return AccFactor * (t - TimeShift); };
}

void
NetworkPropagationModel2
::ComputeInterpolationCoefficients(const Realizations& R) 
{
    VectorType Delta(m_NbControlPoints, 0.0);
    
    int i = 1;
    for(auto it = Delta.begin() + 1; it != Delta.end(); ++it, ++i)
    {
        *it = R.at("Delta#" + std::to_string(i))(0);
    }
    
    m_InterpolationCoefficients = m_InvertKernelMatrix * Delta;
}

void
NetworkPropagationModel2
::ComputeOrthonormalBasis(const Realizations& R) 
{
    /// Initialization
    std::shared_ptr<ExponentialCurveManifold> CastedManifold = std::static_pointer_cast<ExponentialCurveManifold>(m_Manifold);
    double T0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetMean();
    double KsiMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    VectorType V0(1, exp(KsiMean));
    VectorType P0(1, exp(R.at("P0")(0)));
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    
    /// Compute the transformation to do the Householder reflection in a Euclidean space
    VectorType U = CastedManifold->GetVelocityTransformToEuclideanSpace(P0, T0, V0, PropagationCoefficients);

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
    // TODO : mettre un std::move(m_OrthogonalBasis, Q) ...
    m_OrthogonalBasis = Q;
}

void 
NetworkPropagationModel2
::ComputeAMatrix(const Realizations& R) 
{
    MatrixType AMatrix(m_Manifold->GetDimension(), m_NbIndependentComponents);
    
    for(int i = 0; i < m_NbIndependentComponents ; ++i)
    {
        VectorType Beta(m_Manifold->GetDimension() - 1);
        for(int j = 0; j < m_Manifold->GetDimension() - 1 ; ++j)
        {
            std::string Number = std::to_string(int(j + i*(m_Manifold->GetDimension() - 1)));
            Beta(j) = R.at( "Beta#" + Number)(0);
        }
        
        VectorType V = LinearCombination(Beta, m_OrthogonalBasis)   ;     
        AMatrix.set_column(i, V);
    }
    
    m_AMatrix = AMatrix;
    
    /// TESTS NEEDED
    // TODO : Mettre des tests
    /// END TESTS
}

void 
NetworkPropagationModel2
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

    m_SpaceShifts = SpaceShifts;
    
    /// TESTS NEEDED
    // TODO : Mettre des tests
    /// END TESTS
}