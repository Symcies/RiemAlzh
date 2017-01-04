#include "UnivariateModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

UnivariateModel
::UnivariateModel(std::shared_ptr<AbstractBaseManifold>& M) 
{
    m_OutputParameters.open("Parameters.txt", std::ofstream::out | std::ofstream::trunc);
    m_BaseManifold = M;
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
::Initialize(const Data& D)
{
    
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.8, 0.01*0.01 );
    auto Tau = std::make_shared< GaussianRandomVariable >(60.0, 7.0*7.0);
    auto Ksi = std::make_shared< GaussianRandomVariable >(log(0.05), 1.0*1.0); 
    
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.1 * 0.1);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
}


void 
UnivariateModel
::UpdateModel(const Realizations &R, const std::vector<std::string> Names) 
{
    // TODO : Check if something has to be added
}


UnivariateModel::SufficientStatisticsVector
UnivariateModel
::GetSufficientStatistics(const Realizations& R, const Data& D) 
{
    /// Get the data to compute the geometric operation on the manifold
    double T0 = 0.0;
    double V0 = 1.0;
    double P0 = R.at("P0")(0);
    unsigned int NumberOfSubjects = (int)D.size();
    
    /// Compute S0 <-- Sum(y_ijk)
    double SumData = 0.0;
    unsigned int K = 0;
    for(auto it = D.begin(); it != D.end(); ++it)
    {
        K += it->size();
        for(auto it2 = it->begin(); it2 != it->end(); ++it2)
        {
            SumData += it2->first.squared_magnitude();
        }
    }
    VectorType S0(1, SumData);
    
    
    /// Compute S1 <- [Data_ij * ParallelCurve_ij]  and  S2 <- [ParrallelCurve_ij ^2]
    VectorType S1(K, 0), S2(K, 0);
    
    int i = 0;
    auto IterS1 = S1.begin(), IterS2 = S2.begin();
    for(auto IterData = D.begin(); IterData != D.end() && IterS1 != S1.end() && IterS2 != S2.end(); ++i, ++IterData)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        for(auto IterIndivData = IterData->begin(); IterIndivData != IterData->end(); ++IterIndivData)
        {
            double TimePoint = SubjectTimePoint(IterIndivData->second);
            
            auto ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
            auto Observation = IterIndivData->first[0];
            
            *IterS1 = ParallelCurve * Observation;
            *IterS2 = ParallelCurve * ParallelCurve;
            ++IterS1, ++IterS2;
        }
    }
    
    /// Compute S3 <- P0(k)
    VectorType S3(1, R.at("P0")(0));
    
    /// Compute S4 <- tau_i  and  S5 <- tau_i * tau_i
    VectorType S4(NumberOfSubjects), S5(NumberOfSubjects);
    auto IterTau = R.at("Tau").begin();
    auto IterS4 = S4.begin();
    auto IterS5 = S5.begin();
    for(    ; IterS4 != S4.end() && IterS5 != S5.end() && IterTau != R.at("Tau").end(); ++IterS4, ++IterS5, ++IterTau)
    {
        *IterS4 = *IterTau;
        *IterS5 = *IterTau * *IterTau;
    }
    
    /// Compute S5 <- ksi_i and 6 <- [ksi_i * ksi_i]
    VectorType S6(NumberOfSubjects, 0), S7(NumberOfSubjects, 0);
    auto IterKsi = R.at("Ksi").begin();
    auto IterS6 = S6.begin();
    auto IterS7 = S7.begin();
    for(    ; IterS6 != S6.end() && IterS7 != S7.end() && IterKsi != R.at("Ksi").end(); ++IterS6, ++IterS7, ++IterKsi)
    {
        *IterS6 = *IterKsi;
        *IterS7 = *IterKsi * *IterKsi;
    }
    
    SufficientStatisticsVector S = {S0, S1, S2, S3, S4, S5, S6, S7 };
    
    /// Tests
    /// There are two way to compute the logLikelihood. Check if they are equal
    std::function<double()> f1 = [=, &D, &R] () { return this->ComputeLogLikelihood(R, D); };
    std::function<double()> f2 = [=, &D, &S] ()
    {
        double Calculation = S[0](0);
        for(auto it : S[1])
        {
            Calculation += -2 * it;
        }
        for(auto it : S[2])
        {
            Calculation += it;
        }
        Calculation /= -2*this->m_Noise->GetVariance();
        Calculation -= K * log(sqrt( 2 * M_PI * m_Noise->GetVariance() ));        
        return Calculation;
    };
    
    TestAssert::WarningEquality_Function(f1, f2, "Likelihood not correct. UnivariateModel > GetSufficientStatistics");
    
    return S;
}

void
UnivariateModel
::UpdateRandomVariables(const SufficientStatisticsVector &StochSufficientStatistics,
                        const Data& D) 
{
    double NumberOfSubjects = D.size();
    
    /// Update P0
    auto AbstractP0 = m_PopulationRandomVariables.at("P0");
    auto P0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractP0 );
    P0->SetMean( StochSufficientStatistics[3](0));
    
    
   /// Update mean and variance of ksi
    double KsiMean = 0, KsiVar = 0;
    for(auto IterV0 = StochSufficientStatistics[6].begin(); IterV0 != StochSufficientStatistics[6].end(); ++IterV0)
    {
        KsiMean += *IterV0;
    }
    KsiMean /= NumberOfSubjects;
    
    for(auto IterKsi = StochSufficientStatistics[7].begin(); IterKsi != StochSufficientStatistics[7].end(); ++IterKsi)
    {
        KsiVar += *IterKsi;
    }
    
    KsiVar -= NumberOfSubjects * KsiMean * KsiMean;
    KsiVar /= NumberOfSubjects;
    
    auto AbstractKsi = m_IndividualRandomVariables.at("Ksi");
    auto Ksi = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractKsi );
    
    Ksi->SetMean(KsiMean); 
    Ksi->SetVariance(KsiVar);
    
    /// Update Mean and Variance of Tau
    double TauMean = 0, TauVar = 0;
    for(auto IterT0 = StochSufficientStatistics[4].begin(); IterT0 != StochSufficientStatistics[4].end(); ++IterT0)
    {
        TauMean += *IterT0;
    }
    
    TauMean /= NumberOfSubjects;
    
    for(auto IterTau = StochSufficientStatistics[5].begin(); IterTau != StochSufficientStatistics[5].end(); ++IterTau)
    {
        TauVar += *IterTau;
    }
    
    TauVar -= NumberOfSubjects * TauMean * TauMean;
    TauVar /= NumberOfSubjects;
    
    auto AbstractTau = m_IndividualRandomVariables.at("Tau");
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>( AbstractTau );
    
    Tau->SetMean(TauMean);
    Tau->SetVariance(TauVar);
    
    
    /// Update Uncertainty variance
    double K = 0.0;
    for(auto it = D.begin(); it != D.end(); ++it)
    {
        K += it->size();
    }
    
    /// Sum Yijk²
    double NoiseVariance = StochSufficientStatistics[0](0);


    /// Sum -2 S1 + S2
    auto IterS1 = StochSufficientStatistics[1].begin();
    auto IterS2 = StochSufficientStatistics[2].begin();
    for( ; IterS1 != StochSufficientStatistics[1].end() && IterS2 != StochSufficientStatistics[2].end(); ++IterS1, ++IterS2)
    {
        NoiseVariance += - 2 * *IterS1 + *IterS2;
    }

    /// Divide by K
    NoiseVariance /= K;
    m_Noise->SetVariance(NoiseVariance);
    
}


double 
UnivariateModel
::ComputeLogLikelihood(const Realizations& R, const Data& D) 
{
    /// Get the data
    double T0 = 0.0;
    double P0 = R.at("P0")(0);
    double V0 = 1.0; 
    
    /// Compute the loglikelihood
    double LogLikelihood = 0, K = 0;
    int i = 0;
    for(auto IterData = D.begin(); IterData != D.end(); ++i, ++IterData)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        K += IterData->size();
        for(const auto& IterIndivData : *IterData)
        {
            double TimePoint = SubjectTimePoint(IterIndivData.second);
            auto ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
            double Norm = (IterIndivData.first[0] - ParallelCurve);
            LogLikelihood += Norm * Norm ;
        }
    }
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    
    LogLikelihood -= K * log(sqrt( 2 * M_PI * m_Noise->GetVariance() ));
    
    return LogLikelihood;
}

double 
UnivariateModel
::ComputeIndividualLogLikelihood(const Realizations& R,
                                 const Data& D, const int SubjectNumber) 
{
    /// Initialize the individual parameters
    double T0 = 0.0;
    double P0 = R.at("P0")(0);
    double V0 = 1.0;
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    
    /// Compute the loglikelihood
    double LogLikelihood = 0;
    double k = D.at(SubjectNumber).size();
    int i = 0;
    for(auto IterData = D.at(SubjectNumber).begin(); IterData != D.at(SubjectNumber).end(); ++IterData, ++i)
    {
        double TimePoint = SubjectTimePoint(IterData->second);
        auto ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
        double Norm = (IterData->first[0] - ParallelCurve);
        LogLikelihood += Norm*Norm;
    }

    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= k * log(sqrt( 2 * m_Noise->GetVariance() * M_PI));
    return LogLikelihood;
}


UnivariateModel::Data
UnivariateModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    /// Simulate Realizations
    auto R = SimulateRealizations(NumberOfSubjects);
    
    /// Initialize
    std::ofstream IOTimePoints, IOSimulatedData;
    IOTimePoints.open("TimePoints.txt", std::ofstream::out | std::ofstream::trunc);
    IOSimulatedData.open("SimulatedData.txt", std::ofstream::out | std::ofstream::trunc);
    
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::uniform_real_distribution<double> ObsDistrib(50, 95);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    Data D;
    double T0 = 0.0;
    double P0 = R.at("P0")(0);
    double V0 = 1.0;
    
    
    /// Simulate the data
    double q = 0, SumNoise;
    
    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        int NbObservations = Uni(RNG);
        IndividualData InDa;
        
        /// Generate time points
        std::vector<double> TimePoints;
        for(int j = 0; j < NbObservations; ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));   
        }
        std::sort(TimePoints.begin(), TimePoints.end());
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        
        /// Generate observations corresponding to the time points
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it)
        {
            q += 1;
            double TimePoint = SubjectTimePoint(*it);
            double Noise = NoiseDistrib(RNG);
            SumNoise += Noise;
            double ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
            VectorType Scores(1, ParallelCurve + Noise);
            
            InDa.push_back( std::pair< VectorType, double> (Scores, *it));
            
            IOTimePoints << *it;
            IOSimulatedData << ParallelCurve + Noise;
            if(it != TimePoints.end() -1)
            {
                IOTimePoints << " ";
                IOSimulatedData << " ";
            }
        }
        IOTimePoints << std::endl;
        IOSimulatedData << std::endl;
        
        D.push_back(InDa);
    }
    double A = ComputeLogLikelihood(R, D) ; 
    std::cout << "Likelihood : " << A << std::endl;
    
    /*
    for(double P = 0.36; P < 0.43; P += 0.005)
    {
        R->at("P0")(0) = P;
        std::cout << "LogLikelihood with P0 = " << P << " : " << ComputeLogLikelihood(R, std::make_shared<Data>(D)) << std::endl;
    }
     */
    
    std::cout << "Real Noise : " << SumNoise/q << std::endl;
    std::cout << "P0 / T0 / V0 : " << P0 << " / " << T0 << " / " << V0 << std::endl;
    
    return D;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
UnivariateModel
::DisplayOutputs() 
{
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    double Sigma = m_Noise->GetVariance();
        
    std::cout << "Parameters : P0: " << P0->GetMean() << ". T0: " << Tau->GetMean() << ". Var(Tau): " << Tau->GetVariance();
    std::cout << ". V0 : " << exp(Ksi->GetMean()) << ". Var(Ksi): " << Ksi->GetVariance() << ". Sigma: " << Sigma << ". " << std::endl;
    
    
}



void 
UnivariateModel
::SaveData(unsigned int IterationNumber, const Realizations& R) 
{
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
UnivariateModel
::InitializeFakeRandomVariables() 
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.24, 0.01*0.01 );
    auto Tau = std::make_shared< GaussianRandomVariable >(70.0, 7.0*7.0);
    auto Ksi = std::make_shared< GaussianRandomVariable >(log(0.034), 0.5*0.5); 
    
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.01 * 0.01);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::function<double(double)>  
UnivariateModel
::GetSubjectTimePoint(const int SubjectNumber, const Realizations& R)
{
    double AccFactor = exp(R.at("Ksi")(SubjectNumber));
    double TimeShift = R.at("Tau")(SubjectNumber);
    
    return [AccFactor, TimeShift](double t) { return AccFactor * (t - TimeShift); };
}



