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
::Initialize()
{
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.5, 0.00001 );
    auto Tau = std::make_shared< GaussianRandomVariable >(72.0, 2*2);
    auto Ksi = std::make_shared< GaussianRandomVariable >(-2.5, 0.020); 
    
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.00001);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
}


void 
UnivariateModel
::UpdateParameters(const std::shared_ptr<MultiRealizations> &R, const std::vector<std::string> Names) 
{
    // TODO : Check if something has to be added
}


std::map< std::string, double >
UnivariateModel
::GetParameters() 
{
    auto P0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    
    std::map<std::string, double> Parameters;
    
    Parameters["P0"] = P0->GetMean();
    Parameters["V0"] = Ksi->GetMean();
    Parameters["Ksi"] = Ksi->GetVariance();
    Parameters["T0"] = Tau->GetMean();
    Parameters["Tau"] = Tau->GetVariance();
    Parameters["NoiseVariance"] = m_Noise->GetVariance();
    
    return Parameters;
}


UnivariateModel::SufficientStatisticsVector
UnivariateModel
::GetSufficientStatistics(const std::shared_ptr<MultiRealizations> &R, const std::shared_ptr<Data> &D) 
{
    /// Get the data to compute the geometric operation on the manifold
    double T0 = GetInitialTime();
    double V0 = GetInitialVelocity();
    double P0 = GetInitialPosition(R);
    unsigned int NumberOfSubjects = (int)D->size();
    
    /// Compute S0 <-- Sum(y_ijk)
    double SumData = 0.0;
    unsigned int K = 0;
    for(auto it = D->begin(); it != D->end(); ++it)
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
    for(auto IterData = D->begin(); IterData != D->end() && IterS1 != S1.end() && IterS2 != S2.end(); ++i, ++IterData)
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
    
    /// Compute S3 <- [ksi_i * ksi_i]  and  S4 <- ksi_i
    VectorType S3(NumberOfSubjects, 0), S4(NumberOfSubjects, 0);
    auto IterKsi = R->at("Ksi").begin();
    auto IterS3 = S3.begin();
    auto IterS4 = S4.begin();
    for(    ; IterS3 != S3.end() && IterS4 != S4.end() && IterKsi != R->at("Ksi").end(); ++IterS3, ++IterS4, ++IterKsi)
    {
        *IterS3 = *IterKsi * *IterKsi;
        *IterS4 = *IterKsi;
    }
    
    /// Compute S5 <- [tau_i * tau_i]  and  S6 <- tau_i
    VectorType S5(NumberOfSubjects), S6(NumberOfSubjects);
    auto IterTau = R->at("Tau").begin();
    auto IterS5 = S5.begin();
    auto IterS6 = S6.begin();
    for(    ; IterS5 != S5.end() && IterS6 != S6.end() && IterTau != R->at("Tau").end(); ++IterS5, ++IterS6, ++IterTau)
    {
        *IterS5 = *IterTau * *IterTau;
        *IterS6 = *IterTau;
    }
    
    SufficientStatisticsVector S = {S0, S1, S2, S3, S4, S5, S6, VectorType(1, P0) };
    
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
                        const std::shared_ptr<Data> &D) 
{
    double NumberOfSubjects = D->size();
    
    /// Update P0
    auto AbstractP0 = m_PopulationRandomVariables.at("P0");
    auto P0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractP0 );
    P0->SetMean( StochSufficientStatistics[7](0));
    
   /// Update mean and variance of ksi
    double KsiMean = 0, KsiVar = 0;
    for(auto IterV0 = StochSufficientStatistics[4].begin(); IterV0 != StochSufficientStatistics[4].end(); ++IterV0)
    {
        KsiMean += *IterV0;
    }
    KsiMean /= NumberOfSubjects;
    
    for(auto IterKsi = StochSufficientStatistics[3].begin(); IterKsi != StochSufficientStatistics[3].end(); ++IterKsi)
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
    for(auto IterT0 = StochSufficientStatistics[6].begin(); IterT0 != StochSufficientStatistics[6].end(); ++IterT0)
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
    for(const auto& it : *D)
    {
        K += it.size();
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
::ComputeLogLikelihood(const std::shared_ptr<MultiRealizations> &R, const std::shared_ptr<Data> &D) 
{
    /// Get the data
    double T0 = GetInitialTime();
    double P0 = GetInitialPosition(R);
    double V0 = GetInitialVelocity();
    
    //std::cout << "T0/P0/V0 : " << T0 << " / " << P0 << " / " << V0 << std::endl;
    
    /// Compute the loglikelihood
    double LogLikelihood = 0, K = 0;
    int i = 0;
    for(auto IterData = D->begin(); IterData != D->end(); ++i, ++IterData)
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
    
    //std::cout << "Noise (param) : " << m_Noise->GetVariance() << std::endl;
    LogLikelihood /= -2*m_Noise->GetVariance();
    
    LogLikelihood -= K * log(sqrt( 2 * M_PI * m_Noise->GetVariance() ));
    
    return LogLikelihood;
}

double 
UnivariateModel
::ComputeIndividualLogLikelihood(const std::shared_ptr<MultiRealizations> &R,
                                 const std::shared_ptr<Data> &D, const int SubjectNumber) 
{
    /// Initialize the individual parameters
    double T0 = GetInitialTime();
    double P0 = GetInitialPosition(R);
    double V0 = GetInitialVelocity();
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    
    /// Compute the loglikelihood
    double LogLikelihood = 0;
    double k = D->at(SubjectNumber).size();
    int i = 0;
    for(auto IterData = D->at(SubjectNumber).begin(); IterData != D->at(SubjectNumber).end(); ++IterData, ++i)
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
    auto R = std::make_shared<MultiRealizations>( SimulateRealizations(NumberOfSubjects) );
    
    /// Initialize
    std::ofstream IOTimePoints, IOSimulatedData;
    IOTimePoints.open("TimePoints.txt", std::ofstream::out | std::ofstream::trunc);
    IOSimulatedData.open("SimulatedData.txt", std::ofstream::out | std::ofstream::trunc);
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::uniform_real_distribution<double> ObsDistrib(45, 95);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    Data D;
    double T0 = GetInitialTime();
    double P0 = GetInitialPosition(R);
    double V0 = GetInitialVelocity();
    
    
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
    double A = ComputeLogLikelihood(R, std::make_shared<Data>(D)) ; 
    std::cout << "Real Noise : " << SumNoise/q << std::endl;
    std::cout << "P0 / T0 / V0 : " << P0 << " / " << T0 << " / " << V0 << std::endl;
    std::cout << "Likelihood : " << A << std::endl;
    
    
    
    
    
    return D;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
UnivariateModel
::ComputeOutputs() 
{
    double MeanP0 = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"))->GetMean();
    double TauMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetMean();
    double TauVar = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetVariance();
    double KsiMean = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetMean();
    double KsiVar = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetVariance();
    double Sigma = m_Noise->GetVariance();
        
    std::cout << "Parameters : P0: " << MeanP0 << ". TauMean: " << TauMean << ". TauVar: " << TauVar << ". KsiMean: " << KsiMean << ". KsiVar: " << KsiVar << ". Sigma: " << Sigma << ". " << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
UnivariateModel
::InitializeFakeRandomVariables() 
{
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.4, 0.001 );
    auto Tau = std::make_shared< GaussianRandomVariable >(75.0, 10);
    auto Ksi = std::make_shared< GaussianRandomVariable >(-3, 0.3); 
    
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.000001);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::function<double(double)>  
UnivariateModel
::GetSubjectTimePoint(const int SubjectNumber, const std::shared_ptr<MultiRealizations>& R)
{
    double AccFactor = exp(R->at("Ksi")(SubjectNumber));
    double V0 = GetInitialVelocity();
    double TimeShift = R->at("Tau")(SubjectNumber);
    double T0 = GetInitialTime();
    
    return [AccFactor, TimeShift, T0, V0](double t) { return AccFactor / V0 * (t - TimeShift) + T0; };
}

double
UnivariateModel
::GetInitialTime() 
{
    auto Tau = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    return Tau->GetMean();
}

double 
UnivariateModel
::GetInitialPosition(const std::shared_ptr<MultiRealizations>& R)
{
    double P0 = R->at("P0")(0);
    return P0;
    
    // It's not a logit function anymore
    double LogitP0 = 1.0 / (1.0 + exp(-P0));
    
    return LogitP0;
}

double 
UnivariateModel
::GetInitialVelocity() 
{
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    return exp(Ksi->GetMean());
}





