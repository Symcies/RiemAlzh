#include "UnivariateModel.h"
#include "../RandomVariables/ConstantRandomVariable.h"


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
::InitializeRandomVariables() 
{
    /// Initialization
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    m_ManifoldRandomVariables.clear();
    
    /// Initial Position
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.35, 0.01 );
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    
    /// Initial Time
    auto T0 = std::make_shared<ConstantRandomVariable>( 70.0);
    m_PopulationRandomVariables.insert( RandomVariable("T0", T0));
    
    /// Initial Velocity
    auto V0 = std::make_shared<GaussianRandomVariable>( 0.06, 0.001 );
    m_PopulationRandomVariables.insert( RandomVariable("V0", V0) );
    
    /// Noise
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.05);
    
    /// Individual Time Shift
    auto Tau = std::make_shared< GaussianRandomVariable >(0.0, 3.0);
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    
    /// Individual pace / Preaceleration factor
    auto Ksi = std::make_shared< GaussianRandomVariable >(0.0, 1.5); 
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
}


void 
UnivariateModel
::UpdateParameters(const std::shared_ptr<Realizations> &R, std::string Name) 
{
    // TODO : Check if something has to be added
}


std::vector< std::vector< double >> 
UnivariateModel
::GetSufficientStatistics(const std::shared_ptr<Realizations> &R, const std::shared_ptr<Data> &D) 
{
   
    /// Get the data to compute the geometric operation on the manifold
    double T0 = R->at("T0")[0];
    double V0 = R->at("V0")[0];
    double P0 = R->at("P0")[0];
    
    /// Compute S1 <- [Data_ij * ParallelCurve_ij]  and  S2 <- [ParrallelCurve_ij ^2]
    std::vector<double> S1, S2;
    
    int i =0;
    for(auto IterData : *D)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        for(auto IterIndivData : IterData)
        {
            double TimePoint = SubjectTimePoint(IterIndivData.second);
            double SpaceShift = 0;
            
            auto ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint);
            auto Observation = IterIndivData.first;
            
            S1.push_back( ParallelCurve * Observation[0] );
            S2.push_back( ParallelCurve * ParallelCurve );
            
        }
    }
    
    /// Compute S3 <- [ksi_i * ksi_i]  and  S4 <- [tau_i * tau_i]
    std::vector<double> S3, S4;
    for(int i = 0; i < D->size(); ++i)
    {
        double Ksi = R->at("Ksi")[i];
        double Tau = R->at("Tau")[i];
        
        S3.push_back(Ksi * Ksi);
        S4.push_back(Tau * Tau);
    }
    
    /// Compute S5 <- P0  and  S6 <- T0  and S7 <- V0
    std::vector<double> S5 = {P0};
    std::vector<double> S6 = {T0};
    std::vector<double> S7 = {V0};
    
    SufficientStatisticsVector SufficientStatistics = {S1, S2, S3, S4, S5, S6, S7};
    
    return SufficientStatistics;
}

void
UnivariateModel
::UpdateRandomVariables(const SufficientStatisticsVector &StochSufficientStatistics,
                        const std::shared_ptr<Data> &D) 
{
    
    double NumberOfSubject = D->size();
    
    auto AbstractP0 = m_PopulationRandomVariables.at("P0");
    auto P0 = std::dynamic_pointer_cast< GaussianRandomVariable> (AbstractP0);
    P0->SetMean( StochSufficientStatistics[4][0]);
    
    auto AbstractT0 = m_PopulationRandomVariables.at("T0");
    auto T0 = std::dynamic_pointer_cast< ConstantRandomVariable >( AbstractT0 );
    T0->SetMean(StochSufficientStatistics[5][0]); 
    
    auto AbstractV0 = m_PopulationRandomVariables.at("V0");
    auto V0 = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractV0 );
    V0->SetMean(StochSufficientStatistics[6][0]); 
    
    /// Update Ksi
    double VarianceKsi = 0;
    for(auto it : StochSufficientStatistics[2])
    {
        VarianceKsi += it;
    }
    VarianceKsi /= NumberOfSubject;
    auto AbstractKsi = m_IndividualRandomVariables.at("Ksi");
    auto Ksi = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractKsi );
    Ksi->SetVariance(VarianceKsi);

    /// Update Tau
    double VarianceTau = 0;
    //std::cout << "Tau & Tau^2 : "; 
    for(auto it : StochSufficientStatistics[3])
    {
        VarianceTau += it;
        //std::cout << sqrt(it) << "/" << it << " ## ";
    }
    VarianceTau /= NumberOfSubject;
    //std::cout << std::endl;
    auto AbstractTau = m_IndividualRandomVariables.at("Tau");
    auto Tau = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractTau );
    Tau->SetVariance(VarianceTau);  

    
    /// Update Uncertainty variance

    /// Initialization
    double K = 0;

    /// Sum Yijk²
    double NoiseVariance = 0;
    for(auto it : *D)
    {
        K += it.size();
        for(auto it2 : it)
        {
            for(auto it3 : it2.first)
            {
                NoiseVariance += it3 * it3;
            }
        }
    }

    /// Sum -2 S1
    for(auto it : StochSufficientStatistics[0])
    {
        NoiseVariance -= 2* it;
    }


    /// Sum S2
    for(auto it : StochSufficientStatistics[1])
    {
        NoiseVariance += it;
    }


    // Divide by N*K, then take the square root
    NoiseVariance /= K;

    m_Noise->SetVariance(NoiseVariance);
    
    

    ComputeOutputs();

}


double 
UnivariateModel
::ComputeLikelihood(const std::shared_ptr<Realizations> &R, const std::shared_ptr<Data> &D,
                    const std::pair<std::string, int> NameRandomVariable) 
{
    double LogLikelihood = ComputeLogLikelihood(R, D, NameRandomVariable);
    return exp(LogLikelihood);
}

double 
UnivariateModel
::ComputeLogLikelihood(const std::shared_ptr<Realizations> &R, const std::shared_ptr<Data> &D,
                       const std::pair<std::string, int> NameRandomVariable) 
{
    /// Get the name of the realization, and its number (in case it is subject specific)
    std::string Name = NameRandomVariable.first.substr(0, NameRandomVariable.first.find_first_of("#"));
    int SubjectNumber = NameRandomVariable.second;
    
    bool PreviousEqualCurrentRealizations = (*R == std::get<2>(m_LastLogLikelihood));
    bool CurrentIsGeneric = !(Name == "Tau" or Name == "Ksi");
    bool PreviousIsGeneric = std::get<0>(m_LastLogLikelihood);
    
    /// Compute LogLikelihood
    double LogLikelihood;
    if(PreviousEqualCurrentRealizations && CurrentIsGeneric && PreviousIsGeneric)
    {
        LogLikelihood = std::get<1>(m_LastLogLikelihood);
    }
    else if(!CurrentIsGeneric)
    {
        LogLikelihood = ComputeLogLikelihoodIndividual(R, D, SubjectNumber);
    }
    else
    {
        LogLikelihood = ComputeLogLikelihoodGeneric(R, D);
    }
    
    m_LastLogLikelihood = std::tuple<bool, double, Realizations>(CurrentIsGeneric, LogLikelihood, *R);
    return LogLikelihood;
}

std::vector< std::vector< std::pair< std::vector<double>, double> > >
UnivariateModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    /// Simulate Realizations
    auto R = std::make_shared<Realizations>( SimulateRealizations(NumberOfSubjects) );
    
    /// Initialize
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::normal_distribution<double> ObsDistrib(70.0, sqrt(3.0));
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    Data D;
    double T0 = R->at("T0")[0];
    double P0 = R->at("P0")[0];
    double V0 = R->at("V0")[0];
    
    /// Simulate the data
    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        IndividualData InDa;
        
        /// Generate time points
        std::vector<double> TimePoints;
        for(int j = 0; j < Uni(RNG); ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));   
        }
        std::sort(TimePoints.begin(), TimePoints.end());
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        
        /// Generate observations corresponding to the time points
        for(auto it : TimePoints)
        {
            double TimePoint = SubjectTimePoint(it);
            double Noise = NoiseDistrib(RNG);
            double ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
            auto Scores = { ParallelCurve + Noise};
            
            InDa.push_back(std::pair<std::vector<double>, double> (Scores, it));
        }
        
        D.push_back(InDa);
    }
    
    return D;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
UnivariateModel
::InitializeFakeRandomVariables() 
{
 
    /// Initial Position
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.3, 0.01 );
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    
    /// Initial Time
    auto T0 = std::make_shared<ConstantRandomVariable>( 70.0 );
    m_PopulationRandomVariables.insert( RandomVariable("T0", T0));
    
    /// Initial Velocity
    auto V0 = std::make_shared<GaussianRandomVariable>( 0.06, 0.001 );
    m_PopulationRandomVariables.insert( RandomVariable("V0", V0) );
    
    /// Noise
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.05);
    
    /// Individual Time Shift
    auto Tau = std::make_shared< GaussianRandomVariable >(0.0, 3.0);
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    
    /// Individual pace / Preaceleration factor
    auto Ksi = std::make_shared< GaussianRandomVariable >(0.0, 1.5); 
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
    
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::function<double(double)>  
UnivariateModel
::GetSubjectTimePoint(const int SubjectNumber, const std::shared_ptr<Realizations>& R)
{
    double AccFactor = exp(R->at("Ksi")[SubjectNumber]);
    double TimeShift = R->at("Tau")[SubjectNumber];
    double T0 = R->at("T0")[0];
    
    return [AccFactor, TimeShift, T0](double t) { return AccFactor * (t - TimeShift - T0) + T0; };
}


double 
UnivariateModel
::ComputeLogLikelihoodGeneric(const std::shared_ptr<Realizations> &R,
                              const std::shared_ptr<Data> &D) 
{
    /// Get the data
    double T0 = R->at("T0")[0];
    double P0 = R->at("P0")[0];
    double V0 = R->at("V0")[0];
    
    /// Compute the loglikelihood
    double LogLikelihood = 0, K = 0;
    int i = 0;
    for(auto IterData = D->begin(); IterData != D->end(); ++i, ++IterData)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        K += IterData->size();
        
        for(auto it : *IterData)
        {
            double TimePoint = SubjectTimePoint(it.second);
            auto ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
            double Norm = (it.first[0] - ParallelCurve);
            LogLikelihood += Norm * Norm ;
        }
    }
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    
    LogLikelihood -= K * log(sqrt( 2 * M_PI * m_Noise->GetVariance() ));
    
    return LogLikelihood;
}


double
UnivariateModel
::ComputeLogLikelihoodIndividual(const std::shared_ptr<Realizations> &R,
                                 const std::shared_ptr<Data> &D, const int SubjectNumber) 
{
    /// Initialize the individual parameters
    double T0 = R->at("T0")[0];
    double P0 = R->at("P0")[0];
    double V0 = R->at("V0")[0];
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    
    /// Compute the loglikelihood
    double LogLikelihood = 0;
    double k = D->at(SubjectNumber).size();
    int i = 0;
    for(auto IterData = D->at(SubjectNumber).begin(); IterData != D->at(SubjectNumber).end(); ++IterData, ++i)
    {
        double TimePoint = SubjectTimePoint(IterData->second);
        auto ParallelCurve = m_BaseManifold->ComputeParallelCurve(P0, T0, V0, 0.0, TimePoint);
        double Norm = (IterData->first[0] - ParallelCurve) * (IterData->first[0] - ParallelCurve);
        LogLikelihood += Norm*Norm;
    }

    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= k * log(sqrt( 2 * m_Noise->GetVariance() * M_PI));
    return LogLikelihood;
    
}


void
UnivariateModel
::ComputeOutputs() 
{
    double MeanP0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"))->GetMean();
    double MeanT0 = std::dynamic_pointer_cast<ConstantRandomVariable>(m_PopulationRandomVariables.at("T0"))->GetMean();
    double MeanV0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("V0"))->GetMean();
    double SigmaKsi = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetVariance();
    double SigmaTau = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetVariance();
    double Sigma = m_Noise->GetVariance();

    m_OutputParameters << MeanP0 << ", " << MeanT0 << ", " << MeanV0 << ", ";
    m_OutputParameters << SigmaKsi << ", " << SigmaTau << ", " << Sigma << ", " << std::endl;

    
    std::cout << "Parameters : P0: " << MeanP0 << ". T0: " << MeanT0 << ". V0: " << MeanV0 << ". Ksi: " << SigmaKsi << ". Tau: " << SigmaTau << ". Sigma: " << Sigma << ". " << std::endl;
}






