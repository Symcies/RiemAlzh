#include "TestModel.h"
#include "../RandomVariables/GaussianRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


TestModel
::TestModel() 
{
    
}


TestModel
::~TestModel()
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
TestModel
::Initialize() 
{
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto A = std::make_shared<GaussianRandomVariable>(0, 1.1);
    auto B = std::make_shared<GaussianRandomVariable>(1.2, 0.3);
    auto C = std::make_shared<GaussianRandomVariable>(-1.2, 0.00001);
    m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.001);
    
    m_IndividualRandomVariables.insert(RandomVariable("A", A));
    m_IndividualRandomVariables.insert(RandomVariable("B", B));
    m_PopulationRandomVariables.insert(RandomVariable("C", C));
}


void 
TestModel
::UpdateParameters(const std::shared_ptr <MultiRealizations> &R, const std::vector<std::string> Names) 
{
    // TODO : Check if something has to be added
}


std::map< std::string, double >
TestModel
::GetParameters() 
{
    auto A = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("A"));
    auto B = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("B"));
    auto C = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("C"));
    std::map<std::string, double> Parameters;
    Parameters["A_Mean"] = A->GetMean();
    Parameters["A_Var"] = A->GetVariance();
    Parameters["B_Mean"] = B->GetMean();
    Parameters["B_Var"] = B->GetVariance();
    Parameters["C_Mean"] = C->GetMean();
    Parameters["C_Var"] = C->GetVariance();
    return Parameters;
}

TestModel::SufficientStatisticsVector
TestModel
::GetSufficientStatistics(const std::shared_ptr<MultiRealizations> &R,
                          const std::shared_ptr<Data> &D) 
{
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
    auto IterS1 = S1.begin();
    auto IterS2 = S2.begin();
    int i = 0;
    for(auto IterData = D->begin(); IterData != D->end() &&  IterS1 != S1.end() && IterS2 != S2.end(); ++i, ++IterData) 
    {
        for(auto IterIndivData = IterData->begin(); IterIndivData != IterData->end(); ++IterIndivData)
        {
            double TimePoint = IterIndivData->second;
            double at = GetA(i, R) * TimePoint + GetB(i, R) + GetC(R) * TimePoint * TimePoint;
            double y = IterIndivData->first[0];
            
            *IterS1 = at * y;
            *IterS2 = at * at;
            ++IterS1, ++IterS2;
        }
    }
    
    /// Compute S3 and S4
    VectorType S3(NumberOfSubjects), S4(NumberOfSubjects);
    auto IterS3 = S3.begin();
    auto IterS4 = S4.begin();
    i = 0;
    for(    ; IterS3 != S3.end() && IterS4 != S4.end(); ++IterS3, ++IterS4, ++i)
    {
        double a = GetA(i, R);
        *IterS3 = a;
        *IterS4 = a*a;
    }
    
    /// Compute S5 and S6
    VectorType S5(NumberOfSubjects), S6(NumberOfSubjects);
    auto IterS5 = S5.begin();
    auto IterS6 = S6.begin();
    i = 0;
    for(    ; IterS5 != S5.end() && IterS6 != S6.end(); ++IterS5, ++IterS6, ++i)
    {
        double b = GetB(i, R);
        *IterS5 = b;
        *IterS6 = b*b;
    }
    
    
    /// Compute S7 
    VectorType S7(1, R->at("C")(0));
    
    
    SufficientStatisticsVector S = {S0, S1, S2, S3, S4, S5, S6, S7};
    
    /// Test
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
    
    TestAssert::WarningEquality_Function(f1, f2, "Likelihood not correct. TestModel > GetSufficientStatistics");
    
    return S;
}


void 
TestModel
::UpdateRandomVariables(const SufficientStatisticsVector &StochSufficientStatistics,
                        const std::shared_ptr <Data> &D) 
{
    double NumberOfSubjects = D->size();
    
    /// Update a mean and variance
    double AMean = 0, AVariance = 0;
    for(auto IterA = StochSufficientStatistics[3].begin(); IterA != StochSufficientStatistics[3].end(); ++IterA)
    {
        AMean += *IterA;
    }
    AMean /= NumberOfSubjects;
    
    for(auto IterA_ = StochSufficientStatistics[4].begin(); IterA_ != StochSufficientStatistics[4].end(); ++IterA_)
    {
        AVariance += *IterA_;
    }
    AVariance -= NumberOfSubjects * AMean * AMean;
    AVariance /= NumberOfSubjects;
    
    auto AbstractA = m_IndividualRandomVariables.at("A");
    auto A = std::static_pointer_cast<GaussianRandomVariable>(AbstractA);
    
    A->SetMean(AMean);
    A->SetVariance(AVariance);
    
    /// Update b mean and variance
    double BMean = 0, BVariance = 0;
    for(auto IterB = StochSufficientStatistics[5].begin(); IterB != StochSufficientStatistics[5].end(); ++IterB)
    {
        BMean += *IterB;
    }
    BMean /= NumberOfSubjects;
    
    for(auto IterB_ = StochSufficientStatistics[6].begin(); IterB_ != StochSufficientStatistics[6].end(); ++IterB_)
    {
        BVariance += *IterB_;
    }
    BVariance -= NumberOfSubjects * BMean * BMean;
    BVariance /= NumberOfSubjects;
    
    auto AbstractB = m_IndividualRandomVariables.at("B");
    auto B = std::static_pointer_cast<GaussianRandomVariable>(AbstractB);
    
    B->SetMean(BMean);
    B->SetVariance(BVariance);
    
    
    /// Update c mean
    double CMean = StochSufficientStatistics[7](0);
    
    auto AbstractC = m_PopulationRandomVariables.at("C");
    auto C = std::static_pointer_cast<GaussianRandomVariable>(AbstractC);
    
    C->SetMean(CMean);
    
    /// Update noise
    double K = 0.0;
    for(const auto& it : *D)
    {
        K += it.size();
    }
    
    double NoiseVariance = StochSufficientStatistics[0](0);
    
    auto IterS1 = StochSufficientStatistics[1].begin();
    auto IterS2 = StochSufficientStatistics[2].begin();
    for(    ; IterS1 != StochSufficientStatistics[1].end() && IterS2 != StochSufficientStatistics[2].end(); ++IterS1, ++IterS2)
    {
        NoiseVariance += -2 * *IterS1 + *IterS2;
    }
    
    NoiseVariance /= K;
    m_Noise->SetVariance(NoiseVariance);
}


double
TestModel
::ComputeLogLikelihood(const std::shared_ptr <MultiRealizations> &R,
                       const std::shared_ptr <Data> &D) 
{
    double LogLikelihood = 0.0, K = 0.0;
    int i = 0;
    for(auto IterData = D->begin(); IterData != D->end(); ++IterData, ++i)
    {
        K += IterData->size();
        for(auto IterIndivData = IterData->begin(); IterIndivData != IterData->end(); ++IterIndivData)
        {
            double y = IterIndivData->first(0);
            double t = IterIndivData->second ;
            double at = GetA(i, R) * t + GetB(i, R) + GetC(R)*t*t;
            double Norm = (y - at);
            LogLikelihood += Norm * Norm;
        }
    }
    
    LogLikelihood /= -2 * m_Noise->GetVariance();
    
    LogLikelihood -= K * log(sqrt( 2 * M_PI * m_Noise->GetVariance()));
    
    return LogLikelihood;
}

double 
TestModel
::ComputeIndividualLogLikelihood(const std::shared_ptr <MultiRealizations> &R,
                                 const std::shared_ptr <Data> &D, const int SubjectNumber) 
{
    double LogLikelihood = 0.0;
    for(auto IterData = D->at(SubjectNumber).begin(); IterData != D->at(SubjectNumber).end(); ++IterData)
    {
        double y = IterData->first(0);
        double t = IterData->second;
        double at = GetA(SubjectNumber, R) * t + GetB(SubjectNumber, R) + GetC(R) * t * t;
        double Norm = y - at;
        LogLikelihood += Norm * Norm;
    }
    double K = D->at(SubjectNumber).size();
    LogLikelihood /= -2 * m_Noise->GetVariance();
    LogLikelihood -= K *log( sqrt( 2 * M_PI * m_Noise->GetVariance() ));
    
    return LogLikelihood;
}

TestModel::Data
TestModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    auto R = std::make_shared<MultiRealizations>( SimulateRealizations(NumberOfSubjects));
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::uniform_real_distribution<double> TimeObs(-5, 5);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    
    Data D;
    double q = 0.0, SumNoise = 0.0;
    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        int NbObservations = Uni(RNG);
        IndividualData InDa;
        std::vector<double> TimePoints;
        for(int j = 0; j < NbObservations; ++j)
        {
            TimePoints.push_back(TimeObs(RNG));
        }
        std::sort(TimePoints.begin(), TimePoints.end());
        
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it)
        {
            double at = GetA(i, R) * *it + GetB(i, R) + GetC(R) * *it * *it;
            double Noise = NoiseDistrib(RNG);
            SumNoise += Noise;
            q += 1;
            InDa.push_back(std::pair<VectorType, double>( VectorType(1, at + Noise), *it));
        }
        D.push_back(InDa);
    }
   
    
    double P = ComputeLogLikelihood(R, std::make_shared<Data>(D)) ; 
    std::cout << "Real Noise : " << SumNoise/q << std::endl;
    std::cout << "Likelihood : " << P << std::endl;
    
    return D;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
TestModel
::ComputeOutputs()
{
    auto A = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("A"));
    auto B = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("B"));
    auto C = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("C"));
    double AMean = A->GetMean();
    double AVar = A->GetVariance();
    double BMean = B->GetMean();
    double BVar = B->GetVariance();
    double CMean = C->GetMean();
    double Sigma = m_Noise->GetVariance();
    
    std::cout << "Parameters : AMean " << AMean << ". AVar : " << AVar << ". BMean : " << BMean;
    std::cout << ". BVar : " << BVar << ". CMean : " << CMean <<  ". Sigma : " << Sigma << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
TestModel
::InitializeFakeRandomVariables()
{
    auto A = std::make_shared<GaussianRandomVariable>(2, 0.2);
    auto B = std::make_shared<GaussianRandomVariable>(-2, 0.4);
    auto C = std::make_shared<GaussianRandomVariable>(0.5, 0.00001);
    m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.00001);
    
    m_IndividualRandomVariables.insert(RandomVariable("A", A));
    m_IndividualRandomVariables.insert(RandomVariable("B", B));
    m_PopulationRandomVariables.insert(RandomVariable("C", C));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double 
TestModel
::GetA(int i, const std::shared_ptr<MultiRealizations>& R) 
{
    return R->at("A")(i);
}

double 
TestModel
::GetB(int i, const std::shared_ptr<MultiRealizations>& R) 
{
    return R->at("B")(i);
}


double 
TestModel
::GetC(const std::shared_ptr<MultiRealizations> R)
{
    return R->at("C")(0);
}