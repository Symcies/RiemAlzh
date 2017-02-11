#include <stdexcept>
#include "LongitudinalModel.h"
#include "PropagationManifold.h"
#include "TestAssert.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LongitudinalModel
::LongitudinalModel(const unsigned int NbIndependentComponents, std::shared_ptr<AbstractManifold>& M)
{
    m_NbIndependentComponents = NbIndependentComponents;
    std::shared_ptr<PropagationManifold> Manifold = std::dynamic_pointer_cast<PropagationManifold>(M);
    // m_Manifold = Manifold;
    m_ManifoldDimension = (int)Manifold->GetDimension();
}

LongitudinalModel
::~LongitudinalModel()
{ 
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::Initialize(const Data& D)
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto P0 = std::make_shared<GaussianRandomVariable>(0.3, 0.0001 * 0.0001);
    auto Ksi = std::make_shared<GaussianRandomVariable>(-4.6, 0.4);
    auto Tau = std::make_shared<GaussianRandomVariable>(65, 1.0);
    m_Noise = std::make_shared<GaussianRandomVariable>(0.0, 0.0000001);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0) );
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    
    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension-1) ; ++i)
    {
        auto Beta = std::make_shared<GaussianRandomVariable>(1 , 0.001 * 0.001);
        m_PopulationRandomVariables.insert(RandomVariable("Beta#" + std::to_string(i), Beta));
    }
    
    for(int i = 1; i < m_ManifoldDimension; ++i)
    {
        auto Delta = std::make_shared<GaussianRandomVariable>(0, 0.001 * 0.001);
        m_PopulationRandomVariables.insert( RandomVariable("Delta#" + std::to_string(i), Delta) );
    }

    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    { 
        auto S = std::make_shared<GaussianRandomVariable>(0.0, 1.0);
        m_IndividualRandomVariables.insert(RandomVariable("S#" + std::to_string(i), S));
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
    
    /// Tests
    TestAssert::WarningInequality_GreaterThan(1.0, P0->GetMean(), 0.0, "LongitudinalModel>Initialize : wrong P0");
    TestAssert::WarningInequality_GreaterThan(m_Noise->GetVariance(), 0.0, "LongitudinalModel>Initialize : wrong noise");
}

void 
LongitudinalModel
::UpdateModel(const Realizations &R, int Type,
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
        else if(Name == "P0" or Name == "Delta" or Name == "All")
        {
            UpdateCase = 4;
            break;
        }
        else
        {
            UpdateCase = 4;
            std::cout << "Should be" << Name << "be in LongitudinalModel > Update Parameters?" << std::endl;
            break;
        }
    }
    
    /// Update Case in fonction of the the vector of names
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
        default:
            std::cout << "Error? LongitudinalModel > UpdateModel";
            break;
    }
}


LongitudinalModel::SufficientStatisticsVector
LongitudinalModel
::GetSufficientStatistics(const Realizations& R, 
                          const Data& D)
{
    /////////////////////////
    /// Initialization
    /////////////////////////
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    double NumberOfSubjects = R.at("Ksi").size();
    
    
    /////////////////////////
    /// Compute S1 and S2
    /////////////////////////
    VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
    int i = 0;
    auto IterS1 = S1.begin(), IterS2 = S2.begin();
    for(auto Iter = D.begin(); Iter != D.end(); ++Iter, ++i)
    {
        /// Given a particular subject, get its attributes, then, loop over its observation
        auto TimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        for(auto it : *Iter)
        {
            double Time = TimePoint(it.second);
            VectorType ParallelCurve = ComputeGeodesic(R.at("P0")(0), Time, PropagationCoefficients, SpaceShift);

            *IterS1 = dot_product(ParallelCurve, it.first);
            *IterS2 = ParallelCurve.squared_magnitude();
            
            ++IterS1, ++IterS2;
        }

    }

    /////////////////////////
    /// Compute S3 and S4
    /////////////////////////
    VectorType S3(NumberOfSubjects), S4(NumberOfSubjects);
    auto itKsi = R.at("Ksi").begin();
    auto itS3 = S3.begin(), itS4 = S4.begin();
    for(    ; itKsi != R.at("Ksi").end() ; ++itKsi, ++itS3, ++itS4)
    {
        *itS3 = *itKsi;
        *itS4 = *itKsi * *itKsi;
    }
    
    /////////////////////////
    /// Compute S5 and S6
    /////////////////////////
    VectorType S5(NumberOfSubjects), S6(NumberOfSubjects);
    auto itTau = R.at("Tau").begin();
    auto itS5 = S5.begin(), itS6 = S6.begin();
    for(    ; itTau != R.at("Tau").end(); ++itTau, ++itS5, ++itS6)
    {
        *itS5 = *itTau;
        *itS6 = *itTau * *itTau;
    }
    
    /////////////////////////
    /// Compute S7
    /////////////////////////
    VectorType S7(1, R.at("P0")(0));
    
    /////////////////////////
    /// Compute S8
    /////////////////////////
    VectorType S8((m_ManifoldDimension - 1)*m_NbIndependentComponents);
    i = 0;
    for(auto it = S8.begin(); it != S8.end(); ++it, ++i)
    {
        *it = R.at("Beta#" + std::to_string(i))(0);
    }
    
    /////////////////////////
    /// Compute S9
    /////////////////////////
    VectorType S9(m_ManifoldDimension - 1);
    i = 1;
    for(auto it = S9.begin(); it != S9.end(); ++it, ++i)
    {
        *it = R.at("Delta#" + std::to_string(i))(0);
    }
    
    SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9};
    
    /*
    /////////////////////////
    /// Tests
    /////////////////////////
    /// There are two ways to compute the logLikelihood. Check if they are equal
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
        Calculation -= m_NbTotalOfObservations * log(sqrt( 2 * M_PI * m_Noise->GetVariance() ));  
            
        return Calculation;
    };
    
    TestAssert::WarningEquality_Function(f1, f2, "Likelihood not correct. LongitudinalModel > GetSufficientStatistics");
    */
     
    return S;
    
}


void
LongitudinalModel
::UpdateRandomVariables(const SufficientStatisticsVector& SS, const Data& D)
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
}


double
LongitudinalModel
::ComputeLogLikelihood(const Realizations& R, const Data& D) 
{
    /*
    /// Get the data
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = GetInitialTime();
    VectorType P0(1, R.at("P0")(0) );
    VectorType V0(1, R.at("V0")(0) );
    VectorType Delta = GetPropagationCoefficients(R);
    
    /// Compute the likelihood
    double LogLikelihood = 0, K = 0;
    int i = 0;
    for(auto IterData = D.begin(); IterData != D.end(); ++IterData, ++i)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        K += IterData->size();
        
        for(auto it = IterData->begin(); it != IterData->end(); ++it)
        {
            double TimePoint = SubjectTimePoint(it->second);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Delta);
            LogLikelihood  += (it->first - ParallelCurve).squared_magnitude();
        }
    }
    /// Tests
    TestAssert::WarningInequality_GreaterThan(LogLikelihood, 0.0, "LongitudinalModel>ComputeIndividualLogLikelihood ; wrong Likelihood");
    
    LogLikelihood  /= -2*m_Noise->GetVariance();
    LogLikelihood -= K*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
    */
    double a = 0;
    return 0 ;
}

double
LongitudinalModel
::ComputeIndividualLogLikelihood(const Realizations& R, const Data& D, const int SubjectNumber) 
{
    // TODO : send only D(i) to the function
    
    /// Get the data
    double P0 = R.at("P0")(0);
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(SubjectNumber));
    VectorType PropagationCoefficients = GetPropagationCoefficients(R);
    
    /// Compute the likelihood
    double LogLikelihood = 0;
    double N = D.at(SubjectNumber).size();

#pragma omp parallel for reduction(+:LogLikelihood)
    for(size_t i = 0; i < N; ++i)
    {
        auto& it = D.at(SubjectNumber).at(i);
        double TimePoint = SubjectTimePoint(it.second);
        VectorType ParallelCurve = ComputeGeodesic(P0, TimePoint, PropagationCoefficients, SpaceShift);
        LogLikelihood += (it.first - ParallelCurve).squared_magnitude();
    }
    
    TestAssert::WarningInequality_GreaterThan(LogLikelihood, 0.0, "LongitudinalModel>ComputeIndividualLogLikelihood ; wrong Likelihood");
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= N * log( 2 * m_Noise->GetVariance() * M_PI) / 2.0;

    return LogLikelihood;
}


LongitudinalModel::Data
LongitudinalModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs)
{
    /*
    /// Simulate realizations
    auto R = std::make_shared<Realizations>( SimulateRealizations(NumberOfSubjects) );
    
    /// Initialize model attributes
    ComputeOrthonormalBasis(*R);
    ComputeAMatrix(*R);
    ComputeSpaceShifts(*R);

    /// Initialize
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::normal_distribution<double> ObsDistrib(70.0, sqrt(3.0));
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = GetInitialTime();
    VectorType P0(1, R->at("P0")(0) );
    VectorType V0(1, R->at("V0")(0) );
    VectorType Delta = GetPropagationCoefficients(*R);
    
    Data D;
    
    /// Simulate the data
    for(int i = 0; i < NumberOfSubjects; ++i) 
    {
        IndividualData InDa;
        
        /// Generate the time points of the subjets
        std::vector<double> TimePoints;
        for(int j = 0; j < Uni(RNG); ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));
        }
        std::sort(TimePoints.begin(), TimePoints.end());
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, *R);
        auto SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        /// Generate observations corresponding to the time points
        for(auto it : TimePoints)
        {
            VectorType Scores(SpaceShift.size());
            
            double TimePoint = SubjectTimePoint(it);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Delta);
            int i = 0;
            for(auto it2 = ParallelCurve.begin(); it2 != ParallelCurve.end(); ++it2, ++i)
            {
                Scores(i) = *it2 + NoiseDistrib(RNG);
            }
            InDa.push_back(std::pair<VectorType, double> (Scores, it));
        }
        
        D.push_back(InDa);
    }
    */
    Data D;
    return D;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::DisplayOutputs(const Realizations& R)
{
    auto P0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"));
    auto Tau = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    auto Ksi = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    double Sigma = m_Noise->GetVariance();
    
    std::cout << "Parameters : P0: " << P0->GetMean() << ". T0: " << Tau->GetMean() << ". Var(Tau): " << Tau->GetVariance();
    std::cout << ". V0: " << exp(Ksi->GetMean()) << ". Sigma: " << Sigma << std::endl;

}


void
LongitudinalModel
::SaveData(unsigned int IterationNumber, const Realizations &R) 
{
    unsigned int NumberOfSubjects = R.at("Tau").size();
    std::ofstream Outputs;
    std::string FileName = "Outputs/MultivariateModel_Parameters_Iteration" + std::to_string(IterationNumber) + ".txt";
    //std::string FileName = "Outputs/Multivariate/Parameters.txt";
    Outputs.open(FileName, std::ofstream::out | std::ofstream::trunc);
    
    /// Save Number of subject, Dimension and Number of Sources
    Outputs << NumberOfSubjects << ", " << m_ManifoldDimension << ", " << m_NbIndependentComponents << std::endl;
    
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
    for(size_t i = 1; i < m_ManifoldDimension; ++i)
    {
        Outputs << R.at("Delta#" + std::to_string(i))(0);
        if(i != m_ManifoldDimension - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Delta_k)
    Outputs << 0 << ", ";
    for(size_t i = 1; i < m_ManifoldDimension; ++i)
    {
        std::string Name = "Delta#" + std::to_string(i);
        auto Delta = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at(Name));
        Outputs << Delta->GetMean();
        if(i != m_ManifoldDimension - 1) { Outputs << ", "; }
    }
    Outputs << std::endl;
    
    /// Save (Beta_k)
    for(size_t i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension-1); ++i) 
    {
        Outputs << R.at("Beta#" + std::to_string(i))(0);
        if(i != m_NbIndependentComponents*(m_ManifoldDimension-1) - 1) { Outputs << ", "; }
    }
    
    /// Save (Beta_mean_k)
    for(size_t i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension-1); ++i)
    {
        std::string Name = "Beta#" + std::to_string(i);
        auto Beta = std::static_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at(Name));
        Outputs << Beta->GetMean();
        if(i != m_NbIndependentComponents*(m_ManifoldDimension-1) - 1) { Outputs << ", "; }
    }
    
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
    for(size_t i = 0; i < NumberOfSubjects; ++i)
    {
        VectorType W = m_SpaceShifts.at("W" + std::to_string(i));
        for(auto it = W.begin(); it != W.end(); ++it)
        {
            Outputs << *it;
            if(i != W.size() - 1) { Outputs << ", "; }
        }
        Outputs << std::endl;
    }
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::InitializeFakeRandomVariables()
{
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    /////////////////////////////
    /// Population Parameters ///
    /////////////////////////////
    auto P0 = std::make_shared<GaussianRandomVariable>(0.35, 0.000001);
    auto V0 = std::make_shared<GaussianRandomVariable>(0.06, 0.000001);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    m_PopulationRandomVariables.insert(RandomVariable("V0", V0));

    m_Noise = std::make_shared< GaussianRandomVariable >(0.0, 0.0001);

    for(int i = 0; i < m_NbIndependentComponents*(m_ManifoldDimension-1) ; ++i)
    {
        auto Beta = std::make_shared< GaussianRandomVariable> ((double)i/5.0, 0.0001);
        std::string name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(name, Beta));
    }
    
    for(int i = 1; i < m_ManifoldDimension ; ++i)
    {
        auto Delta = std::make_shared< GaussianRandomVariable >(0.5 + (double)i, 0.0001);
        m_PopulationRandomVariables.insert(RandomVariable("Delta#" + std::to_string(i), Delta));
    }


    /////////////////////////////
    /// Individual Parameters ///
    /////////////////////////////
    double T0 = 70.0;
    auto Ksi = std::make_shared< GaussianRandomVariable >(0.0, 0.5*0.5);
    auto Tau = std::make_shared< GaussianRandomVariable >(T0, 3.0*3.0);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        auto S = std::make_shared< GaussianRandomVariable >(0.0, 1.0/2.0);
        m_IndividualRandomVariables.insert(RandomVariable("S#" + std::to_string(i), S));
    }

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double 
LongitudinalModel
::GetInitialTime() 
{
    auto T0 = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"));
    return T0->GetMean();
}


LongitudinalModel::VectorType
LongitudinalModel
::GetPropagationCoefficients(const Realizations& R)
{
    VectorType Delta(m_ManifoldDimension, 0.0);

    int i = 1;
    for(auto it = Delta.begin() + 1; it != Delta.end(); ++it, ++i)
    {
        *it = R.at("Delta#" + std::to_string(i))(0);
    }

    return Delta;
}

std::function<double(double)>  
LongitudinalModel
::GetSubjectTimePoint(const int SubjectNumber, const Realizations& R)
{
    double AccFactor = exp(R.at("Ksi")(SubjectNumber));
    double TimeShift = R.at("Tau")(SubjectNumber);
    
    return [AccFactor, TimeShift](double t) { return AccFactor * (t - TimeShift); };
}

void
LongitudinalModel
::ComputeOrthonormalBasis(const Realizations& R)
{
    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    auto Ksi = std::static_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"));
    double V0 = exp(Ksi->GetMean());
    double T0 = GetInitialTime();
    double P0 = R.at("P0")(0);
    VectorType Delta = GetPropagationCoefficients(R);
    
    /// Compute the transformation to do the Householder reflection in a Euclidean space
    VectorType U = ComputeGeodesicTransformation(P0, T0, V0, Delta);
    
    /// Compute the initial pivot vector U
    double Norm = U.magnitude();
    U(0) += copysign(1, -U(0))*Norm;
    double NormU = U.squared_magnitude();

    /// Compute Q = I(N) - 2U . Ut / (Ut . U)
    std::vector<VectorType> Q;
    for(auto it = U.begin(); it != U.end(); ++it)
    {
        VectorType Coordinate(U.size());
        int i = 0;
        for(auto it2 = U.begin(); it2 != U.end(); ++it2, ++i)
        {
            Coordinate(i) =  - 2 * *it * *it2 / NormU;
        }
        Q.push_back( Coordinate );
    }

    for(int i = 0; i < Q.size() ; ++i)
    {
        Q[i](i) += 1;
    }
    

    /*
    /// TESTS
    
    /// Test colinearity between the first vector and Q[0]
    
    
    /// Orthogonality between the basis vectors and the velocity
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    for(auto it = Q.begin() + 1; it != Q.end(); ++it)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(*it, InitialVelocity, InitialPosition);
        };
        std::function<double()> f2 = []() { return 0;};
        TestAssert::WarningEquality_Function(f1, f2, "Basis vector B not orthogonal to the velocity. Longotudinal Model > ComputeOrthonormalBasis");
    }   
    
    /// END TEST
    
    /// Drop the first vector which is colinear to gamma_derivative(t0)
    */
    Q.erase(Q.begin());
    m_OrthogonalBasis = Q;
}

void
LongitudinalModel
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
    
    m_AMatrix = AMatrix;
    
    
    /*
    /// TESTS
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    
    /// Check if the basis is orthogonal to the velocity
    for(int i = 0; i < m_AMatrix.columns(); ++i)
    {
        auto it = m_AMatrix.get_column(i);
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        std::function<double()> f2 = []() { return 0;};
        TestAssert::WarningEquality_Function(f1, f2, "Basis vector B not orthogonal to the velocity. LongitudinalModel > ComputeAMatrix");
    }
    
    /// Check if the matrix column is orthogonal to the velocity
    for(auto it : m_OrthogonalBasis)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        std::function<double()> f2 = []() { return 0;};
        TestAssert::WarningEquality_Function(f1, f2, "A column not orthogonal to the velocity. LongitudinalModel > ComputeAMatrix");
    }
    /// END TESTS
    */
}

void
LongitudinalModel
::ComputeSpaceShifts(const Realizations& R)
{
    std::map< std::string, VectorType> SpaceShifts;
    int NumberOfSubjects = (int)R.at("Tau").size();
    if(NumberOfSubjects != R.at("Ksi").size())
    {
        throw std::invalid_argument("Not the same number of realization in Tau and in Ksi. Whereas its equal to the number of subjects");
    }

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


    ///////////////////////////////
    //// Debugging : Unit tests ///
    ///////////////////////////////
    
    /*
    /// Get Data
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = GetInitialTime();
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    auto Delta = GetPropagationCoefficients(R);
    std::function<double()> NullFunction = []() { return 0;};
    
    /// DEBUG 1 : The orthonormal basis is orthogonal to the velocity
    for(int i = 0; i < m_AMatrix.columns(); ++i)
    {
        auto it = m_AMatrix.get_column(i);
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        
        TestAssert::WarningEquality_Function(f1, NullFunction, "Basis vector B not orthogonal to the velocity. LongitudinalModel > ComputeSpaceShifts");
    }
    
    /// DEBUG 2 : Check if the matrix column is orthogonal to the velocity
    for(auto it : m_OrthogonalBasis)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        TestAssert::WarningEquality_Function(f1, NullFunction, "A column not orthogonal to the velocity. LongitudinalModel > ComputeSpaceShifts");
    }
    
    /// DEBUG 3 : Check if the space shifts are orthogonal to the velocity
    for(auto it : m_SpaceShifts)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it.second, InitialVelocity, InitialPosition);
        };
        TestAssert::WarningEquality_Function(f1, NullFunction, "Space Shifts not orthogonal to the velocity. Longitudinal Model > ComputeSpaceShifts");
    }
    
    /// DEBUG 4 : <diff(g(t)) , ParallelTransport>|g(t) = 0
    /// DEBUG 5 : <ParallelTransport, ParallelTransport>g(t) = Constant
    auto ParallelTransport = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, m_SpaceShifts["W0"], 60.0, Delta);
    auto position = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, 60.0, Delta);
    double Constant = m_Manifold->ComputeScalarProduct(ParallelTransport, ParallelTransport, position);
    std::function<double()> ConstantFunction = [Constant] () { return Constant; };
    
    
    for(double t = 61; t < 80; t += 0.5)
    {
        
        std::function<double()> f1 = [&CastedManifold, &InitialPosition, &InitialVelocity, &T0, &Delta, &SpaceShifts, t] () 
        {
            auto a = CastedManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, t, Delta);
            auto b = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, SpaceShifts["W0"], t, Delta);
            auto pos = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
            return CastedManifold->ComputeScalarProduct(a, b, pos);
        };
    
        std::function<double()> f2= [&SpaceShifts, &CastedManifold, &InitialPosition, &InitialVelocity, &T0, &Delta, t] () 
        {
            auto a = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, SpaceShifts["W0"], t, Delta);
            auto pos = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
            return CastedManifold->ComputeScalarProduct(a, a, pos);
        };
    
        
        
        
        TestAssert::WarningEquality_Function(f1, NullFunction, "Parallel transport not orthogonal to velocity at " + std::to_string(t) + ". LongitudinalModel > ComputeSpaceShifts");
        TestAssert::WarningEquality_Function(f2, ConstantFunction, "Norm of parallel transport not constant over time. LongiudinalModel > ComputeSpaceShifts");
    }
    
    
    /// DEBUG 6 : <ParallelTransport, ParallelTransport>g(t) = Constant
    
    for(double t = 62; t < 75; t += 0.5)
    {
        std::vector<double> Geodesic = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
        std::vector<double> GeodesicDerivative = CastedManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, t, Delta);
        std::vector<double> ParallelTransport = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, m_SpaceShifts["W0"], t, Delta);
        double ScalarProd = m_Manifold->ComputeScalarProduct(ParallelTransport, GeodesicDerivative, Geodesic);
        std::cout << "<diff(g(t)), ParallelTransport(W0)>|g(t)> at " << t << " : " << ScalarProd << std::endl;
        double Norm = m_Manifold->ComputeScalarProduct(ParallelTransport, ParallelTransport, Geodesic);
        std::cout << "Norm at " << t << " : " << Norm << std::endl;
    }
    */
    
    ///////////////////////////////
    ////End Unit tests ///
    ///////////////////////////////
    
}



LongitudinalModel::VectorType
LongitudinalModel
::ComputeGeodesic(double P0, double TimePoint, VectorType Delta, VectorType SpaceShift) 
{
    auto N = Delta.size();
    VectorType ParallelCurve(N);
    
    ScalarType * curve = ParallelCurve.memptr();
    ScalarType * d = Delta.memptr();
    ScalarType * s = SpaceShift.memptr();

#pragma omp simd
    for(size_t i = 0; i < N; ++i)
    {
        double Val = exp(-d[i] / (P0 * (1-P0)));
        Val = TimePoint + d[i] + s[i] * (P0 + (1.-P0)*Val) * (P0 + (1.-P0)*Val) / Val;
        Val = 1 + (1./P0 - 1.) * exp( - Val / (P0 * (1.-P0)));
        curve[i] = 1. /Val;
    }
    
    return ParallelCurve;
}

LongitudinalModel::VectorType
LongitudinalModel
::ComputeGeodesicTransformation(double P0, double T0, double V0, VectorType Delta) 
{
    VectorType TransformedVelocity(Delta.size());
    auto itTransformed = TransformedVelocity.begin();
    
    for(auto itProp = Delta.begin(); itProp != Delta.end(); ++itProp, ++itTransformed)
    {
        double Val = exp( - *itProp / (P0 * (1-P0)));
        double Num = V0 * ( 1 + (1./P0 - 1) * Val) * ( 1 + (1./P0 - 1) * Val);
        double Denom = (1 - P0) * (1 - P0) * Val * Val;
        *itTransformed = Num/Denom;
    }
    
    return TransformedVelocity;
}



