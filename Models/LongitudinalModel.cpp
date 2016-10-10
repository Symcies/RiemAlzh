#include <stdexcept>
#include "LongitudinalModel.h"
#include "../Manifolds/PropagationManifold.h"
#include "../Tests/TestAssert.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LongitudinalModel
::LongitudinalModel(const unsigned int NbIndependentComponents, std::shared_ptr<AbstractManifold>& M)
{

    m_OutputParameters.open("Parameters.txt", std::ofstream::out | std::ofstream::trunc);
    m_NbIndependentComponents = NbIndependentComponents;
    std::shared_ptr<PropagationManifold> Manifold = std::dynamic_pointer_cast<PropagationManifold>(M);
    m_Manifold = Manifold;
}

LongitudinalModel
::~LongitudinalModel()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::InitializeRandomVariables()
{
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    m_ManifoldRandomVariables.clear();
    InitializePopulationRandomVariables();
    InitializeIndividualRandomVariables();
    InitializeManifoldRandomVariables();
    
    m_OutputParameters << "P0, T0, V0, sigmaKsi, sigmaTau, sigma, " << std::endl;
    ComputeOutputs();

}

void 
LongitudinalModel
::UpdateParameters(std::shared_ptr<Realizations> &R, std::string Name) 
{
    /// Preprocess the name to delete the # if any
    Name = Name.substr(0, Name.find_first_of("#"));
    
    
    if(Name == "P0" or Name == "T0" or Name == "V0" or Name == "Delta")
    {   
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }
    else if(Name == "Beta")
    {
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }
    else if(Name == "S")
    {
        ComputeSpaceShifts(R);
    }
    else if(Name == "Ksi" or Name == "Tau")
    {
        // Nothing to do
    }
    else if(Name == "All")
    {
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    } 
    else
    {
        std::cout << "Should be name in LongitudinalModel > Update Parameters?" << std::endl;
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }

}

std::vector< std::vector< double >>
LongitudinalModel
::GetSufficientStatistics(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D)
{
    SufficientStatisticsVector SufficientStatistics;

    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    auto Delta = GetPropagationCoefficients(R);
    
    /////////////////////////
    /// Compute S1 and S2
    /////////////////////////
    std::vector< double > S1, S2;


    for(std::pair<Data::const_iterator, int> i(D->begin(), 0) ;
        i.first != D->end() && i.second < D->size() ;
        ++i.first, ++i.second)
    {
        /// Given a particular subject, get its attributes, then, loop over its observation
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i.second, R);
        auto SpaceShift = m_SpaceShifts.at("W" + std::to_string(i.second));
        
        for(auto it : *i.first)
        {
            double TimePoint = SubjectTimePoint(it.second);

            std::vector<double> ParallelCurve = CastedManifold->ComputeParallelCurve(InitialPosition, T0, InitialVelocity, SpaceShift, TimePoint, Delta);
            std::vector<double> Observation = it.first;

            S1.push_back( ComputeEuclideanScalarProduct(ParallelCurve, Observation ) );
            S2.push_back( ComputeEuclideanScalarProduct(ParallelCurve, ParallelCurve) );
            
        }

    }

    SufficientStatistics.push_back(S1);
    SufficientStatistics.push_back(S2);
    

    /////////////////////////
    /// Compute S3 and S4
    /////////////////////////
    std::vector< double > S3, S4;

    for(int i = 0; i < D->size() ; ++i )
    {
        double Ksi = R->at("Ksi")[i];
        double Tau = R->at("Tau")[i];

        S3.push_back(Ksi * Ksi);
        S4.push_back(Tau * Tau); 
    }

    SufficientStatistics.push_back(S3);
    SufficientStatistics.push_back(S4);

    /////////////////////////
    /// Compute S5, S6 and S7
    /////////////////////////
    std::vector< double > S5, S6, S7;
    S5.push_back( R->at("P0")[0] );
    S6.push_back( R->at("T0")[0] );
    S7.push_back( R->at("V0")[0] );

    SufficientStatistics.push_back(S5);
    SufficientStatistics.push_back(S6);
    SufficientStatistics.push_back(S7);


    /////////////////////////
    /// Compute S8
    /////////////////////////
    SufficientStatistics.push_back(Delta);


    /////////////////////////
    /// Compute S9
    /////////////////////////
    std::vector< double > S9;
    for(int i = 0; i < (m_Manifold->GetDimension() - 1) * m_NbIndependentComponents; ++i)
    {
        double Beta = R->at("Beta#" + std::to_string(i))[0];
        S9.push_back(Beta);
    }

    SufficientStatistics.push_back(S9);
    
    
    /// Tests
    /// There are two way to compute the logLikelihood. Check if they are equal
    double LogLikelihood = ComputeLogLikelihoodGeneric(R, D);
    double Calculation = 0;
    for(auto it : *D)
    {
        for(auto it2 : it)
        {
            for(auto it3 : it2.first)
            {
                Calculation += it3*it3;
            }
        }
    }
    for(auto it : SufficientStatistics[0])
    {
        Calculation += -2 * it;
    }
    for(auto it : SufficientStatistics[1])
    {
        Calculation += it;
    }
    Calculation /= -2*m_Noise->GetVariance();
    std::cout << "Likelihood difference  " << LogLikelihood - Calculation << std::endl;
    /// End tests
        
    return SufficientStatistics;
    
}


void
LongitudinalModel
::UpdateRandomVariables(const SufficientStatisticsVector& StochSufficientStatistics, const std::shared_ptr<Data>& D)
{
    
    
    
    
    /// Update P0(mean)
    auto AbstractP0 = m_PopulationRandomVariables.at("P0");
    auto P0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractP0 );
    P0->SetMean(StochSufficientStatistics[4][0]);     //TODO : Maybe store P0...

    /// Update T0(mean)
    auto AbstractT0 = m_PopulationRandomVariables.at("T0");
    auto T0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractT0 );
    T0->SetMean(StochSufficientStatistics[5][0]);    //TODO : Maybe store T0...

    /// Update V0(mean)
    auto AbstractV0 = m_PopulationRandomVariables.at("V0");
    auto V0 = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractV0 );
    V0->SetMean(StochSufficientStatistics[6][0]);    //TODO : Maybe store V0...

    /// Update Delta(k)(mean)
    for(std::pair< std::vector< double >::const_iterator, int> i(StochSufficientStatistics[7].begin(), 0) ;
        i.first != StochSufficientStatistics[7].end() && i.second < StochSufficientStatistics[7].size() ;
        ++i.first, ++i.second)
    {
        std::string Name = "Delta#" + std::to_string(i.second);
        auto AbstractDelta = m_ManifoldRandomVariables.at(Name);
        auto Delta = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractDelta ) ;
        Delta->SetMean(*i.first) ;              //TODO : Maybe store P0...
    }

    /// Update Beta(k)(mean)
    for(std::pair< std::vector< double >::const_iterator, int> i(StochSufficientStatistics[8].begin(), 0) ;
        i.first != StochSufficientStatistics[8].end() && i.second < StochSufficientStatistics[8].size() ;
        ++i.first, ++i.second )
    {
        std::string Name = "Beta#" + std::to_string(i.second);
        auto AbstractBeta = m_PopulationRandomVariables.at(Name);
        auto Beta = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractBeta );
        Beta->SetMean(*i.first);                //TODO : Maybe store P0...
    }

    /// Update Ksi
    double VarianceKsi = 0;
    for(auto it : StochSufficientStatistics[2])
    {
        VarianceKsi += it;
    }
    VarianceKsi /= StochSufficientStatistics[2].size(); // TODO : CHECK IF SIZE EQUAL TO NUMBER OF PEOPLE
    //VarianceKsi = VarianceKsi * VarianceKsi;
    auto AbstractKsi = m_IndividualRandomVariables.at("Ksi");
    auto Ksi = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractKsi );
    Ksi->SetVariance(VarianceKsi);                  //TODO : Maybe store P0...


    /// Update Tau
    double VarianceTau = 0;
    for(auto it : StochSufficientStatistics[3])
    {
        VarianceTau += it;
    }
    VarianceTau /= StochSufficientStatistics[3].size(); // TODO : CHECK IF SIZE EQUAL TO NUMBER OF PEOPLE
    //VarianceTau = VarianceTau * VarianceTau;
    
    auto AbstractTau = m_IndividualRandomVariables.at("Tau");
    auto Tau = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractTau );
    Tau->SetVariance(VarianceTau);                 //TODO : Maybe store P0...


    /// Update Uncertainty variance

    /// Initialization
    double N = m_Manifold->GetDimension();
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
    NoiseVariance /= N*K;

    m_Noise->SetVariance(NoiseVariance);
    
    

    ComputeOutputs();


}


double
LongitudinalModel
::ComputeLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D,
                    const std::pair<std::string, int> NameRandomVariable)
{
    double LogLikelihood = ComputeLogLikelihood(R, D, NameRandomVariable);
    return exp(LogLikelihood);
}


double
LongitudinalModel
::ComputeLogLikelihood(const std::shared_ptr<Realizations> &R, const std::shared_ptr<Data> &D,
                       const std::pair<std::string, int> NameRandomVariable) 
{
    /// Get the name of the realization, and its number (in case it is subject specific)
    std::string Name = NameRandomVariable.first.substr(0, NameRandomVariable.first.find_first_of("#"));
    int SubjectNumber = NameRandomVariable.second;


    bool PreviousEqualCurrentRealizations = (*R == std::get<2>(m_LastLogLikelihood));
    bool CurrentIsGeneric = !(Name == "Tau" or Name == "Ksi");
    bool PreviousIsGeneric = std::get<0>(m_LastLogLikelihood);

    /// COMPUTE LIKELIHOOD
    double LogLikelihood;
    if (PreviousEqualCurrentRealizations && CurrentIsGeneric && PreviousIsGeneric)
    {
        LogLikelihood = std::get<1>(m_LastLogLikelihood);
    }
    else if (!CurrentIsGeneric)
    {
        LogLikelihood = ComputeLogLikelihoodIndividual(R, D, SubjectNumber);
    }
    else
    {
        /// Here it can be assumed that CurrentIsGeneric / (CurrentIs Generic && PreviousIsGeneric && PreviousEqualCurrentRealizations)
        if(Name == "S")
        {
            ComputeSpaceShifts(R);
        }
        else if(Name == "Beta")
        {
            ComputeAMatrix(R);
            ComputeSpaceShifts(R);
        }
        else
        {
            // The else condition contains the P0, T0, V0 and Delta random variables
            ComputeOrthonormalBasis(R);
            ComputeAMatrix(R);
            ComputeSpaceShifts(R);
        }
        LogLikelihood = ComputeLogLikelihoodGeneric(R, D);
    }

    m_LastLogLikelihood = std::tuple<bool, double, Realizations>(CurrentIsGeneric, LogLikelihood, *R);
    return LogLikelihood;
}

std::vector< std::vector< std::pair< std::vector<double>, double> > >
LongitudinalModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs)
{



    //////////////////////////
    // Simulate Time Point ///
    //////////////////////////
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::normal_distribution<double> Normal(70.0, sqrt(3.0));

    std::vector< std::vector< double >> TimePoint;
    for(int i = 0; i<NumberOfSubjects; ++i)
    {
        std::vector<double> SubjectTimePoint;
        for(int j = 0; j < Uni(RNG) ; ++j)
        {
            SubjectTimePoint.push_back( Normal(RNG));
        }
        std::sort(SubjectTimePoint.begin(), SubjectTimePoint.end());
        TimePoint.push_back(SubjectTimePoint);
    }

    ///////////////////
    // Realizations ///
    ///////////////////
    auto R = std::make_shared<Realizations>( SimulateRealizations(NumberOfSubjects) );


    /////////////////////////////
    // Compute some functions ///
    /////////////////////////////
    ComputeOrthonormalBasis(R);
    ComputeAMatrix(R);
    ComputeSpaceShifts(R);


    /////////////////////////
    /// Get the data
    /////////////////////////
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    std::vector<double> InitialPosition = GetInitialPosition(R);
    std::vector<double> InitialVelocity = GetInitialVelocity(R);
    std::vector<double> Delta = GetPropagationCoefficients(R);


    //////////////////////////////////////////////////////////////////////////
    // Simulate Data
    //////////////////////////////////////////////////////////////////////////
    Data D;



    std::normal_distribution<double> Normal2(0.0, sqrt(m_Noise->GetVariance()));
    int i = 0;

    for(auto it : TimePoint)
    {
        /// Given a particular subject, get its attributes, then, loop over its observation
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        auto SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));

        IndividualData ID;
        for(auto it2 : it)
        {
            std::vector<double> Observation;

            double Time = SubjectTimePoint(it2);
            auto ParallelCurve = CastedManifold->ComputeParallelCurve(InitialPosition, T0, InitialVelocity, SpaceShift, Time, Delta);
            for(auto it3 : ParallelCurve)
            {
                double Gen = Normal2(RNG);
                Observation.push_back(it3 + Gen);
            }
            ID.push_back(std::pair<std::vector<double>, double> (Observation, it2));

        }
        D.push_back(ID);
        i+= 1;
    }

    return D;

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::InitializeFakeRandomVariables()
{
    /////////////////////////////
    /// Population Parameters ///
    /////////////////////////////

    // Initial Position
    double P0Mean = 0.3;
    double P0Variance = 0.01;
    auto P0 = std::make_shared<GaussianRandomVariable>(P0Mean, P0Variance);
    RandomVariable P0_("P0", P0);
    m_PopulationRandomVariables.insert(P0_);

    // Initial Time
    double T0Mean = 70;
    double T0variance = 0.1;
    auto T0 = std::make_shared<GaussianRandomVariable>(T0Mean, T0variance);
    RandomVariable T0_("T0", T0);
    m_PopulationRandomVariables.insert(T0_);

    // Initial Velocity
    double V0Mean = 0.01;
    double V0Variance = 0.001;
    auto V0 = std::make_shared<GaussianRandomVariable>(V0Mean, V0Variance);
    RandomVariable V0_("V0", V0);
    m_PopulationRandomVariables.insert(V0_);

    // Initial Beta coefficient
    double BetaVariance = 0.01;
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        double BetaMean = (double)i;
        auto Beta = std::make_shared< GaussianRandomVariable> (BetaMean, BetaVariance);
        std::string name = "Beta#" + std::to_string(i);
        RandomVariable Beta_(name, Beta);
        m_PopulationRandomVariables.insert(Beta_);
    }

    // Noise
    double NoiseMean = 0.0;
    double NoiseVariance = 0.05;
    m_Noise = std::make_shared< GaussianRandomVariable >(NoiseMean, NoiseVariance);

    /////////////////////////////
    ///  Manifold Parameters  ///
    /////////////////////////////

    std::map<std::string, double> ManifoldParameters;

    double DeltaVariance = 0.001;
    for(int i = 0; i < m_Manifold->GetDimension() - 1 ; ++i)
    {
        double DeltaMean = 0.5 + (double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        std::string name = "Delta#" + std::to_string(i);
        RandomVariable Delta_(name, Delta);
        ManifoldParameters.insert(std::pair<std::string, double> (name, Delta->Sample()));
        m_ManifoldRandomVariables.insert(Delta_);
    }


    /////////////////////////////
    /// Individual Parameters ///
    /////////////////////////////

    // Initial pre acceleration factor
    double KsiMean = 0.0;
    double KsiVariance = 0.5;
    auto Ksi = std::make_shared< GaussianRandomVariable >(KsiMean, KsiVariance);
    RandomVariable Ksi_("Ksi", Ksi);
    m_IndividualRandomVariables.insert(Ksi_);

    // Initial Time Shift
    double TauMean = 0.0;
    double TauVariance = 3.0;
    auto Tau = std::make_shared< GaussianRandomVariable >(TauMean, TauVariance);
    RandomVariable Tau_("Tau", Tau);
    m_IndividualRandomVariables.insert(Tau_);

    // Initial Space shift coefficient
    double SLocation = 0.0;
    double SScale = 1.0/2.0;
    auto S = std::make_shared< LaplaceRandomVariable >(SLocation, SScale);
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        RandomVariable S_("S#" + std::to_string(i), S);
        m_IndividualRandomVariables.insert(S_);
    }

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::InitializePopulationRandomVariables()
{
    // Initial Position
    double P0Mean = 0.35;
    double P0Variance = 0.01;
    auto P0 = std::make_shared<GaussianRandomVariable>(P0Mean, P0Variance);
    RandomVariable P0_("P0", P0);
    m_PopulationRandomVariables.insert(P0_);

    // Initial Time
    double T0Mean = 70;
    double T0variance = 0.1;
    auto T0 = std::make_shared<GaussianRandomVariable>(T0Mean, T0variance);
    RandomVariable T0_("T0", T0);
    m_PopulationRandomVariables.insert(T0_);

    // Initial Velocity
    double V0Mean = 0.06;
    double V0Variance = 0.001;
    auto V0 = std::make_shared<GaussianRandomVariable>(V0Mean, V0Variance);
    RandomVariable V0_("V0", V0);
    m_PopulationRandomVariables.insert(V0_);

    // Initial Beta coefficient
    double BetaVariance = 0.01;
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        double BetaMean = (double)i + 0.2;
        auto Beta = std::make_shared< GaussianRandomVariable> (BetaMean, BetaVariance);
        std::string name = "Beta#" + std::to_string(i);
        RandomVariable Beta_(name, Beta);
        m_PopulationRandomVariables.insert(Beta_);
    }


    // Noise
    double NoiseMean = 0.2;
    double NoiseVariance = 0.02;
    m_Noise = std::make_shared< GaussianRandomVariable >(NoiseMean, NoiseVariance);
    
}

void
LongitudinalModel
::InitializeIndividualRandomVariables()
{
    // Initial pre acceleration factor
    double KsiMean = 0.0;
    double KsiVariance = 0.5;
    auto Ksi = std::make_shared< GaussianRandomVariable >(KsiMean, KsiVariance);
    RandomVariable Ksi_("Ksi", Ksi);
    m_IndividualRandomVariables.insert(Ksi_);

    // Initial Time Shift
    double TauMean = 0.0;
    double TauVariance = 3.0;
    auto Tau = std::make_shared< GaussianRandomVariable >(TauMean, TauVariance);
    RandomVariable Tau_("Tau", Tau);
    m_IndividualRandomVariables.insert(Tau_);

    // Initial Space shift coefficient
    double SLocation = 0.0;
    double SScale = 1.0/2.0;
    auto S = std::make_shared< LaplaceRandomVariable >(SLocation, SScale);
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        RandomVariable S_("S#" + std::to_string(i), S);
        m_IndividualRandomVariables.insert(S_);
    }

}

void
LongitudinalModel
::InitializeManifoldRandomVariables()
{

    // Initial Propagation coefficient
    double DeltaVariance = 0.001;
    for(int i = 0; i < m_Manifold->GetDimension() - 1; ++i)
    {
        double DeltaMean = 0.5 + (double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        std::string name = "Delta#" + std::to_string(i);
        RandomVariable Delta_(name, Delta);
        m_ManifoldRandomVariables.insert(Delta_);
    }
}

std::vector<double>
LongitudinalModel
::GetInitialPosition(const std::shared_ptr<Realizations> &R)
{
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);

    double T0 = R->at("T0")[0];
    std::vector<double> P0, V0, Delta;
    P0.push_back( R->at("P0")[0] );
    V0.push_back( R->at("V0")[0] );
    for(int i = 0; i < m_Manifold->GetDimension() - 1; ++i)
    {
        Delta.push_back( R->at("Delta#" + std::to_string(i))[0]);
    }

    auto InitialPosition = CastedManifold->ComputeGeodesic(P0, T0, V0, T0, Delta);
    
    /// Tests
    std::function<double()> f1 = [&R]() { return R->at("P0")[0]; };
    std::function<double()> f2 = [&InitialPosition]() { return InitialPosition[0]; };
    TestAssert::WarningEquality_Function(f1, f2, "P0 != InitialPosition[0]. LongitudinalModel > GetInitialPosition");
    
    return InitialPosition;
}

std::vector<double>
LongitudinalModel
::GetInitialVelocity(const std::shared_ptr<Realizations> &R)
{
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);

    double T0 = R->at("T0")[0];
    std::vector<double> P0, V0, Delta;
    P0.push_back( R->at("P0")[0] );
    V0.push_back( R->at("V0")[0] );
    for(int i = 0; i < m_Manifold->GetDimension() - 1; ++i)
    {
        Delta.push_back( R->at("Delta#" + std::to_string(i))[0]);
    }

    auto InitialVelocity = CastedManifold->ComputeGeodesicDerivative(P0, T0, V0, T0, Delta);
    
    /// Tests 
    std::function<double()> f1 = [&R]() { return R->at("V0")[0]; };
    std::function<double()> f2 = [&InitialVelocity]() { return InitialVelocity[0]; };
    TestAssert::WarningEquality_Function(f1, f2, "P0 != InitialVelocity[0]. LongitudinalModel > GetInitialVelocity");
    
    return InitialVelocity;
    
}

std::vector<double>
LongitudinalModel
::GetPropagationCoefficients(const std::shared_ptr<Realizations> &R)
{
    std::vector<double> Delta;

    for(int i = 0; i < m_Manifold->GetDimension() - 1; ++i)
    {
        Delta.push_back( R->at("Delta#" + std::to_string(i))[0]);
    }

    return Delta;
}

std::function<double(double)>  
LongitudinalModel
::GetSubjectTimePoint(const int SubjectNumber, const std::shared_ptr<Realizations>& R)
{
    double AccFactor = exp(R->at("Ksi")[SubjectNumber]);
    double TimeShift = R->at("Tau")[SubjectNumber];
    double T0 = R->at("T0")[0];
    
    return [AccFactor, TimeShift, T0](double t) { return AccFactor * (t - TimeShift - T0) + T0; };
}

void
LongitudinalModel
::ComputeOrthonormalBasis(const std::shared_ptr<Realizations>& R)
{
    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    std::vector<double> InitialPosition = GetInitialPosition(R);
    std::vector<double> InitialVelocity = GetInitialVelocity(R);
    std::vector<double> Delta = GetPropagationCoefficients(R);

    /// Compute the transformation to do the Householder reflection in a Euclidean space
    std::vector<double> U = CastedManifold->GetVelocityTransformToEuclideanSpace(InitialPosition, T0, InitialVelocity, Delta);

    /// Compute the initial pivot vector U
    double Norm = 0;
    for(auto it : U)
    {
        Norm += it * it;
    }
    Norm = sqrt(Norm);

    U[0] += copysign(1, -U[0])*Norm;

    double NormU = 0;
    for(auto it : U)
    {
        NormU += it * it;
    }

    /// Compute Q = I(N) - 2U . Ut / (Ut . U)
    std::vector< std::vector< double > > Q;
    for(auto it : U)
    {
        std::vector<double> Coordinate;
        for(auto it2 : U)
        {
            Coordinate.push_back( - 2 * it * it2 / NormU);
        }
        Q.push_back( Coordinate );
    }

    for(int i = 0; i < Q.size() ; ++i)
    {
        Q[i][i] += 1;
    }
    
    /// TESTS
    
    /// Test colinearity between the first vector and Q[0]
    

    /// Orthogonality between the basis vectors and the velocity
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
    Q.erase(Q.begin());
    m_OrthogonalBasis = Q;
}

void
LongitudinalModel
::ComputeAMatrix(const std::shared_ptr<Realizations>& R)
{
    std::vector<std::vector<double>> AMatrix;
    for(int i = 0; i < m_NbIndependentComponents ; ++i)
    {
        std::vector<double> Beta;
        for(int j = 0; j < m_Manifold->GetDimension() - 1 ; ++j)
        {
            std::string Number = std::to_string(int(j + i*(m_Manifold->GetDimension() - 1)));
            Beta.push_back( R->at( "Beta#" + Number)[0]);

        }
        std::vector<double> AColumn = LinearCombination(Beta, m_OrthogonalBasis);
        AMatrix.push_back( AColumn );
    }
    
    m_AMatrix = AMatrix;
    
    
    /// TESTS
    std::vector<double> InitialPosition = GetInitialPosition(R);
    std::vector<double> InitialVelocity = GetInitialVelocity(R);
    
    /// Check if the basis is orthogonal to the velocity
    for(auto it : m_AMatrix)
    {
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
}

void
LongitudinalModel
::ComputeSpaceShifts(const std::shared_ptr<Realizations>& R)
{
    std::map< std::string, std::vector<double>> SpaceShifts;
    int NumberOfSubjects = (int)R->at("Tau").size();
    if(NumberOfSubjects != R->at("Ksi").size())
    {
        throw std::invalid_argument("Not the same number of realization in Tau and in Ksi. Whereas its equal to the number of subjects");
    }

    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        std::vector<double> Si;
        for(int j = 0; j < m_NbIndependentComponents; ++j)
        {
            double Realization = R->at("S#" + std::to_string(j))[i];
            Si.push_back(Realization);
        }


        std::pair< std::string, std::vector<double> > SpaceShift("W"+std::to_string(i), LinearCombination(Si, m_AMatrix) );
        SpaceShifts.insert( SpaceShift);
    }

    m_SpaceShifts = SpaceShifts;


    ///////////////////////////////
    //// Debugging : Unit tests ///
    ///////////////////////////////
    
    
    /// Get Data
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    std::vector<double> InitialPosition = GetInitialPosition(R);
    std::vector<double> InitialVelocity = GetInitialVelocity(R);
    std::vector<double> Delta = GetPropagationCoefficients(R);
    std::function<double()> NullFunction = []() { return 0;};
    
    /// DEBUG 1 : The orthonormal basis is orthogonal to the velocity
    for(auto it : m_AMatrix)
    {
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


    /// DEBUG 4 : ParallelTransport(0) = W(0). Is it true?
    /*
    for(auto it : m_SpaceShifts)
    {
        auto ParallelTransport = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, it.second, T0, Delta);
        TestAssert::WarningEquality_Object(it.second, ParallelTransport, "Parallel Transport(T0) != SpaceShift. Longitudinal Model > ComputeSpaceShift");
    }
     */


    /// DEBUG 5 : <diff(g(t)) , ParallelTransport>|g(t) = 0
    /// DEBUG 6 : <ParallelTransport, ParallelTransport>g(t) = Constant
    auto ParallelTransport = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, m_SpaceShifts["W0"], 55, Delta);
    auto position = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, 55, Delta);
    double Constant = m_Manifold->ComputeScalarProduct(ParallelTransport, ParallelTransport, position);
    std::function<double()> ConstantFunction = [Constant] () { return Constant; };
    
    
    for(double t = 61; t < 80; t += 0.5)
    {
        
        std::function<double()> f1 = [&CastedManifold, &InitialPosition, &InitialVelocity, &T0, &Delta, &SpaceShifts, t] () 
        {
            auto a = CastedManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, t, Delta);
            auto b = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, SpaceShifts["W0"], t, Delta);
            auto velocity = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
            return CastedManifold->ComputeScalarProduct(a, b, velocity);
        };
    
        std::function<double()> f2= [&SpaceShifts, &CastedManifold, &InitialPosition, &InitialVelocity, &T0, &Delta, t] () 
        {
            auto a = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, SpaceShifts["W0"], t, Delta);
            auto position = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
            return CastedManifold->ComputeScalarProduct(a, a, position);
        };
    
        
        
        
        TestAssert::WarningEquality_Function(f1, NullFunction, "Parallel transport not orthogonal to velocity at any t. LongitudinalModel > ComputeSpaceShifts");
        TestAssert::WarningEquality_Function(f2, ConstantFunction, "Norm of parallel transport not constant over time. LongiudinalModel > ComputeSpaceShifts");
    }
    




    /// DEBUG 6 : <ParallelTransport, ParallelTransport>g(t) = Constant
    /*
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

double
LongitudinalModel
::ComputeLogLikelihoodGeneric(const std::shared_ptr<Realizations> &R, const std::shared_ptr<Data> &D)
{

    /// Get the data
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    std::vector<double> InitialPosition = GetInitialPosition(R);
    std::vector<double> InitialVelocity = GetInitialVelocity(R);
    std::vector<double> Delta = GetPropagationCoefficients(R);
    
    /// Test the data 
    std::function<double()> f1 = [&InitialPosition]() { return InitialPosition[0];};
    std::function<double()> f2 = [&R]() { return R->at("P0")[0];};
    TestAssert::WarningEquality_Function(f1, f2, "P0 != gamma(t0) in Longitudinal model -> ComputeLogLikelihoodGeneric");

    /// Compute the likelihood
    double LogLikelihood = 0;

    for(std::pair<Data::const_iterator, int> i(D->begin(), 0);
        i.first != D->end() && i.second < D->size();
        ++i.first, ++i.second)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i.second, R);
        std::vector<double> SpaceShift = m_SpaceShifts.at("W" + std::to_string(i.second));

        for(auto it : *i.first)
        {

            
            double TimePoint = SubjectTimePoint(it.second);

            std::vector<double> ParallelCurve = CastedManifold->ComputeParallelCurve(InitialPosition, T0, InitialVelocity, SpaceShift, TimePoint, Delta);
            std::vector<double> Observation = it.first;
            
            LogLikelihood  += NormOfVectorDifference(Observation, ParallelCurve);
        }
    }
    
    LogLikelihood  /= -2*m_Noise->GetVariance();

    return LogLikelihood ;
}


double
LongitudinalModel
::ComputeLogLikelihoodIndividual(const std::shared_ptr<Realizations> &R, const std::shared_ptr<Data>& D,
                              const int SubjectNumber)
{
    /// Initialize the individual parameters
    double AccFactor = exp(R->at("Ksi")[SubjectNumber]);
    double TimeShift = R->at("Tau")[SubjectNumber];
    std::vector<double> SpaceShift = m_SpaceShifts.at("W" + std::to_string(SubjectNumber));

    /// Get the global parameters
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    std::vector<double> InitialPosition = GetInitialPosition(R);
    std::vector<double> InitialVelocity = GetInitialVelocity(R);
    std::vector<double> Delta = GetPropagationCoefficients(R);

    /// Compute the likelihood
    double LogLikelihood = 0;
    for(auto it = D->at(SubjectNumber).begin(); it != D->at(SubjectNumber).end(); ++it)
    {
        std::vector<double> Observation = it->first;
        double Tij = it->second;
        double TimePoint = AccFactor * (Tij - T0 - TimeShift) + T0;

        std::vector<double> ParallelCurve = CastedManifold->ComputeParallelCurve(InitialPosition, T0, InitialVelocity, SpaceShift, TimePoint, Delta);

        LogLikelihood += NormOfVectorDifference(Observation, ParallelCurve);
    }

    LogLikelihood  /= -2*m_Noise->GetVariance();
    return LogLikelihood;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::ComputeOutputs()
{
    double MeanP0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"))->GetMean();
    double MeanT0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("T0"))->GetMean();
    double MeanV0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("V0"))->GetMean();
    double SigmaKsi = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetVariance();
    double SigmaTau = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetVariance();
    double Sigma = m_Noise->GetVariance();

    m_OutputParameters << MeanP0 << ", " << MeanT0 << ", " << MeanV0 << ", ";
    m_OutputParameters << SigmaKsi << ", " << SigmaTau << ", " << Sigma << ", " << std::endl;

    
    std::cout << "Parameters : P0: " << MeanP0 << ". T0: " << MeanT0 << ". V0: " << MeanV0 << ". Ksi: " << SigmaKsi << ". Tau: " << SigmaTau << ". Sigma: " << Sigma << ". " << std::endl;
}