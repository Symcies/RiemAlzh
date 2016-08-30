#include "LongitudinalModel.h"



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LongitudinalModel
::LongitudinalModel(unsigned int NbIndependentComponents)
{
    m_NbIndependentComponents = NbIndependentComponents;
}

LongitudinalModel
::~LongitudinalModel()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::SetRealization(std::string Name, double Realization)
{
    // TODO : CHECK IF THIS FUNCTION IS NEEDED AND HOW IT INTERACTS WITH INITIALITE_MODEL_PARAMETERS!
    // TODO : Especially when one knows if the realizations are / aren't in the model
    if(m_IndividualRandomVariables.count(Name))
    {

    }
    else if(m_PopulationRandomVariables.count(Name))
    {

    }
    else if(m_ManifoldRandomVariables.count(Name))
    {
        m_Manifold->SetParameter(Name, Realization);
    }
    else
    {
        std::cout << "OWN WARNING : The key does not exist in the map" << std::endl;
    }
}

std::shared_ptr< AbstractRandomVariable >
LongitudinalModel
::GetCandidateRandomVariable(const std::string Name, const double Realization)
{
    /// Initialization
    std::map<std::string, double> CandidateVariance;
    CandidateVariance.insert(std::pair<std::string, double> ("P0", 0.3));
    CandidateVariance.insert(std::pair<std::string, double> ("V0", 0.01));
    CandidateVariance.insert(std::pair<std::string, double> ("T0", 0.01));

    for(int i = 0; i < m_Manifold->GetDimension() - 1; ++i)
    {
        CandidateVariance.insert(std::pair<std::string, double> ("Delta" + std::to_string(i), 0.01));
    }

    for(int i = 0; i < (m_Manifold->GetDimension() - 1) * m_NbIndependentComponents; ++i)
    {
        CandidateVariance.insert(std::pair<std::string, double> ("Beta" + std::to_string(i), 0.01));
    }

    for(int i = 0; i < m_NbIndependentComponents; ++i)
    {
        CandidateVariance.insert(std::pair<std::string, double> ("S" + std::to_string(i), 0.01));
    }

    CandidateVariance.insert(std::pair<std::string, double> ("Ksi", 0.01));
    CandidateVariance.insert(std::pair<std::string, double> ("Tau", 0.01));

    /// Return part
    std::shared_ptr<AbstractRandomVariable> CandidateRandomVariable = std::make_shared<GaussianRandomVariable>( Realization, CandidateVariance.at(Name));

    return CandidateRandomVariable;
}

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


}

void
LongitudinalModel
::InitializeModelParameters(const std::shared_ptr<Realizations> &R)
{

    // Initialize the manifold parameters, if any
    std::map<std::string, double> Parameters;
    for(RandomVariableMap::iterator it = m_ManifoldRandomVariables.begin(); it != m_ManifoldRandomVariables.end(); ++it)
    {
        Parameters.insert(std::pair<std::string, double> (it->first, R->at(it->first)[0]));
    }
    m_Manifold->InitializeParameters(Parameters);

    // TODO : Here is the first place in the model where we have access to the realizations
    // TODO : Thus, if the candidate random variables depends on the realization, then it might be the right place
    // TODO : CHECK IF THIS IS THE RIGHT PLACE TO INITIALIZE THE PROPOSITION LAW
    //InitializeCandidateRandomVariables(R);
}

std::vector< std::vector< double >>
LongitudinalModel
::GetSufficientStatistics(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D)
{
    SufficientStatisticsVector SufficientStatistics;

    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    std::vector<double> P0, V0;
    for(int i = 0; i < m_Manifold->GetDimension(); ++i)
    {
        P0.push_back(R->at("P0")[0]);
        V0.push_back(R->at("V0")[0]);
    }

    /////////////////////////
    /// Compute S1 and S2
    /////////////////////////
    std::vector< double > S1, S2;

    /// Get population wide attribute , then, Loop over all the subjects
    double T0 = R->at("T0")[0];

    for(std::pair<Data::const_iterator, int> i(D->begin(), 0) ;
        i.first != D->end() && i.second < D->size() ;
        ++i.first, ++i.second)
    {
        /// Given a particular subject, get its attributes, then, loop over its observation
        double Ksi = R->at("Ksi")[i.second];
        double Tau = R->at("Tau")[i.second];
        std::vector< double > SpaceShift = m_SpaceShifts.at("W" + std::to_string(i.second));

        for(IndividualData::const_iterator it = i.first->begin() ; it != i.first->end() ; ++it)
        {
            double TimePoint = exp(Ksi) * ( it->second - T0 - Tau) + T0;

            std::vector<double> ParallelCurve = m_Manifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint);
            std::vector<double> Observation = it->first;

            S1.push_back( ComputeEuclideanScalarProduct(ParallelCurve, Observation ) );
            S2.push_back( ComputeEuclideanScalarProduct(ParallelCurve, ParallelCurve) );
        }

    }

    SufficientStatistics.push_back(S1);
    SufficientStatistics.push_back(S2);

    //std::cout << "1. Calculation : " << SufficientStatistics[0][0] << " & " << SufficientStatistics[0][1] << std::endl;

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
    std::vector< double > S8;
    for(int i = 0; i<m_Manifold->GetDimension() - 1; ++i)
    {
        double Delta = R->at("Delta" + std::to_string(i))[0];
        S8.push_back(Delta);
    }
    
    SufficientStatistics.push_back(S8);


    /////////////////////////
    /// Compute S9
    /////////////////////////
    std::vector< double > S9;
    for(int i = 0; i < (m_Manifold->GetDimension() - 1) * m_NbIndependentComponents; ++i)
    {
        double Beta = R->at("Beta" + std::to_string(i))[0];
        S9.push_back(Beta);
    }

    SufficientStatistics.push_back(S9);


    return SufficientStatistics;
    
}


void
LongitudinalModel
::UpdateRandomVariables(const SufficientStatisticsVector& SufficientStatistics, const std::shared_ptr<Data>& D)
{

    /// Update P0(mean)
    auto AbstractP0 = m_PopulationRandomVariables.at("P0");
    auto P0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractP0 );
    P0->SetMean(SufficientStatistics[4][0]);     //TODO : Maybe store P0...

    /// Update T0(mean)
    auto AbstractT0 = m_PopulationRandomVariables.at("T0");
    auto T0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractT0 );
    T0->SetMean(SufficientStatistics[5][0]);    //TODO : Maybe store T0...

    /// Update V0(mean)
    auto AbstractV0 = m_PopulationRandomVariables.at("V0");
    auto V0 = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractV0 );
    V0->SetMean(SufficientStatistics[6][0]);    //TODO : Maybe store V0...

    /// Update Delta(k)(mean)
    for(std::pair< std::vector< double >::const_iterator, int> i(SufficientStatistics[7].begin(), 0) ;
        i.first != SufficientStatistics[7].end() && i.second < SufficientStatistics[7].size() ;
        ++i.first, ++i.second)
    {
        std::string Name = "Delta" + std::to_string(i.second);
        auto AbstractDelta = m_ManifoldRandomVariables.at(Name);
        auto Delta = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractDelta ) ;
        Delta->SetMean(*i.first) ;              //TODO : Maybe store P0...
    }

    /// Update Beta(k)(mean)
    for(std::pair< std::vector< double >::const_iterator, int> i(SufficientStatistics[8].begin(), 0) ;
        i.first != SufficientStatistics[8].end() && i.second < SufficientStatistics[8].size() ;
        ++i.first, ++i.second )
    {
        std::string Name = "Beta" + std::to_string(i.second);
        auto AbstractBeta = m_PopulationRandomVariables.at(Name);
        auto Beta = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractBeta );
        Beta->SetMean(*i.first);                //TODO : Maybe store P0...
    }

    /// Update Ksi
    double VarianceKsi = 0;
    for(std::vector< double >::const_iterator it = SufficientStatistics[2].begin() ; it != SufficientStatistics[2].end() ; ++it)
    {
        VarianceKsi += *it;
    }
    VarianceKsi /= SufficientStatistics[2].size(); // TODO : CHECK IF SIZE EQUAL TO NUMBER OF PEOPLE
    auto AbstractKsi = m_IndividualRandomVariables.at("Ksi");
    auto Ksi = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractKsi );
    Ksi->SetVariance(VarianceKsi);                  //TODO : Maybe store P0...


    /// Update Tau
    double VarianceTau = 0;
    for(std::vector<double>::const_iterator it = SufficientStatistics[3].begin() ; it != SufficientStatistics[3].end() ; ++it)
    {
        VarianceTau += *it;
    }
    VarianceTau /= SufficientStatistics[3].size(); // TODO : CHECK IF SIZE EQUAL TO NUMBER OF PEOPLE

    auto AbstractTau = m_IndividualRandomVariables.at("Tau");
    auto Tau = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractTau );
    Tau->SetVariance(VarianceTau);                 //TODO : Maybe store P0...


    /// Update Uncertainty variance

    /// Initialization
    double N = m_Manifold->GetDimension();
    double K = 0;

    /// Sum Yijk²
    double NoiseVariance;
    for(Data::const_iterator it = D->begin() ; it != D->end() ; ++it)
    {
        K += it->size();
        for(IndividualData::const_iterator it2 = it->begin() ; it2 != it->end() ; ++it2)
        {
            for(std::vector<double>::const_iterator it3 = it2->first.begin() ; it3 != it2->first.end() ; ++it3)
            {
                NoiseVariance += *it3 * *it3;
            }
        }
    }

    //std::cout << "1. Noise Var : " << NoiseVariance << std::endl;

    /// Sum -2 S1
    for(std::vector< double >::const_iterator it = SufficientStatistics[0].begin() ; it != SufficientStatistics[0].end() ; ++it)
    {
        NoiseVariance -= 2* *it;
        //std::cout << *it << ". ";
    }
    //std::cout << std::endl;

    //std::cout << "2. Noise Var : " << NoiseVariance << std::endl;

    /// Sum S2
    for(std::vector< double >::const_iterator it = SufficientStatistics[1].begin() ; it != SufficientStatistics[1].end(); ++it)
    {
        NoiseVariance += *it;
    }


    // Divide by N*K, then take the square root
    NoiseVariance /= N*K;
    NoiseVariance = sqrt(NoiseVariance);

    m_Noise->SetVariance(NoiseVariance);


}


double
LongitudinalModel
::ComputeLikelihood(const std::shared_ptr<Realizations>& R, const std::shared_ptr<Data>& D)
{
    // TODO : Choose where the following function should be used
    ComputeOrthonormalBasis(R);
    ComputeAMatrix(R);
    ComputeSpaceShifts(R, D->size());

    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    std::vector<double> P0, V0;
    for(int i = 0; i < m_Manifold->GetDimension(); ++i)
    {
        P0.push_back(R->at("P0")[0]);
        V0.push_back(R->at("V0")[0]);
    }

    double Likelihood = 0;

    double T0 = R->at("T0")[0];

    for(std::pair<Data::const_iterator, int> i(D->begin(), 0);
        i.first != D->end() && i.second < D->size();
        ++i.first, ++i.second)
    {
        double AccFactor = exp( R->at("Ksi")[i.second]);
        double TimeShift = R->at("Tau")[i.second];
        std::vector<double> SpaceShift = m_SpaceShifts.at("W" + std::to_string(i.second));

        //std::cout << "1. Likelihood -> ParallelCurve : " << std::endl;
        for(IndividualData::const_iterator it2 = i.first->begin(); it2 != i.first->end(); ++it2)
        {
            std::vector<double> Observation = it2->first;

            double Tij = it2->second;
            double TimePoint = AccFactor * (Tij - T0 - TimeShift) + T0;

            std::vector<double> ParallelCurve = m_Manifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint);
            /*
            if(isnan(ParallelCurve[0]))
            {
                std::cout << "Parallel Curve 0 : " << ParallelCurve[0] << std::endl;
                std::cout << "Args : " << P0[0] << " & " << V0[0] << " & " << SpaceShift[0] << std::endl;
            }
            */


            Likelihood += NormOfVectorDifference(Observation, ParallelCurve);
        }
    }
    //std::cout << "2. Likelihood : " << Likelihood << std::endl;


    //std::cout << "1. Likelihood & Variance : " << Likelihood << " & " << m_Noise->GetVariance() << std::endl;
    Likelihood /= -2*m_Noise->GetVariance();
    //std::cout << "2. Likelihood : " << Likelihood << std::endl;
    Likelihood = exp(Likelihood);
    //std::cout << "3. Likelihood : " << Likelihood << std::endl;

    return Likelihood;
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
    std::normal_distribution<double> Normal(70.0, 3.0);

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
    Realizations Real = SimulateRealizations(NumberOfSubjects);
    std::shared_ptr<Realizations> R = std::make_shared<Realizations>(Real);


    /////////////////////////////
    // Compute some functions ///
    /////////////////////////////
    // TODO : Choose where the following function should be used
    ComputeOrthonormalBasis(R);
    ComputeAMatrix(R);
    ComputeSpaceShifts(R, NumberOfSubjects);


    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    std::vector<double> P0, V0;
    for(int i = 0; i < m_Manifold->GetDimension(); ++i)
    {
        P0.push_back(R->at("P0")[0]);
        V0.push_back(R->at("V0")[0]);
    }

    //////////////////////////////////////////////////////////////////////////
    // Simulate Data
    //////////////////////////////////////////////////////////////////////////
    Data D;

    double T0 = R->at("T0")[0];
    std::normal_distribution<double> Normal2(0.0, m_Noise->Sample());
    int i = 0;

    for(std::vector<std::vector<double>>::iterator it = TimePoint.begin() ; it != TimePoint.end() ; ++it)
    {
        double Tau = R->at("Tau")[i];
        IndividualData ID;
        for(std::vector<double>::iterator it2 = it->begin() ; it2 != it->end(); ++it2)
        {
            std::vector<double> Observation;

            double Time = exp(R->at("Ksi")[i]) * (*it2 - T0 - Tau) + T0;
            std::vector<double> SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
            std::vector<double> ParallelCurve = m_Manifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, Time);
            for(std::vector<double>::iterator it3 = ParallelCurve.begin() ; it3 != ParallelCurve.end(); ++it3)
            {
                double Gen = Normal2(RNG);
                Observation.push_back(*it3 + Gen);
            }
            ID.push_back(std::pair<std::vector<double>, double> (Observation, *it2));

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
    double P0Mean = 0.5;
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
    double V0Variance = 0.01;
    auto V0 = std::make_shared<GaussianRandomVariable>(V0Mean, V0Variance);
    RandomVariable V0_("V0", V0);
    m_PopulationRandomVariables.insert(V0_);

    // Initial Beta coefficient
    double BetaVariance = 0.1;
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        double BetaMean = (double)i;
        auto Beta = std::make_shared< GaussianRandomVariable> (BetaMean, BetaVariance);
        std::string name = "Beta" + std::to_string(i);
        RandomVariable Beta_(name, Beta);
        m_PopulationRandomVariables.insert(Beta_);
    }

    // Noise
    double NoiseMean = 0.0;
    double NoiseVariance = 0.01;
    m_Noise = std::make_shared< GaussianRandomVariable >(NoiseMean, NoiseVariance);

    /////////////////////////////
    ///  Manifold Parameters  ///
    /////////////////////////////

    std::map<std::string, double> ManifoldParameters;

    double DeltaVariance = 0.1;
    for(int i = 0; i < m_Manifold->GetDimension() - 1 ; ++i)
    {
        double DeltaMean = 0.5*(double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        std::string name = "Delta" + std::to_string(i);
        RandomVariable Delta_(name, Delta);
        ManifoldParameters.insert(std::pair<std::string, double> (name, Delta->Sample()));
        m_ManifoldRandomVariables.insert(Delta_);
    }
    m_Manifold->InitializeParameters(ManifoldParameters);


    /////////////////////////////
    /// Individual Parameters ///
    /////////////////////////////

    // Initial pre acceleration factor
    double KsiMean = 0.0;
    double KsiVariance = 0.1;
    auto Ksi = std::make_shared< GaussianRandomVariable >(KsiMean, KsiVariance);
    RandomVariable Ksi_("Ksi", Ksi);
    m_IndividualRandomVariables.insert(Ksi_);

    // Initial Time Shift
    double TauMean = 0.0;
    double TauVariance = 0.1;
    auto Tau = std::make_shared< GaussianRandomVariable >(TauMean, TauVariance);
    RandomVariable Tau_("Tau", Tau);
    m_IndividualRandomVariables.insert(Tau_);

    // Initial Space shift coefficient
    double SLocation = 0.0;
    double SScale = 1.0/2.0;
    auto S = std::make_shared< LaplaceRandomVariable >(SLocation, SScale);
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        RandomVariable S_("S" + std::to_string(i), S);
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
    double P0Mean = 0.5;
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
    double V0Variance = 0.01;
    auto V0 = std::make_shared<GaussianRandomVariable>(V0Mean, V0Variance);
    RandomVariable V0_("V0", V0);
    m_PopulationRandomVariables.insert(V0_);

    // Initial Beta coefficient
    double BetaVariance = 0.1;
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        double BetaMean = (double)i;
        auto Beta = std::make_shared< GaussianRandomVariable> (BetaMean, BetaVariance);
        std::string name = "Beta" + std::to_string(i);
        RandomVariable Beta_(name, Beta);
        m_PopulationRandomVariables.insert(Beta_);
    }


    // Noise
    double NoiseMean = 0.0;
    double NoiseVariance = 0.01;
    m_Noise = std::make_shared< GaussianRandomVariable >(NoiseMean, NoiseVariance);
    
}

void
LongitudinalModel
::InitializeIndividualRandomVariables()
{
    // Initial pre acceleration factor
    double KsiMean = 0.0;
    double KsiVariance = 0.1;
    auto Ksi = std::make_shared< GaussianRandomVariable >(KsiMean, KsiVariance);
    RandomVariable Ksi_("Ksi", Ksi);
    m_IndividualRandomVariables.insert(Ksi_);

    // Initial Time Shift
    double TauMean = 0.0;
    double TauVariance = 0.1;
    auto Tau = std::make_shared< GaussianRandomVariable >(TauMean, TauVariance);
    RandomVariable Tau_("Tau", Tau);
    m_IndividualRandomVariables.insert(Tau_);

    // Initial Space shift coefficient
    double SLocation = 0.0;
    double SScale = 1.0/2.0;
    auto S = std::make_shared< LaplaceRandomVariable >(SLocation, SScale);
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        RandomVariable S_("S" + std::to_string(i), S);
        m_IndividualRandomVariables.insert(S_);
    }

}

void
LongitudinalModel
::InitializeManifoldRandomVariables()
{

    // Initial Propagation coefficient
    double DeltaVariance = 0.1;
    for(int i = 0; i < m_Manifold->GetDimension() - 1; ++i)
    {
        double DeltaMean = 0.5*(double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        std::string name = "Delta" + std::to_string(i);
        RandomVariable Delta_(name, Delta);
        m_ManifoldRandomVariables.insert(Delta_);
    }
}


void
LongitudinalModel
::ComputeOrthonormalBasis(const std::shared_ptr<Realizations>& R)
{
    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    double T0 = R->at("T0")[0];
    std::vector<double> P0, V0;
    for(int i = 0; i < m_Manifold->GetDimension(); ++i)
    {
        P0.push_back(R->at("P0")[0]);
        V0.push_back(R->at("V0")[0]);
    }

    /// Compute the transformation to do the Householder reflection in a Euclidean space
    std::vector<double> U = m_Manifold->GetVelocityTransformToEuclideanSpace(P0, T0, V0);

    /// Compute the initial pivot vector U
    double Norm = 0;
    for(std::vector<double>::iterator it = U.begin() ; it != U.end() ; ++it)
    {
        Norm += *it * *it;
    }
    Norm = sqrt(Norm);

    U[0] += copysign(1, -U[0])*Norm;

    double NormU = 0;
    for(std::vector<double>::iterator it = U.begin(); it != U.end(); ++it)
    {
        NormU += *it * *it;
    }

    /// Compute Q = I(N) - 2U . Ut / (Ut . U)
    std::vector< std::vector< double > > Q;
    for(std::vector<double>::iterator it = U.begin(); it != U.end(); ++it)
    {
        std::vector<double> Coordinate;
        for(std::vector<double>::iterator it2 = U.begin(); it2 != U.end(); ++it2)
        {
            Coordinate.push_back( - 2 * *it * *it2 / NormU);
        }
        Q.push_back( Coordinate );
    }

    for(int i = 0; i < Q.size() ; ++i)
    {
        Q[i][i] += 1;
    }

    /// COLINEARITY BETWEEN THE FIRST VECTOR AND Q[0] HAS BEEN CHECKED

    /// ORTHOGONALITY BETWEEN BASIS VECTORS HAS BEEN CHECKED


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
            Beta.push_back( R->at( "Beta" + Number)[0]);

        }
        std::vector<double> AColumn = LinearCombination(Beta, m_OrthogonalBasis);
        AMatrix.push_back( AColumn );

        /*
        if(isnan(LinearCombination(Beta, m_OrthogonalBasis)[0]))
        {
            std::cout << "Compute A Matrix : ";
            for(int i = 0; i < m_OrthogonalBasis.size(); ++i)
            {
                for (int j = 0; j < m_OrthogonalBasis[0].size(); ++j)
                {
                    std::cout << m_OrthogonalBasis[i][j] << ". ";

                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
         */


        /*
        /// DEBUGGING METHOD : CHECK IF A COLUMN ORTHOGONAL TO THE GEO DERIVATIVE AT gamma(T0)
        // TODO : the manifold does not return the geodesic or the geodesic derivative!
        double T0 = R->at("T0")[0];
        std::vector<double> Geo = m_Manifold->GetGeodesic(T0, R);
        std::vector<double> GeoDeriv = m_Manifold->GetGeodesicDerivative(T0, R);
        std::cout << "Should be 0 (Ortho to B0) : " << m_Manifold->ComputeScalarProduct(m_OrthogonalBasis[0], GeoDeriv, Geo) << std::endl;
        std::cout << "Should be 0 (Ortho to A ) : " << m_Manifold->ComputeScalarProduct(AColumn, GeoDeriv, Geo) << std::endl;
        */

    }

    m_AMatrix = AMatrix;
}

void
LongitudinalModel
::ComputeSpaceShifts(const std::shared_ptr<Realizations>& R, const int NumberOfSubjects)
{
    std::map< std::string, std::vector<double>> SpaceShifts;

    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        std::vector<double> Si;
        for(int j = 0; j < m_NbIndependentComponents; ++j)
        {
            double Realization = R->at("S" + std::to_string(j))[i];
            Si.push_back(Realization);
        }

        /*
        if(isnan(LinearCombination(Si, m_AMatrix)[0]))
        {
            std::cout << "Linear combin : ";
            for(int i = 0; i < m_AMatrix.size(); ++i)
            {
                for(int j = 0; j < m_AMatrix[0].size(); ++j)
                {
                    std::cout << m_AMatrix[i][j] << ". ";
                }
            }
            std::cout << std::endl;
        }
         */

        std::pair< std::string, std::vector<double> > SpaceShift("W"+std::to_string(i), LinearCombination(Si, m_AMatrix) );
        SpaceShifts.insert( SpaceShift);
    }

    m_SpaceShifts = SpaceShifts;


    ///////////////////////////////
    //// Debugging : Unit tests ///
    ///////////////////////////////
    /*
    /// Get Data
    double T0 = R->at("T0")[0];
    std::vector<double> P0, V0;
    for(int i = 0; i < m_Manifold->GetDimension(); ++i)
    {
        P0.push_back(R->at("P0")[0]);
        V0.push_back(R->at("V0")[0]);
    }

    /// Compute Norm of the parallel transport
    //std::cout << "Args : " << P0[0] << ". " << T0 << ". " << V0[0] << std::endl;
    std::vector<double> Geodesic0 = m_Manifold->ComputeGeodesic(P0, T0, V0, T0);
    std::vector<double> GeodesicDerivative0 = m_Manifold->ComputeGeodesicDerivative(P0, T0, V0, T0);
    std::vector<double> OK = m_Manifold->ComputeParallelTransport(P0, T0, V0, m_SpaceShifts["W0"], T0);

    double ScalarProduct = m_Manifold->ComputeScalarProduct(OK, GeodesicDerivative0, Geodesic0);
    std::cout << "<diff(g(0)) , W0>|g(0) : " << ScalarProduct << std::endl;

    for(double t = 67; t < 73; t += 1.0)
    {
        std::vector<double> Geodesic = m_Manifold->ComputeGeodesic(P0, T0, V0, t);
        std::vector<double> GeodesicDerivative = m_Manifold->ComputeGeodesicDerivative(P0, T0, V0, t);
        std::vector<double> ParallelTransport = m_Manifold->ComputeParallelTransport(P0, T0, V0, m_SpaceShifts["W0"], t);
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