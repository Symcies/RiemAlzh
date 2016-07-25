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
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
LongitudinalModel
::InitializeRandomVariables()
{
    InitializePopulationRandomVariables();
    InitializeIndividualRandomVariables();
    InitializePropositionRandomVariables();
}

SufficientStatistics
LongitudinalModel
::GetSufficientStatistics(const Realizations& R, const Data& D)
{
    std::vector< std::vector< double >> SufficientStatistics;

    /////////////////////////
    /// Compute S1 and S2
    /////////////////////////
    std::vector< double > S1, S2;
    std::pair<Data, int> 

    /// Get population wide attribute , then, Loop over all the subjects
    double T0 = R.at("T0");

    for(std::pair<Data, int> i(D.begin(), 0) ;
        i.first != D.end() && i.second < D.size() ;
        ++i.first, ++i.second)
    {
        /// Given a particular subject, get its attributes, then, loop over its observation
        std::string KsiI = "Ksi" + i.second; 
        double Ksi = R.at(KsiI);    // TODO : Check if cannot be "Ksi" + (string)i.second; 

        std::string TauI = "Tau" + i.second;
        double Tau = R.at(TauI);
        std::vector< double > SpaceShift = m_SpaceShifts[i.second]

        for(IndividualData::iterator it = i.first.begin() ; it != i.first.end() ; ++it)
        {
            double TimePoint = exp(Ksi) * ( *it.second - T0 - Tau) + T0;

            std::vector<double> ParallelCurve = m_Manifold->ComputeParallelCurve(TimePoint, SpaceShift, R);
            std::vector<double> Observation = *it.first;

            S1.push_back( ComputeEuclidianScalarProduct(ParallelCurve, Observation ) );
            S2.push_back( ComputeEuclidianScalarProduct(ParallelCurve, ParallelCurve) );
        }
    }

    SufficientStatistics.push_back(S1);
    SufficientStatistics.push_back(S2);


    /////////////////////////
    /// Compute S3 and S4
    /////////////////////////
    std::vector< double > S3, S4;

    for(int i = 0; i < D.size() ; ++i )
    {
        std::string KsiI = "Ksi" + i;
        Ksi = R.at(KsiI);   // TODO : Check if cannot be "Ksi" + (string)i.second; 
        std::string TauI = "Tau" + i;
        Tau = R.at(TauI);

        S3.push_back(Ksi * Ksi);
        S4.push_back(Tau * Tau); 
    }

    SufficientStatistics.push_back(S3);
    SufficientStatistics.push_back(S4);

    /////////////////////////
    /// Compute S5, S6 and S7
    /////////////////////////
    std::vector< double > S5, S6, S7;
    S5.push_back( S.at("P0") );
    S6.push_back( S.at("T0") );
    S7.push_back( S.at("V0") );

    SufficientStatistics.push_back(S5);
    SufficientStatistics.push_back(S6);
    SufficientStatistics.push_back(S7);


    /////////////////////////
    /// Compute S8
    /////////////////////////
    std::vector< double > S8;
    for(int i = 0; i<m_Manifold->GetDimension() ; ++i)
    {
        string DeltaI = "Delta" + i; 
        double Delta = R.at(DeltaI) // TODO : try at(Delta + (string)i)
        S8.push_back(Delta)
    }
    
    SufficientStatistics.push_back(S8);


    /////////////////////////
    /// Compute S9
    /////////////////////////
    std::vector< double > S9;
    for(int i = 0; i < (m_Manifold->GetDimension()) * m_NbIndependentComponents; ++i)
    {
        string = BetaI = "Beta" + i;
        double Beta = R.at(BetaI); // TODO : try strings
        S9.push_back(Beta);
    }

    SufficientStatistics.push_back(S9);


    return SufficientStatistics;
    
}



void
LongitudinalModel
::UpdateRandomVariables(const std::vector< std::vector< double >>& SufficientStatistics, const Data& D)
{

    /// Update P0(mean)
    auto P0 = m_PopulationRandomVariables.at("P0");
    P0->SetMean(SufficientStatistics[4][0]);

    /// Update T0(mean)
    auto T0 = m_PopulationRandomVariables.at("T0");
    T0->SetMean(SufficientStatistics[5][0]);

    /// Update V0(mean)
    auto V0 = m_PopulationRandomVariables.at("V0");
    V0->SetMean(SufficientStatistics[6][0]);

    /// Update Delta(k)(mean)
    for(std::pair< std::vector< double >, int> i(SufficientStatistics[7].begin(), 0) ;
        i.first != SufficientStatistics[7].end() && i.second < SufficientStatistics[7].size() ;
        ++i.first, ++i.second)
    {
        string DeltaI = "Delta" + i;
        auto Delta = m_PopulationRandomVariables.at(DeltaI) ;
        Delta->SetMean(*i.first) ;
    }

    /// Update Beta(k)(mean)
    for(std::pair< std::vector< double >, int> i(SufficientStatistics[8].begin(), 0) ;
        i.first != SufficientStatistics[8].end() && i.second < SufficientStatistics.[8].size() ;
        ++i.first, ++i.second )
    {
        string BetaI = "Beta" + i;
        auto Beta = m_PopulationRandomVariables.at(BetaI);
        Beta->SetMean(*i.first);
    }

    /// Update Ksi
    double VarianceKsi = 0;
    for(std:vector< double >::iterator it = SufficientStatistics[2].begin() ; it != SufficientStatistics[2].end() : ++it)
    {
        VarianceKsi += *it;
    }
    VarianceKsi /= SufficientStatistics[2].size(); // TODO : CHECK IF SIZE EQUAL TO NUMBER OF PEOPLE

    auto Ksi = m_IndividualRandomVariables.at("Ksi");
    Ksi->SetVariance(VarianceKsi)


    /// Update Tau
    double VarianceTau = 0;
    for(std::vector<double>::iterator it = SufficientStatistics[3].begin() ; it != SufficientStatistics[3].end() ; ++it)
    {
        VarianceTau += *it;
    }
    VarianceTau /= SufficientStatistics[3].size(); // TODO : CHECK IF SIZE EQUAL TO NUMBER OF PEOPLE

    auto Tau = m_IndividualRandomVariables.at("Tau");
    Tau->SetVariance(VarianveTau0);


    /// Update Uncertainty variance

    /// Sum Yijk²
    double NoiseVariance;
    for(Data::iterator it = D.begin() ; it != D.end() ; ++it)
    {
        for(IndividualData::iterator it2 = it.begin() ; it2 != it.end() ; ++it2)
        {
            for(std::vector<double>::iterator it3 = it2.first.begin() ; it3 != it2.second.begin() ; ++it3)
            {
                NoiseVariance += *it3 * *it3;
            }
        }
    }

    /// Sum -2 S1
    for(std::vector< double >::iterator it = SufficientStatistics[0].begin() ; it != SufficientStatistics[0].end() ; ++it)
    {
        NoiseVariance -= 2* *it;
    }

    /// Sum S2
    for(std::vector< double >::iterator it = SufficientStatistics[1].begin() ; it != SufficientStatistics[1].end(); ++it)
    {
        NoiseVariance += *it;
    }


    m_Noise->SetVariance(NoiseVariance);


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
    double V0Mean = 2.0;
    double V0Variance = 0.1;
    auto V0 = std::make_shared<GaussianRandomVariable>(V0Mean, V0Variance);
    RandomVariable V0_("V0", V0);
    m_PopulationRandomVariables.insert(V0_);

    // Initial Beta coefficient
    double BetaVariance = 0.1;
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        double BetaMean = (double)i;
        auto Beta = std::make_shared< GaussianRandomVariable> (BetaMean, BetaVariance);
        std::string name = "Beta" + i;
        RandomVariable Beta_(name, Beta);
        m_PopulationRandomVariables.insert(Beta_);
    }


    // Initialize Random variables related to the manifold
    RandomVariableMap ManifoldRandomVariables = m_Manifold->GetManifoldRandomVariables();
    for(RandomVariableMap::iterator it = ManifoldRandomVariables.begin() ; it != ManifoldRandomVariables.end(); ++it)
    {
        m_PopulationRandomVariables.insert(*it);
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
    m_PopulationRandomVariables.insert(Ksi_);

    // Initial Time Shift
    double TauMean = 0.0;
    double TauVariance = 0.1;
    auto Tau = std::make_shared< GaussianRandomVariable >(TauMean, TauVariance);
    RandomVariable Tau_("Tau", Tau);
    m_PopulationRandomVariables.insert(Tau_);

    // Initial Space shift coefficient
    double SLocation = 0.0;
    double SScale = 1/2;
    auto S = std::make_shared< LaplaceRandomVariable >(SLocation, SScale);
    RandomVariable S_("S", S);
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        m_PopulationRandomVariables.insert(S_);
    }

}