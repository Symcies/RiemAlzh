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

    // Initial uncertainty
    double UncertaintyMean = 0.0;
    double UncertaintyVariance = 0.01;
    auto Epsilon = std::make_shared< GaussianRandomVariable >(UncertaintyMean, UncertaintyVariance);
    RandomVariable Epsilon_("Epsilon", Epsilon);
    m_PopulationRandomVariables.insert(Epsilon_);

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

