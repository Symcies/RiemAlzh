#include "NetworkPropagationModel.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
NetworkPropagationModel
::Initialize(const std::shared_ptr<ControlPoints> &P, const std::shared_ptr<Data> &D)
{
    /// Initialization
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    m_ManifoldRandomVariables.clear();
    
     /// Initial Position
    auto P0 = std::make_shared<GaussianRandomVariable>( 0.35, 0.000001 );
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    
    /// Initial Time
    auto T0 = std::make_shared<GaussianRandomVariable>( 72.0, 0.0001);
    m_PopulationRandomVariables.insert( RandomVariable("T0", T0));
    
    /// Initial Velocity
    auto V0 = std::make_shared<GaussianRandomVariable>( 0.03, 0.000001 );
    m_PopulationRandomVariables.insert( RandomVariable("V0", V0) );
    
    /// Noise
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.005);
        
    /// Individual Time Shift
    auto Tau = std::make_shared< GaussianRandomVariable >(0.0, 3*3);
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
    
    /// Individual pace / Preaceleration factor
    auto Ksi = std::make_shared< GaussianRandomVariable >(0.0, 0.5); 
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
    
    // Delta - temporal translation
    for(int i = 0; i < P->size(); ++i)
    {
        
    }
}
