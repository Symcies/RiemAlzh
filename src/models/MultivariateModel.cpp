#include "MultivariateModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


MultivariateModel
::MultivariateModel(io::ModelSettings &MS) 
{
    m_ManifoldDimension = MS.GetManifoldDimension();
    m_NbIndependentSources = MS.GetNumberOfIndependentSources();
    
}

MultivariateModel
::~MultivariateModel() 
{
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultivariateModel
::Initialize(const OldData &D) 
{

}

ScalarType
MultivariateModel
::InitializePropositionDistributionVariance(std::string Name) 
const 
{
    
}

void
MultivariateModel
::UpdateModel(const Realizations &R, int Type,
              const std::vector<std::string, std::allocator<std::string>> Names)
{
    
}

AbstractModel::SufficientStatisticsVector
MultivariateModel
::GetSufficientStatistics(const Realizations &AR, const OldData &D) 
{
    
}


void 
MultivariateModel
::UpdateRandomVariables(const SufficientStatisticsVector &StochSufficientStatistics,
                        const OldData &D) 
{
    
}


double
MultivariateModel
::ComputeLogLikelihood(const OldData &D) 
{
    
}


double
MultivariateModel
::ComputeIndividualLogLikelihood(const OldData &D, const int SubjectNumber) 
{
    
}


AbstractModel::OldData
MultivariateModel
::SimulateData(io::DataSettings &DS) 
{
    
}

std::vector<AbstractModel::SamplerBlock>
MultivariateModel
::GetSamplerBlocks() 
const 
{
    
}

AbstractModel::VectorType
MultivariateModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MultivariateModel
::DisplayOutputs(const Realizations &AR) 
{
    
}

void 
MultivariateModel
::SaveData(unsigned int IterationNumber, const Realizations &R) 
{
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////


void
MultivariateModel
::InitializeFakeRandomVariables() 
{
    
}