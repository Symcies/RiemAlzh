#include "MeshworkModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


MeshworkModel
::MeshworkModel(const ModelSettings& MS) 
{
    m_ManifoldDimension = MS.GetManifoldDimension();
    m_NbIndependentSources = MS.GetNumberOfIndependentSources();
    
    
    std::string KernelMatrixPath = MS.GetInvertKernelPath();
    std::string InterpolationMatrixPath = MS.GetInterpolationKernelPath();
    
    m_InvertKernelMatrix = ReadData::OpenKernel(KernelMatrixPath).transpose();
    m_InterpolationMatrix = ReadData::OpenKernel(InterpolationMatrixPath);
    
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_Thicknesses.set_size(m_ManifoldDimension);
    m_Deltas.set_size(m_ManifoldDimension);
    
}

MeshworkModel
::~MeshworkModel() 
{
    
}


 ////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MeshworkModel
::Initialize(const Data& D) 
{
    
    
    /// 
    
}

void
MeshworkModel
::UpdateModel(const Realizations &R, int Type, const std::vector<std::string> Names) 
{
    bool ComputeThickness = false;
    bool ComputeDeltas = false;
    bool ComputeBasis = false;
    bool ComputeA = false;
    bool ComputeSpaceShifts = false;
    bool ComputeBlock = false;
    bool ComputeIndividual = false;
    
    
}

MeshworkModel::Data
MeshworkModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    
}

double 
MeshworkModel
::ComputeLogLikelihood(const Realizations &R, const Data &D) 
{
    
}

double 
MeshworkModel
::ComputeIndividualLogLikelihood(const Realizations &R, const Data &D, const int SubjectNumber) 
{
    
}


MeshworkModel::SufficientStatisticsVector
MeshworkModel
::GetSufficientStatistics(const Realizations &R, const Data &D) 
{
    
}


void
MeshworkModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS, const Data &D) 
{
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
MeshworkModel
::DisplayOutputs(const Realizations &R) 
{
    
}


void
MeshworkModel
::SaveData(unsigned int IterationNumber, const Realizations &R) 
{
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MeshworkModel
::InitializeFakeRandomVariables() 
{
    
}