#include "BlockedGibbsSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


BlockedGibbsSampler
::BlockedGibbsSampler() 
{
    m_CurrentIteration = 0;
    m_CurrentBlockType = -1;
    m_UniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
}

BlockedGibbsSampler
::BlockedGibbsSampler(unsigned int MemorylessSamplingTime, double ExpectedAcceptanceRatio) 
{
    m_UniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
    m_MemorylessSamplingTime = MemorylessSamplingTime;
    m_ExpectedAcceptanceRatio = ExpectedAcceptanceRatio;
    m_CurrentIteration = 0;
    m_CurrentBlockType = -1;
}

BlockedGibbsSampler
::~BlockedGibbsSampler() 
{
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
BlockedGibbsSampler
::InitializeSampler(Realizations& R, AbstractModel &M) 
{
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R, M);
    m_CurrentBlockType = -1;
    m_Blocks = M.GetSamplerBlocks();
    
}


void
BlockedGibbsSampler
::Sample(Realizations& R, AbstractModel& M, const Observations& Obs) 
{
    m_CurrentBlockType = -1;
    ////////////////////////////////////////
    // TODO : Check if the update is needed
    M.UpdateModel(R, -1);
    VectorType LL = ComputeLogLikelihood(M, Obs);
    UpdateLastLogLikelihood(LL);
    ////////////////////////////////////////
    
    
    for (int i = 0; i < m_Blocks.size(); ++i) 
    {
        OneBlockSample(i, R, M, Obs);
    }
        
    ++m_CurrentIteration;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
BlockedGibbsSampler
::OneBlockSample(int BlockNumber, Realizations& R, AbstractModel &M, const Observations& Obs) 
{
    /// Initialization
    SamplerBlock CurrentBlock = m_Blocks.at(BlockNumber);
    m_CurrentBlockType = CurrentBlock.first;
    m_CurrentBlockParameters.clear();
    m_CurrentBlockParametersBIS.clear();
    m_RecoverParameters.clear();
    
    /// Loop over the realizations of the block to update the ratio and the realizations
    double AcceptationRatio = ComputePriorRatioAndUpdateRealizations(R, M, CurrentBlock.second);
    
    /// Compute the previous log likelihood
    AcceptationRatio -= GetPreviousLogLikelihood();
    
    /// Compute the candidate log likelihood
    M.UpdateModel(R, m_CurrentBlockType, m_CurrentBlockParameters);
    VectorType ComputedLogLikelihood = ComputeLogLikelihood(M, Obs);
    AcceptationRatio += ComputedLogLikelihood.sum();
    
    /// Compute the aceceptance ratio
    AcceptationRatio = std::min(AcceptationRatio, 0.0);
    AcceptationRatio = exp(AcceptationRatio);
        
    ///  Rejection : Candidate not accepted
    if(m_UniformDistribution(Generator) > AcceptationRatio)
    {
        for(auto it = m_RecoverParameters.begin(); it != m_RecoverParameters.end(); ++it) 
        {
            std::string Name = std::get<0>(*it);
            unsigned int RealizationNumber = std::get<1>(*it);
            ScalarType PreviousRealization = std::get<2>(*it);
            R.at(Name, RealizationNumber) = PreviousRealization;
        }
           
        M.UpdateModel(R, m_CurrentBlockType, m_CurrentBlockParameters);
    }
        /// Acceptation : Candidate is accepted
    else
    {
        UpdateLastLogLikelihood(ComputedLogLikelihood);
    }
    
    /// Adaptative variances for the realizations
    UpdateBlockRandomVariable(AcceptationRatio, CurrentBlock .second);
    
}

ScalarType 
BlockedGibbsSampler
::ComputePriorRatioAndUpdateRealizations(Realizations& R, const AbstractModel& M, const MiniBlock& Variables) 
{
    double AcceptanceRatio = 0;
    
    for(auto it = Variables.begin(); it != Variables.end(); ++it) 
    {
        /// Initialization
        std::string NameRealization = it->first;
        int Key = R.ReverseNameToKey(NameRealization);
        unsigned int RealizationNumber = it->second;
        m_CurrentBlockParameters.push_back(NameRealization);
        m_CurrentBlockParametersBIS.push_back(Key);
        
        /// Get the current realization and recover it
        auto CurrentRandomVariable = M.GetRandomVariable(Key);
        ScalarType CurrentRealization = R.at(Key, RealizationNumber);
        //m_RecoverParameters[NameRealization] = {RealizationNumber, CurrentRealization};
        m_RecoverParameters.push_back(std::make_tuple(NameRealization, RealizationNumber, CurrentRealization));
        
        /// Get a candidate realization
        auto CandidateRandomVariable = m_CandidateRandomVariables.GetRandomVariable(Key, RealizationNumber);
        CandidateRandomVariable.SetMean(CurrentRealization);
        ScalarType CandidateRealization = CandidateRandomVariable.Sample();
        
        /// Calculate the acceptance ratio
        AcceptanceRatio += CurrentRandomVariable->LogLikelihood(CandidateRealization);
        AcceptanceRatio -= CurrentRandomVariable->LogLikelihood(CurrentRealization);
        
        /// Update the NewRealizations
        R.at(NameRealization, RealizationNumber) = CandidateRealization;
        
    }
    
    return AcceptanceRatio;
    
}



BlockedGibbsSampler::VectorType
BlockedGibbsSampler
::ComputeLogLikelihood(AbstractModel& M, const Observations& Obs) 
{
      
    if(m_CurrentBlockType == -1) 
    {
        VectorType LogLikelihood(Obs.GetNumberOfSubjects(), 0);
        int i = 0;
        for (auto it = LogLikelihood.begin(); it != LogLikelihood.end(); ++it, ++i) 
        {
            *it = M.ComputeIndividualLogLikelihood(Obs.GetSubjectObservations(i) ,i);
        }
        return LogLikelihood;
    }
    else 
    {
        return VectorType(1, M.ComputeIndividualLogLikelihood(Obs.GetSubjectObservations(m_CurrentBlockType), m_CurrentBlockType));
    }
}

double
BlockedGibbsSampler
::GetPreviousLogLikelihood() 
{
    if(m_CurrentBlockType == -1)
        return m_LastLikelihoodComputed.sum();
    else
        return m_LastLikelihoodComputed(m_CurrentBlockType);
}

void 
BlockedGibbsSampler
::UpdateLastLogLikelihood(VectorType& ComputedLogLikelihood) 
{

    if(m_CurrentBlockType == -1)
        m_LastLikelihoodComputed = ComputedLogLikelihood;
    else
        m_LastLikelihoodComputed(m_CurrentBlockType) = ComputedLogLikelihood.sum();
}


void
BlockedGibbsSampler
::UpdateBlockRandomVariable(double AcceptanceRatio, const MiniBlock &Variables) 
{
    double Epsilon = DecreasingStepSize(m_CurrentIteration, m_MemorylessSamplingTime);
    double Denom = m_ExpectedAcceptanceRatio;
    
    if(AcceptanceRatio > m_ExpectedAcceptanceRatio)
        Denom = 1 - m_ExpectedAcceptanceRatio;
    
    for(auto it = Variables.begin(); it != Variables.end(); ++it)
    {
        /// Initialize
        std::string NameRealization = it->first;
        unsigned int RealizationNumber = it->second;
        
        /// Compute newvariance
        GaussianRandomVariable GRV = m_CandidateRandomVariables.GetRandomVariable(NameRealization, RealizationNumber);
        double CurrentVariance = GRV.GetVariance();
        double NewVariance = CurrentVariance * (1 + Epsilon * (AcceptanceRatio - m_ExpectedAcceptanceRatio) / Denom );
        NewVariance = std::max(NewVariance, 0.0000000000000000000000000000001);
        NewVariance = std::max(NewVariance, CurrentVariance / 50);
        NewVariance = std::min(NewVariance, CurrentVariance * 50);
        
        /// Set new variance
        //m_CandidateRandomVariables.UpdatePropositionVariableVariance(NameRealization, RealizationNumber, NewVariance);
        
    }
}