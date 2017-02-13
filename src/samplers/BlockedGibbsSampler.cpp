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
::InitializeSampler(Realizations& R, AbstractModel &M, const Data& D) 
{
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R);
    m_CurrentBlockType = -1;
    m_Blocks = M.GetSamplerBlocks();
    
}


void
BlockedGibbsSampler
::Sample(Realizations& R, AbstractModel& M, const Data &D) 
{
    m_CurrentBlockType = -1;
    ////////////////////////////////////////
    // TODO : Check if the update is needed
    M.UpdateModel(R, -1);
    VectorType LL = ComputeLogLikelihood(M, D);
    UpdateLastLogLikelihood(LL);
    ////////////////////////////////////////
    
    
    for (int i = 0; i < m_Blocks.size(); ++i) 
    {
        OneBlockSample(i, R, M, D);
    }
        
    ++m_CurrentIteration;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
BlockedGibbsSampler
::OneBlockSample(int BlockNumber, Realizations& R, AbstractModel &M, const Data &D) 
{
    /// Initialization
    SamplerBlock CurrentBlock = m_Blocks.at(BlockNumber);
    m_CurrentBlockType = CurrentBlock.first;
    m_CurrentBlockParameters.clear();
    m_RecoverParameters.clear();
    
    /// Loop over the realizations of the block to update the ratio and the realizations
    double AcceptationRatio = ComputePriorRatioAndUpdateRealizations(R, M, CurrentBlock.second);
    
    /// Compute the previous log likelihood
    AcceptationRatio -= GetPreviousLogLikelihood();

    /// Compute the candidate log likelihood
    M.UpdateModel(R, m_CurrentBlockType, m_CurrentBlockParameters);
    VectorType ComputedLogLikelihood = ComputeLogLikelihood(M, D);
    AcceptationRatio += ComputedLogLikelihood.sum();
    
    /// Compute the aceceptance ratio
    AcceptationRatio = std::min(AcceptationRatio, 0.0);
    AcceptationRatio = exp(AcceptationRatio);
    
    
    ///  Rejection : Candidate not accepted
    if(m_UniformDistribution(Generator) > AcceptationRatio)
    {
        for(auto it = m_RecoverParameters.begin(); it != m_RecoverParameters.end(); ++it) 
        {
            R.at(it->first, it->second.first) = it->second.second;
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
::ComputePriorRatioAndUpdateRealizations(Realizations& R, const AbstractModel &M, const MiniBlock &Variables) 
{
    double AcceptanceRatio = 0;
    
    for(auto it = Variables.begin(); it != Variables.end(); ++it) 
    {
        /// Initialization
        std::string NameRealization = it->first;
        unsigned int RealizationNumber = it->second;
        m_CurrentBlockParameters.push_back(NameRealization);
        
        /// Get the current realization and recover it
        auto CurrentRandomVariable = M.GetRandomVariable(NameRealization);
        ScalarType CurrentRealization = R.at(NameRealization, RealizationNumber);
        m_RecoverParameters[NameRealization] = {RealizationNumber, CurrentRealization};
        
        
        /// Get a candidate realization
        auto CandidateRandomVariable = m_CandidateRandomVariables.GetRandomVariable(NameRealization, RealizationNumber);
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
::ComputeLogLikelihood(AbstractModel& M, const Data& D) 
{
      
    if(m_CurrentBlockType == -1) 
    {
        VectorType LogLikelihood(D.size(), 0);
        int i = 0;
        for (auto it = LogLikelihood.begin(); it != LogLikelihood.end(); ++it, ++i) 
        {
            *it = M.ComputeIndividualLogLikelihood(D, i);
        }
        return LogLikelihood;
    }
    else 
    {
        return VectorType(1, M.ComputeIndividualLogLikelihood(D, m_CurrentBlockType));
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
        double NewVariance = CurrentVariance * (1 + Epsilon * 100 * (AcceptanceRatio - m_ExpectedAcceptanceRatio) / Denom );
        NewVariance = std::max(NewVariance, 0.0000001);
        
        // TODO : TO CHANGE LIKE THE MASTER BRANCH
        UpdatePropositionDistributionVariance(GRV, AcceptanceRatio, m_CurrentIteration);
        /// Set new variance
        //m_CandidateRandomVariables.UpdatePropositionVariableVariance(NameRealization, RealizationNumber, NewVariance);
        
    }
}