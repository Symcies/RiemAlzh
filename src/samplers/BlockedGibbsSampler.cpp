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
   
    /// Count number of each random variable    
    unsigned int NbDelta = 0, NbBeta = 0, NbS = 0, NbNu = 0;
    for(auto it = R.begin(); it != R.end(); ++it)
    {
        std::string Name = R.ReverseKeyToName(it->first);
        Name = Name.substr(0, Name.find_first_of("#"));
        if(Name == "Beta")
            NbBeta += 1;
        else if(Name == "Delta")
            NbDelta += 1;
        else if(Name == "S")
            NbS += 1;
        else if(Name == "Nu")
            NbNu += 1;
    }
    
    /// P0 
    auto P0 = std::make_tuple("P0", -1);
    m_Blocks.push_back({P0});
        
    
    /// Beta
    Block BetaPop;
    int SizeOfBetaBlocks = (int)NbBeta;
    for(unsigned int i = 0; i < NbBeta; ++i) 
    {
        auto Beta = std::make_tuple("Beta#" + std::to_string(i), -1);
        
        BetaPop.push_back(Beta);
        if((i == NbBeta - 1 || i%SizeOfBetaBlocks == 0) && i != 0)
        {
            m_Blocks.push_back(BetaPop);
            BetaPop.clear();
        }
    }
    
    
    /// Delta
    Block DeltaPop;
    int SizeOfDeltaBlocks = (int)NbDelta;
    for(unsigned int i = 1; i < NbDelta + 1; ++i) 
    {
        auto Delta = std::make_tuple("Delta#" + std::to_string(i), -1);
        DeltaPop.push_back(Delta);
        
        if(i == NbDelta - 1)
        {
            auto DeltaEnd = std::make_tuple("Delta#" + std::to_string(i + 1), -1);
            DeltaPop.push_back(DeltaEnd);
            m_Blocks.push_back(DeltaPop);
            break;
        }
        
        if(i%SizeOfDeltaBlocks == 0)
        {
            m_Blocks.push_back(DeltaPop);
            DeltaPop.clear();
        }
        
    }

    /// Block Nu for the FastNetwork with Nu as one random variable 
    Block NuPop;
    int SizeOfNuBlocks = (int)NbNu/5;
    for(unsigned int i = 0; i < NbNu; ++i)
    {
        auto Nu = std::make_tuple("Nu#" + std::to_string(i), -1);
        
        NuPop.push_back(Nu);
        if((i == NbNu - 1 || i%SizeOfNuBlocks == 0) && i != 0)
        {
            m_Blocks.push_back(NuPop);
            NuPop.clear();
        }
    }
    
    
    /// Individual Variables
    for(unsigned int i = 0; i < R.at("Ksi").size(); ++i)
    {
        auto Ksi = std::make_tuple("Ksi", i);
        auto Tau = std::make_tuple("Tau", i);
        Block IndividualBlock = {Tau, Ksi};
        for(int j = 0; j < NbS; ++j) 
        {
            auto S = std::make_tuple("S#" + std::to_string(j), i);
            IndividualBlock.push_back(S);
        }
        
        m_Blocks.push_back(IndividualBlock);
    }
    
    M.UpdateModel(R, -1);
    VectorType LL = ComputeLogLikelihood(M, D);
    UpdateLastLogLikelihood(LL);

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
    Block CurrentBlock = m_Blocks[BlockNumber];
    m_CurrentBlockType = TypeRandomVariables(CurrentBlock);
    m_CurrentBlockParameters.clear();
    m_RecoverParameters.clear();
    
    /// Loop over the realizations of the block to update the ratio and the realizations
    double AcceptationRatio = ComputePriorRatioAndUpdateRealizations(R, M, CurrentBlock);
    
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
    UpdateBlockRandomVariable(AcceptationRatio, CurrentBlock);
    
}

ScalarType 
BlockedGibbsSampler
::ComputePriorRatioAndUpdateRealizations(Realizations& R, const AbstractModel &M, const Block &Variables) 
{
    double AcceptanceRatio = 0;
    
    for(auto it = Variables.begin(); it != Variables.end(); ++it) 
    {
        /// Initialization
        std::string NameRealization = std::get<0>(*it);
        unsigned int RealizationNumber = std::max(std::get<1>(*it), 0);
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


int  
BlockedGibbsSampler
::TypeRandomVariables(Block& B) 
{
    int SubjectNumber = std::get<1>(B[0]);
    
    for(auto it = B.begin() + 1; it != B.end(); ++it)
    {
        if(std::get<1>(*it) != SubjectNumber)
        {
            return -1;
        }
    }
    
    return SubjectNumber;
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
::UpdateBlockRandomVariable(double AcceptanceRatio, const Block &Variables) 
{
    double Epsilon = DecreasingStepSize(m_CurrentIteration, m_MemorylessSamplingTime);
    double Denom = m_ExpectedAcceptanceRatio;
    
    if(AcceptanceRatio > m_ExpectedAcceptanceRatio)
        Denom = 1 - m_ExpectedAcceptanceRatio;
    
    for(auto it = Variables.begin(); it != Variables.end(); ++it)
    {
        /// Initialize
        std::string NameRealization = std::get<0>(*it);
        unsigned int RealizationNumber = std::max(std::get<1>(*it), 0);
        
        /// Compute newvariance
        GaussianRandomVariable GRV = m_CandidateRandomVariables.GetRandomVariable(NameRealization, RealizationNumber);
        double CurrentVariance = GRV.GetVariance();
        double NewVariance = CurrentVariance * (1 + Epsilon * 100 * (AcceptanceRatio - m_ExpectedAcceptanceRatio) / Denom );
        NewVariance = std::max(NewVariance, 0.0000001);
        
        // TO CHANGE LIKE THE MASTER BRANCH
        UpdatePropositionDistributionVariance(GRV, AcceptanceRatio, m_CurrentIteration);
        /// Set new variance
        //m_CandidateRandomVariables.UpdatePropositionVariableVariance(NameRealization, RealizationNumber, NewVariance);
        
    }
}