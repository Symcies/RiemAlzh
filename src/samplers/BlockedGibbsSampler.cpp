#include "BlockedGibbsSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


BlockedGibbsSampler
::BlockedGibbsSampler() 
{
    
}

BlockedGibbsSampler
::BlockedGibbsSampler(unsigned int MemorylessSamplingTime, double ExpectedAcceptanceRatio) 
{
    m_MemorylessSamplingTime = MemorylessSamplingTime;
    m_ExpectedAcceptanceRatio = ExpectedAcceptanceRatio;
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
::InitializeSampler(const Realizations &R, AbstractModel &M, const Data& D) 
{
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R);
   
    /// Count number of each random variable    
    unsigned int NbDelta = 0, NbBeta = 0, NbS = 0, NbNu = 0;
    for(auto it = R.begin(); it != R.end(); ++it)
    {
        std::string Name = it->first;
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
    VectorType LL = ComputeLogLikelihood(-1, R, M, D);
    UpdateLastLogLikelihood(-1, LL);
    
    ////////////////////////////////////////
    /// Corresponds to a HM within Gibbs ///
    ////////////////////////////////////////
    /*
    for(auto it = R->begin(); it != R->end(); ++it)
    {
        for(int i = 0; i < it->second.size(); ++i)
        {
            auto Block = std::make_tuple(it->first, i);
            m_Blocks.push_back({Block});
        }
    }
     */
}


BlockedGibbsSampler::Realizations
BlockedGibbsSampler
::Sample(Realizations& R, AbstractModel& M,
         const Data &D, int IterationNumber) 
{
    ////////////////////////////////////////
    // TODO : Check if the update is needed
    M.UpdateModel(R, -1);
    VectorType LL = ComputeLogLikelihood(-1, R, M, D);
    UpdateLastLogLikelihood(-1, LL);
    ////////////////////////////////////////
    
    
    for (int i = 0; i < m_Blocks.size(); ++i) 
    {
        OneBlockSample(i, R, M, D, IterationNumber);
    }
        
    return R;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
BlockedGibbsSampler
::OneBlockSample(int BlockNumber, Realizations& R, AbstractModel &M, 
                 const Data &D, int IterationNumber) 
{
    /// Initialization
    Block CurrentBlock = m_Blocks[BlockNumber];
    double AcceptationRatio = 0;
    m_CurrentBlockParameters.clear();
    m_RecoverParameters.clear();
    
    AcceptationRatio = ComputePriorRatioAndUpdateRealizations(R, M, CurrentBlock);
    
    /// Compute the likelihood
    int Type = TypeRandomVariables(CurrentBlock);

    double PreviousLikelihood = GetPreviousLogLikelihood(Type, M, R, D);
    AcceptationRatio -= PreviousLikelihood;

    M.UpdateModel(R, Type, m_CurrentBlockParameters);
    VectorType ComputedLogLikelihood = ComputeLogLikelihood(Type, R, M, D);
    double NewLikelihood = ComputedLogLikelihood.sum();
    AcceptationRatio += NewLikelihood;
    
    AcceptationRatio = std::min(AcceptationRatio, 0.0);
    AcceptationRatio = exp(AcceptationRatio);
    
    
    /// Adaptative variances for the realizations
    for(auto it = CurrentBlock.begin(); it != CurrentBlock.end(); ++it)
    {
        /// Initialization
        std::string NameRealization = std::get<0>(*it);
        unsigned int SubjectNumber = std::max(std::get<1>(*it), 0);
        ScalarType CurrentRealization = R.at(NameRealization)(SubjectNumber);
                
        /// Update variance
        GaussianRandomVariable GRV = m_CandidateRandomVariables.GetRandomVariable(NameRealization, SubjectNumber);
        GRV.SetMean(CurrentRealization);
        UpdatePropositionDistributionVariance(GRV, AcceptationRatio, IterationNumber);
    }
    
    
    
    /// Return the new realizations
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);
    double UnifSample = Distribution(Generator);
    
    ///  Rejection : Candidate not accepted
    if(UnifSample > AcceptationRatio)
    {
        for(auto it = m_RecoverParameters.begin(); it != m_RecoverParameters.end(); ++it)
        {
            R.at(it->first)(it->second.first) = it->second.second;
        }
           
        M.UpdateModel(R, Type, m_CurrentBlockParameters);
    }
        /// Acceptation : Candidate is accepted
    else
    {
        UpdateLastLogLikelihood(Type, ComputedLogLikelihood);
    }
    
}

ScalarType 
BlockedGibbsSampler
::ComputePriorRatioAndUpdateRealizations(Realizations &R, const AbstractModel &M, const Block &Variables) 
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
        ScalarType CurrentRealization = R.at(NameRealization)(RealizationNumber);
        m_RecoverParameters[NameRealization] = {RealizationNumber, CurrentRealization};
        
        
        /// Get a candidate realization
        auto CandidateRandomVariable = m_CandidateRandomVariables.GetRandomVariable(NameRealization, RealizationNumber);
        CandidateRandomVariable.SetMean(CurrentRealization);
        ScalarType CandidateRealization = CandidateRandomVariable.Sample();
        
        /// Calculate the acceptance ratio
        AcceptanceRatio += CurrentRandomVariable->LogLikelihood(CandidateRealization);
        AcceptanceRatio -= CurrentRandomVariable->LogLikelihood(CurrentRealization);
        
        /// Update the NewRealizations
        R.at(NameRealization)(RealizationNumber) = CandidateRealization;
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
::ComputeLogLikelihood(int Type, const Realizations& R, AbstractModel& M, const Data& D) 
{
      
    if(Type == -1) 
    {
        VectorType LogLikelihood(D.size(), 0);
        int i = 0;
        for (auto it = LogLikelihood.begin(); it != LogLikelihood.end(); ++it, ++i) 
        {
            *it = M.ComputeIndividualLogLikelihood(R, D, i);
        }
        return LogLikelihood;
    }
    else 
    {
        return VectorType(1, M.ComputeIndividualLogLikelihood(R, D, Type));
    }
}

double
BlockedGibbsSampler
::GetPreviousLogLikelihood(int Type, AbstractModel& M, const Realizations& R, const Data& D) 
{
    if(Type == -1)
        return m_LastLikelihoodComputed.sum();
    else
        return m_LastLikelihoodComputed(Type);
}

void 
BlockedGibbsSampler
::UpdateLastLogLikelihood(int Type, VectorType& ComputedLogLikelihood) 
{
    if(Type == -1)
        m_LastLikelihoodComputed = ComputedLogLikelihood;
    else
        m_LastLikelihoodComputed(Type) = ComputedLogLikelihood.sum();
}