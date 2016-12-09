#include "BlockedGibbsSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


BlockedGibbsSampler
::BlockedGibbsSampler() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
BlockedGibbsSampler
::InitializeSampler(const std::shared_ptr<MultiRealizations> &R) 
{
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R);
    
    //////////////////
    ///   Test 1   ///
    //////////////////
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
    
    
    
    
    //////////////////
    ///   Test 2   ///
    //////////////////
    
    unsigned int NbDelta = 0, NbBeta = 0;
    for(auto it = R->begin(); it != R->end(); ++it)
    {
        std::string Name = it->first;
        Name = Name.substr(0, Name.find_first_of("#"));
        if(Name == "Beta")
        {
            NbBeta += 1;
        }
        if(Name == "Delta")
        {
            NbDelta += 1;
        }
    }
        
    
    /// Population Variables
    auto P0 = std::make_tuple("P0", -1);
    m_Blocks.push_back({P0});
    
    /// Beta
    Block BetaPop;
    int SizeOfBetaBlocks = (int)NbBeta/3;
    for(unsigned int i = 0; i < NbBeta; ++i) 
    {
        auto Beta = std::make_tuple("Beta#" + std::to_string(i), -1);
        BetaPop.push_back(Beta);
        if(i == NbBeta - 1 || i%SizeOfBetaBlocks == 0)
        {
            m_Blocks.push_back(BetaPop);
            BetaPop.clear();
        }
    }
    //m_Blocks.push_back(BetaPop);
    
    
    /// Delta
    Block DeltaPop;
    int SizeOfDeltaBlocks = (int)NbDelta/5;
    for(unsigned int i = 1; i < NbDelta + 1; ++i) 
    {
        auto Delta = std::make_tuple("Delta#" + std::to_string(i), -1);
        DeltaPop.push_back(Delta);
        if(i == NbDelta - 1 || i%SizeOfDeltaBlocks == 0)
        {
            m_Blocks.push_back(DeltaPop);
            DeltaPop.clear();
        }
    }
    
    
    
    /// Individual Variables
    for(unsigned int i = 0; i < R->at("Ksi").size(); ++i)
    {
        auto Ksi = std::make_tuple("Ksi", i);
        auto Tau = std::make_tuple("Tau", i);
        auto S0 = std::make_tuple("S#0", i);
        auto S1 = std::make_tuple("S#1", i);
        //auto S2 = std::make_tuple("S#2", i);
        // TODO : add also to next block is uncommented
        Block IndividualBlock = {Tau, Ksi, S0, S1};
        m_Blocks.push_back(IndividualBlock);
    }
    
    for(unsigned int i = 0; i < R->at("Tau").size(); ++i)
    {
        auto Tau = std::make_tuple("Tau", i);
        //Block IndividualBlock = {Tau};
        //m_Blocks.push_back(IndividualBlock);
    }
    
    
    ////////////////////////////////////
    ///   Test 3 : Univariate Model  ///
    ////////////////////////////////////
    /*
    /// Population block
    auto P0 = std::make_tuple("P0", -1);
    m_Blocks.push_back({P0});
    
    
    /// Individual variables block
    for(unsigned int i = 0; i < R->at("Ksi").size(); ++i)
    {
        auto Ksi = std::make_tuple("Ksi", i);
        auto Tau = std::make_tuple("Tau", i);
        
        Block IndividualBlock = {Ksi, Tau};
        m_Blocks.push_back(IndividualBlock);
    }
    */

    
     
    //////////////////
    /// Test Model ///
    //////////////////
    /*
    for(unsigned int i = 0; i < R->at("A").size(); ++i)
    {
        auto A = std::make_tuple("A", i);
        auto B = std::make_tuple("B", i);
        Block Pop = {A, B};
        m_Blocks.push_back(Pop);
    }
    
    auto C = std::make_tuple("C", -1);
    m_Blocks.push_back({C});
    */
    
    /////////////////
    /// End Tests ///
    /////////////////
    
    
    
    
    
}


BlockedGibbsSampler::MultiRealizations
BlockedGibbsSampler
::Sample(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
         const std::shared_ptr<Data> &D, int IterationNumber) 
{
    ////////////////////////////////////////
    // TODO : Check if the update is needed
    M->UpdateParameters(R);
    VectorType LL = ComputeLogLikelihood(-1, R, M, D);
    UpdateLastLogLikelihood(-1, LL);
    ////////////////////////////////////////
    
    
    auto NewRealizations = std::make_shared<MultiRealizations>(*R);
    
    for(int i = 0; i < m_Blocks.size(); ++i)
    {
        NewRealizations = std::make_shared<MultiRealizations>(OneBlockSample(i, NewRealizations, M, D, IterationNumber));
    }
    
        
    return *NewRealizations;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


BlockedGibbsSampler::MultiRealizations
BlockedGibbsSampler
::OneBlockSample(int BlockNumber, const std::shared_ptr<MultiRealizations>& R,
                 std::shared_ptr<AbstractModel> &M, const std::shared_ptr<Data> &D,
                 int IterationNumber) 
{
    /// Initialization
    auto NewRealizations = std::make_shared<MultiRealizations>(*R);
    Block CurrentBlock = m_Blocks[BlockNumber];
    double AcceptationRatio = 0;
    std::vector<std::string> CurrentParameters; 
    
    
    /// Loop over the realizations of the block to update the ratio and the realizations
    for(auto it = CurrentBlock.begin(); it != CurrentBlock.end(); ++it) 
    {
        /// Initialization
        std::string NameRealization = std::get<0>(*it);
        CurrentParameters.push_back(NameRealization);
        unsigned int SubjectNumber = std::max(std::get<1>(*it), 0);
        
        ///Get the current (for the ratio) and candidate random variables (to sample a candidate)
        ScalarType CurrentRealization = R->at(NameRealization)(SubjectNumber);
        ScalarType CandidaRealization = m_CandidateRandomVariables.GetRandomVariable(NameRealization, SubjectNumber, CurrentRealization).Sample();
        //std::cout << "Name: " << NameRealization << ". Current : " << CurrentRealization << ". Candidate : " << CandidaRealization << std::endl;
        
        /// Get the random variable
        auto RandomVariable = M->GetRandomVariable(NameRealization);
        AcceptationRatio += RandomVariable->LogLikelihood(CandidaRealization);
        AcceptationRatio -= RandomVariable->LogLikelihood(CurrentRealization);
        
        /// Update the NewRealizations
        NewRealizations->at(NameRealization)(SubjectNumber) = CandidaRealization;
    }
    
    /// Compute the likelihood
    int Type = TypeRandomVariables(CurrentBlock);

    AcceptationRatio -= GetPreviousLogLikelihood(Type, M, R, D);
    
    M->UpdateParameters(NewRealizations, CurrentParameters);
    VectorType ComputedLogLikelihood = ComputeLogLikelihood(Type, NewRealizations, M, D);
    AcceptationRatio += ComputedLogLikelihood.sum();
    
    AcceptationRatio = std::min(AcceptationRatio, 0.0);
    AcceptationRatio = exp(AcceptationRatio);
    
    /// Adaptative variances for the realizations
    for(auto it = CurrentBlock.begin(); it != CurrentBlock.end(); ++it)
    {
        /// Initialization
        std::string NameRealization = std::get<0>(*it);
        unsigned int SubjectNumber = std::max(std::get<1>(*it), 0);
        ScalarType CurrentRealization = R->at(NameRealization)(SubjectNumber);
                
        /// Update variance
        GaussianRandomVariable& GRV = m_CandidateRandomVariables.GetRandomVariable(NameRealization, SubjectNumber, CurrentRealization);
        UpdatePropositionDistributionVariance(GRV, AcceptationRatio, IterationNumber);
    }
    
    
    
    /// Return the new realizations
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);
    double UnifSample = Distribution(Generator);
    
        ///  Rejection : Candidate not accepted
    if(UnifSample > AcceptationRatio)
    {
        //std::cout << "NO!!" << std::endl;
        M->UpdateParameters(R, CurrentParameters);
        return *R;
    }
        /// Acceptation : Candidate is accepted
    else
    {
        //std::cout << "YES!!" << std::endl;
        UpdateLastLogLikelihood(Type, ComputedLogLikelihood);
        return *NewRealizations;
    }
    
}


int  
BlockedGibbsSampler
::TypeRandomVariables(Block B) 
{
    int SubjectNumber = std::get<1>(B[0]);
    
    for(auto it = B.begin(); it != B.end(); ++it)
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
::ComputeLogLikelihood(int Type, const std::shared_ptr<MultiRealizations> R,
                       const std::shared_ptr<AbstractModel> M, const std::shared_ptr<Data> D) 
{
    if(Type == -1) 
    {
        VectorType LogLikelihood(D->size(), 0);
        int i = 0;
        for (auto it = LogLikelihood.begin(); it != LogLikelihood.end(); ++it, ++i) {
            *it = M->ComputeIndividualLogLikelihood(R, D, i);

        }
        return LogLikelihood;
    }
    else 
    {
        return VectorType(1, M->ComputeIndividualLogLikelihood(R, D, Type));
    }
}

double
BlockedGibbsSampler
::GetPreviousLogLikelihood(int Type, const std::shared_ptr<AbstractModel> M, 
                           const std::shared_ptr<MultiRealizations> R, const std::shared_ptr<Data> D) 
{
    if(Type == -1)
        return m_LastLikelihoodComputed.sum();
    else
        return m_LastLikelihoodComputed(Type);
}

void 
BlockedGibbsSampler
::UpdateLastLogLikelihood(int Type, VectorType ComputedLogLikelihood) 
{
    if(Type == -1)
        m_LastLikelihoodComputed = ComputedLogLikelihood;
    else
        m_LastLikelihoodComputed(Type) = ComputedLogLikelihood.sum();
}