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
    
    // TODO : CHECK IF IT THE ONLY ONE?
    m_IndividualRandomVariables = {"Ksi", "Tau", "S#0", "S#1", "S#2"};
    
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
    
    
    /// Population Variables
    auto P0 = std::make_tuple("P0", 0);
    auto V0 = std::make_tuple("V0", 0);
    auto T0 = std::make_tuple("T0", 0);
    Block Pop = {P0, T0, V0};
    m_Blocks.push_back(Pop);
    
    /// Beta
    Block BetaPop;
    for(unsigned int i = 0; i < 27; ++i) 
    {
        auto Beta = std::make_tuple("Beta#" + std::to_string(i), 0);
        BetaPop.push_back(Beta);
    }
    m_Blocks.push_back(BetaPop);
    
    
    /// Delta
    Block DeltaPop;
    for(unsigned int i = 0; i < 9; ++i) 
    {
        auto Delta = std::make_tuple("Delta#" + std::to_string(i), 0);
        DeltaPop.push_back(Delta);
    }
    m_Blocks.push_back(DeltaPop);
    
    /// Individual Variables
    for(unsigned int i = 0; i < R->at("Ksi").size(); ++i)
    {
        auto Ksi = std::make_tuple("Ksi", i);
        auto Tau = std::make_tuple("Tau", i);
        auto S0 = std::make_tuple("S#0", i);
        auto S1 = std::make_tuple("S#1", i);
        auto S2 = std::make_tuple("S#2", i);
        
        Block IndividualBlock = {Ksi, Tau, S0, S1, S2};
        m_Blocks.push_back(IndividualBlock);
    }
    /*
    /// Other variables (Ksi, Tau and Sji)
    for(auto it = R->begin(); it != R->end(); ++it)
    {
        std::string Name = it->first;
        std::string N = Name.substr(0, Name.find_first_of("#"));
        
        if(N == "P0" || N == "T0" || N == "V0" || N == "Beta" || N == "Delta")
        {
            continue;
        }
        for(int i = 0; i < it->second.size(); ++i)
        {
            auto Mono = std::make_tuple(it->first, i);
            Block MonoBlock = {Mono};
            m_Blocks.push_back( MonoBlock );
        }
    } 
     */
    
    
    //////////////////
    ///   Test 3   ///
    //////////////////
    
    // TODO : Try to group all the population variables
    // TODO : Try to group each individual 
    // TODO : Send the names of the variables to the model
    // TODO : In the model, analyze the names that have been send. Create a switch for each case !
    
    
    
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
    m_LastLikelihoodComputed = 0;
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
    double ComputedLikelihood = 0;
    double AcceptationRatio = 0;
    std::vector<std::string> CurrentParameters; 
    
    /// Loop over the realizations of the block to update the ratio and the realizations
    for(auto it = CurrentBlock.begin(); it != CurrentBlock.end(); ++it) 
    {
        /// Initialization
        std::string NameRealization = std::get<0>(*it); 
        CurrentParameters.push_back(NameRealization);
        unsigned int SubjectNumber = std::get<1>(*it);
        
        ///Get the current (for the ratio) and candidate random variables (to sample a candidate)
        ScalarType CurrentRealization = R->at(NameRealization)(SubjectNumber);
        ScalarType CandidaRealization = m_CandidateRandomVariables.GetRandomVariable(NameRealization, SubjectNumber, CurrentRealization).Sample();
        
        /// Get the random variable
        auto RandomVariable = M->GetRandomVariable(NameRealization);
        AcceptationRatio += RandomVariable->LogLikelihood(CandidaRealization);
        AcceptationRatio -= RandomVariable->LogLikelihood(CurrentRealization);
        
        /// Update the NewRealizations
        NewRealizations->at(NameRealization)(SubjectNumber) = CandidaRealization;
    }
    
    /// Compute the likelihood
    bool Individual = IndividualRandomVariables(CurrentBlock);
    
    if(Individual)
    {
        
        int SubjectNumber = std::get<1>(CurrentBlock[0]);
        AcceptationRatio -= M->ComputeIndividualLogLikelihood(R, D, SubjectNumber);
        
        M->UpdateParameters(NewRealizations, CurrentParameters);
        AcceptationRatio += M->ComputeIndividualLogLikelihood(NewRealizations, D, SubjectNumber);
        
        /*std::cout << "Individual : " << SubjectNumber << " & Parameters ";
        for(auto it = CurrentParameters.begin(); it != CurrentParameters.end(); ++it)
        {
            std::cout << *it << " - ";
        }
        std::cout << std::endl;*/
    } 
    else
    {
        if(m_LastLikelihoodComputed == 0)
        {
            AcceptationRatio -= M->ComputeLogLikelihood(R, D);
            //std::cout << "LENT" << std::endl;
        }
        else
        {
            AcceptationRatio -= m_LastLikelihoodComputed;
            //std::cout << "RAPIDE" << std::endl;
        }
        
        M->UpdateParameters(NewRealizations, CurrentParameters);
        ComputedLikelihood = M->ComputeLogLikelihood(NewRealizations, D); 
        AcceptationRatio += ComputedLikelihood;
    }
    m_LastLikelihoodComputed = ComputedLikelihood;
    
    /// Return the new realizations
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);
    double UnifSample = Distribution(Generator);
    
    if(log(UnifSample) > AcceptationRatio)
    {
        M->UpdateParameters(R, CurrentParameters);
        return *R;
    }
    else
    {
        return *NewRealizations;
    }
    
}


bool 
BlockedGibbsSampler
::IndividualRandomVariables(Block B) 
{
    int SubjectNumber = std::get<1>(B[0]);
    
    for(auto it = B.begin(); it != B.end(); ++it)
    {
        if(std::get<1>(*it) != SubjectNumber)
        {
            return false;
        }
        std::string Name = std::get<0>(*it);
        if(std::find(m_IndividualRandomVariables.begin(), m_IndividualRandomVariables.end(), Name) != m_IndividualRandomVariables.end())
        {
            continue;
        }
        else
        {
            return false;
        }
    }
    
    return true;
}

