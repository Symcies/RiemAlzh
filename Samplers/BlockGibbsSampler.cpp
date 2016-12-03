#include "BlockGibbsSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

BlockGibbsSampler
::BlockGibbsSampler() 
{
    
}


BlockGibbsSampler
::~BlockGibbsSampler() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
BlockGibbsSampler
::InitializeSampler(const std::shared_ptr<MultiRealizations> &R) 
{
    
    m_CandidateRandomVariables.InitializeCandidateRandomVariables(R);
    //////////////////
    /// Test Model ///
    //////////////////    
    for(unsigned int i = 0; i < R->at("A").size(); ++i)
    {
        auto A = std::pair<std::string, int>("A", i);
        auto B = std::pair<std::string, int>("B", i);
        Block Pop = {A, B};
        m_Blocks.push_back(Pop);
    }
    auto C = std::pair<std::string, int>("C", -1);
    m_Blocks.push_back({C});
    
    VectorType LL(R->at("A").size());
    m_LogLikelihood = LL;
}


BlockGibbsSampler::MultiRealizations
BlockGibbsSampler
::Sample(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
         const std::shared_ptr<Data> &D, int IterationNumber) 
{
    /// Update the model parameters
    M->UpdateParameters(R, {"All"});
    
    /// Compute loglikelihood
    int i = 0;
    for(auto it = m_LogLikelihood.begin(); it != m_LogLikelihood.end(); ++it, ++i)
    {
        *it = M->ComputeIndividualLogLikelihood(R, D, i);   
    }
    
    /// Sample block per block
    auto NewRealizations = MultiRealizations(*R);
    for(int i = 0; i < m_Blocks.size(); ++i)
    {
        MultiRealizations&& ret = BlockSample(i, NewRealizations, M, D, IterationNumber);
        NewRealizations = std::move(ret);
    }
    
    return NewRealizations;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


BlockGibbsSampler::MultiRealizations
BlockGibbsSampler
::BlockSample(int BlockNumber, const MultiRealizations& R, std::shared_ptr<AbstractModel> &M,
              const std::shared_ptr<Data> &D, int IterationNumber) 
{
    /// Initialization
    auto NewRealizations = MultiRealizations(R);
    Block CurrentBlock = m_Blocks[BlockNumber];
    double AcceptanceRatio = 0;
    std::vector<std::string> CurrentParameters;
    
    ///Loop over the block realizations = Update the ratio and the realizations
    for(auto it = CurrentBlock.begin(); it != CurrentBlock.end(); ++it)
    {
        /// Initialization 
        std::string NameRealization = std::get<0>(*it);
        CurrentParameters.push_back(NameRealization);
        int SubjectNumber = std::max(0, std::get<1>(*it));
        
        /// Get the current and candidate random variables
        ScalarType CurrentRealization = R.at(NameRealization)(SubjectNumber);
        ScalarType CandidaRealization = m_CandidateRandomVariables.GetRandomVariable(NameRealization, SubjectNumber, CurrentRealization).Sample();
    
        /// Get the random variable
        auto RandomVariable = M->GetRandomVariable(NameRealization);
        AcceptanceRatio += RandomVariable->LogLikelihood(CandidaRealization);
        AcceptanceRatio -= RandomVariable->LogLikelihood(CurrentRealization);
        
        /// Update the NewRealizations
        NewRealizations.at(NameRealization)(SubjectNumber) = CandidaRealization;
    }
    
   /// Compute the type of likelihood to compute 
    int Type = TypeRandomVariables(CurrentBlock);
    
    AcceptanceRatio -= GetPreviousLogLikelihood(Type);
    auto XXX = std::make_shared<MultiRealizations>(R);
    VectorType NewLogLikelihood = ComputeLikelihood(XXX, M, D, Type);
    AcceptanceRatio += NewLogLikelihood.sum();
    
    AcceptanceRatio = std::min(AcceptanceRatio, 0.0);
    AcceptanceRatio = exp(AcceptanceRatio);
    
    /// Adaptative variances for the realizations
    for(auto it = CurrentBlock.begin(); it != CurrentBlock.end(); ++it)
    {
        /// Initialization
        std::string NameRealization = std::get<0>(*it);
        unsigned int SubjectNumber = std::max(0, std::get<1>(*it));
        ScalarType CurrentRealization = R.at(NameRealization)(SubjectNumber);
                
        /// Update variance
        GaussianRandomVariable& GRV = m_CandidateRandomVariables.GetRandomVariable(NameRealization, SubjectNumber, CurrentRealization);
        UpdatePropositionDistributionVariance(GRV, AcceptanceRatio, IterationNumber);
    }
    
    /// Return the new realizations
    std::random_device RD;
    std::mt19937 Generator(RD());
    std::uniform_real_distribution<double> Distribution(0.0, 1.0);
    double UniformDraw = Distribution(Generator);
    
    /// Rejection
    if(UniformDraw > AcceptanceRatio)
    {
        M->UpdateParameters(std::make_shared<MultiRealizations>(R), CurrentParameters);
        return R;
    }
    else
    {
        //UpdateLogLikelihood(NewLogLikelihood);
        return NewRealizations;
    }
}





BlockGibbsSampler::VectorType
BlockGibbsSampler
::ComputeLikelihood(const std::shared_ptr<MultiRealizations> &R, std::shared_ptr<AbstractModel> &M,
                    const std::shared_ptr<Data> &D, int Type) 
{
    if(Type == -1)
    {
        VectorType LogLikelihood(D->size());
        int i = 0;
        for(auto it = LogLikelihood.begin(); it != LogLikelihood.end(); ++it, ++i)
        {
            *it = M->ComputeIndividualLogLikelihood(R, D, i);
        }
        return LogLikelihood;
    }
    else
    {
        VectorType A(1, M->ComputeIndividualLogLikelihood(R, D, Type));
        return A;
    }
}

double 
BlockGibbsSampler
::GetPreviousLogLikelihood(int Type)
{
    if(Type == -1)
    {
        return m_LogLikelihood.sum();
    }
    else
    {
        return m_LogLikelihood(Type);
    }
}

int 
BlockGibbsSampler
::TypeRandomVariables(Block B) 
{
    int Type = std::get<1>(B[0]);
    
    for(auto it = B.begin() + 1; it != B.end(); ++it)
    {
        if(std::get<1>(*it) != Type)
        {
            return -1;
        }
    }
    
    return Type;
    
}