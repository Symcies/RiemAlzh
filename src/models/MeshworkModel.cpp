#include "MeshworkModel.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

MeshworkModel
::MeshworkModel(ModelSettings &MS) 
{
    m_ManifoldDimension = MS.GetManifoldDimension();
    m_NbIndependentSources= MS.GetNumberOfIndependentSources();
    
    std::string KernelMatrixPath = MS.GetInvertKernelPath();
    std::string InterpolationMatrixPath = MS.GetInterpolationKernelPath();

    m_InvertKernelMatrix = ReadData::OpenKernel(KernelMatrixPath).transpose();
    m_InterpolationMatrix = ReadData::OpenKernel(InterpolationMatrixPath);
    
    m_NbControlPoints = m_InvertKernelMatrix.columns();
    m_Thicknesses.set_size(m_ManifoldDimension);
    m_Deltas.set_size(m_ManifoldDimension);
    m_Block1.set_size(m_ManifoldDimension);
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
::Initialize(const Data &D) 
{
    /// Data-related attributes
    m_NumberOfSubjects = D.size();
    std::vector<VectorType> IndividualObservationDate;
    
    double SumObservations = 0.0, K = 0.0;
    unsigned int Indiv = 0;
    for(auto it = D.begin(); it != D.end(); ++it, ++Indiv)
    {
        VectorType IndividualObs(it->size());
        K += it->size();
        int i = 0;
        for(auto it2 = it->begin(); it2 != it->end(); ++it2, ++i)
        {
            SumObservations += it2->first.squared_magnitude();
            IndividualObs(i) = it2->second;
        }
        IndividualObservationDate.push_back(IndividualObs);
    }
    m_IndividualObservationDate = IndividualObservationDate;
    m_SubjectTimePoints = IndividualObservationDate;
    m_SumObservations = SumObservations;
    m_NbTotalOfObservations = K;
    
    
    /// Initialize the realizations
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    
    /// Population variables
    m_Noise = std::make_shared<GaussianRandomVariable>( 0.0, 0.000001 );
    

    m_RandomVariables.AddRandomVariable("P", "Gaussian", {0.1, 0.001 * 0.001});
    m_RealizationsPerRandomVariable["P"] = m_NbControlPoints;
    
    for(size_t i = 1; i < m_NbControlPoints; ++i)
    {
        std::string Name = "Delta#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.0001 * 0.0001});
        m_RealizationsPerRandomVariable[Name] = 1;
    }
    
    for(size_t i = 0; i < m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
    {
        std::string Name = "Beta#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0, 0.0001 * 0.0001});
        m_RealizationsPerRandomVariable[Name] = 1;
    }
    
    /// Individual variables
    m_RandomVariables.AddRandomVariable("Ksi", "Gaussian", {0, 0.000000004});
    m_RealizationsPerRandomVariable.insert({"Ksi", m_NumberOfSubjects});
    
    m_RandomVariables.AddRandomVariable("Tau", "Gaussian", {62, 0.25});
    m_RealizationsPerRandomVariable.insert({"Tau", m_NumberOfSubjects});
        
    for(int i = 0; i < m_NbIndependentSources; ++i)
    {
        std::string Name = "S#" + std::to_string(i);
        m_RandomVariables.AddRandomVariable(Name, "Gaussian", {0.0, 0.5});
        m_RealizationsPerRandomVariable.insert({Name, m_NumberOfSubjects});
    }
    
}


void 
MeshworkModel
::UpdateModel(const Realizations &R, int Type,
              const std::vector<std::string, std::allocator<std::string>> Names) 
{
    bool ComputeThickness = false;
    bool ComputeDelta = false;
    bool ComputeBasis = false;
    bool ComputeA = false;
    bool ComputeSpaceShift = false;
    bool ComputeBlock_ = false;
    
    bool IndividualOnly = (Type > -1);
    
    for(auto it = Names.begin(); it != Names.end(); ++it)
    {
        std::string Name = it->substr(0, it->find_first_of("#"));
        
        if(Name == "None")
        {
            continue;
        }
        else if(Name == "Ksi" or Name == "Tau")
        {
            continue;
        }
        else if(Name == "P")
        {
            IndividualOnly = false;
            ComputeThickness = true;
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_ = true;
            continue;
        }
        else if(Name == "Delta")
        { 
            IndividualOnly = false;
            ComputeDelta = true;
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_ = true;
            continue;
        }
        else if(Name == "Beta")
        {
            IndividualOnly = false;
            ComputeA = true;
            ComputeSpaceShift = true;
            continue;
        }
        else if(Name == "S")
        {
            ComputeSpaceShift = true;
        }
        else if(Name == "All")
        {
            IndividualOnly = false;
            ComputeSubjectTimePoint(R, -1);
            
            ComputeThickness = true;
            ComputeDelta = true;
            ComputeBasis = true;
            ComputeA = true;
            ComputeSpaceShift = true;
            ComputeBlock_ = true;
            break;
        } 
        else
        {
            std::cerr << "The random variable name " << Name << "is unknown to the meshwork model" << std::endl; 
        }
    }
    
    if(IndividualOnly)   ComputeSubjectTimePoint(R, Type);
    if(ComputeThickness) ComputeThicknesses(R);
    if(ComputeDelta)     ComputeDeltas(R);
    if(ComputeBasis)     ComputeOrthonormalBasis();
    if(ComputeA)         ComputeAMatrix(R);
    if(ComputeBlock_)    ComputeBlock();
}


MeshworkModel::SufficientStatisticsVector
MeshworkModel
::GetSufficientStatistics(const Realizations &R, const Data &D) 
{
    /// S1 <- y_ij * eta_ij    &    S2 <- eta_ij * eta_ij
    VectorType S1(m_NbTotalOfObservations), S2(m_NbTotalOfObservations);
    auto itS1 = S1.begin(), itS2 = S2.begin();
    int i = 0;
    for(auto itD = D.begin(); itD != D.end(); ++itD, ++i)
    {        
        int j = 0;
        for(auto itD2 = itD->begin(); itD2 != itD->end(); ++itD2, ++j)
        {
            VectorType P2 = ComputeParallelCurve(i, j);
            *itS1 = dot_product(P2, itD2->first);
            *itS2 = P2.squared_magnitude();
            ++itS1, ++itS2;
        }
    }
    
    /// Sufficient Statistic Ksi and Ksi*Ksi
    VectorType S3 = R.at("Ksi");
    VectorType S4 = R.at("Ksi") % R.at("Ksi");
    
    /// Sufficient statistic Tau and Tau*Tau
    VectorType S5 = R.at("Tau");
    VectorType S6 = R.at("Tau") % R.at("Tau");
    
    /// Sufficient statistic beta_k
    VectorType S7((m_ManifoldDimension-1) * m_NbIndependentSources);
    ScalarType * itS7 = S7.memptr();
    for(size_t i = 0; i < ((m_ManifoldDimension-1) * m_NbIndependentSources); ++i)
        itS7[i] = R.at("Beta#" + std::to_string(i), 0);
    
    /// Sufficient statistic p_k and p_k*p_k
    VectorType S8 = R.at("P");
    VectorType S9 = R.at("P") % R.at("P");
    
    /// Sufficient statistic delta_k
    VectorType S10(m_NbControlPoints - 1);
    ScalarType * itS10 = S10.memptr();
    for(size_t i = 1; i < m_NbControlPoints; ++i)
        itS10[i-1] = R.at("Delta#" + std::to_string(i), 0);
    
    
    /// Return the sufficient statistic vector
    SufficientStatisticsVector S = {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10};
    return S;
}

void 
MeshworkModel
::UpdateRandomVariables(const SufficientStatisticsVector &SS,
                        const Data &D) 
{
    /// Update the noise variance, sigma
    ScalarType NoiseVariance = m_SumObservations;
    const ScalarType * itS1 = SS[0].memptr();
    const ScalarType * itS2 = SS[1].memptr();
    for(size_t i = 0; i < SS[0].size(); ++i)
        NoiseVariance += - 2 * itS1[i] + itS2[i];
    
    NoiseVariance /= m_NbTotalOfObservations * m_ManifoldDimension;
    m_Noise->SetVariance(NoiseVariance);
    
    /// Update Ksi : Mean and Variance
    ScalarType KsiMean = 0.0, KsiVariance = 0.0;
    const ScalarType * itS3 = SS[2].memptr();
    const ScalarType * itS4 = SS[3].memptr();
    
    for(size_t i = 0; i < m_NbControlPoints; ++i) 
    {
        KsiMean     += itS3[i];
        KsiVariance += itS4[i];
    }
        
    KsiMean     /= m_NumberOfSubjects;
    KsiVariance -= m_NumberOfSubjects * KsiMean * KsiMean;
    KsiVariance /= m_NumberOfSubjects;
    
    m_RandomVariables.UpdateRandomVariable("Ksi", {{"Mean", KsiMean}, {"Variance", KsiVariance}});
    
    
    /// Update Tau : Mean and Variance
    ScalarType TauMean = 0.0, TauVariance = 0.0;
    const ScalarType * itS5 = SS[4].memptr();
    const ScalarType * itS6 = SS[5].memptr();
    
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        TauMean     += itS5[i];
        TauVariance += itS6[i];
    }
    
    TauMean     /= m_NumberOfSubjects;
    TauVariance -= m_NumberOfSubjects * TauMean * TauMean;
    TauVariance /= m_NumberOfSubjects;
    
    m_RandomVariables.UpdateRandomVariable("Tau", {{"Mean", TauMean}, {"Variance", TauVariance}});
    
    
    /// Update Beta_k : Mean
    const ScalarType * itS7 = SS[6].memptr();
    for(size_t i = 0; i < SS[6].size(); ++i)
        m_RandomVariables.UpdateRandomVariable("Beta#" + std::to_string(i), {{"Mean", itS7[i]}});
     
    /// Update P_k : Mean and Var
    ScalarType PMean = 0.0, PVariance = 0.0;
    const ScalarType * itS8 = SS[7].memptr();
    const ScalarType * itS9 = SS[8].memptr();
    
    for(size_t i = 0; i < m_NbControlPoints; ++i)
    {
        PMean     += itS8[i];
        PVariance += itS9[i];
    }
    
    PMean     /= m_NumberOfSubjects;
    PVariance -= m_NumberOfSubjects * PMean * PMean;
    PVariance /= m_NumberOfSubjects;
    
    /// Update Delta_k : Mean
    const ScalarType * itS10 = SS[9].memptr();
    for(size_t i = 0; i < m_NbControlPoints - 1; ++i)
        m_RandomVariables.UpdateRandomVariable("Delta#" + std::to_string(i+1), {{"Mean", itS10[i]}});
}

ScalarType 
MeshworkModel
::ComputeLogLikelihood(const Data &D) 
{
    double LogLikelihood = 0;
#pragma omp parallel for reduction(+:LogLikelihood)   
    for(size_t j = 0; j < m_NumberOfSubjects; ++j) 
    {
       double N = D.at(j).size();

        for(size_t i = 0; i < N; ++i)
        {
            auto& it = D.at(j).at(i);
            VectorType ParallelCurve = ComputeParallelCurve(j, i);
            LogLikelihood += (it.first - ParallelCurve).squared_magnitude();
        }
        
    }
    
    LogLikelihood  /= -2*m_Noise->GetVariance();
    LogLikelihood -= m_NbTotalOfObservations*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
    
    return LogLikelihood;
}

ScalarType 
MeshworkModel
::ComputeIndividualLogLikelihood(const Data &D, const int SubjectNumber) 
{
    /// Get the data
    double LogLikelihood = 0;
    auto N = D.at(SubjectNumber).size();
    auto Did = D.at(SubjectNumber);
    
  
    for(size_t i = 0; i < N; ++i)
    {
        auto& it = Did.at(i);
        VectorType P2 = ComputeParallelCurve(SubjectNumber, i);
        LogLikelihood += (it.first - P2).squared_magnitude();
    }
    
    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= N * log(2 * m_Noise->GetVariance() * M_PI) / 2.0;
    
    return LogLikelihood;
}

MeshworkModel::Data
MeshworkModel
::SimulateData(DataSettings &DS) 
{
    typedef std::vector< std::pair< VectorType, double> > IndividualData;
    
    m_NumberOfSubjects = DS.GetNumberOfSimulatedSubjects();
    
    /// Initialize the realizations and simulate them
    m_RealizationsPerRandomVariable["P"] = m_NbControlPoints;
    
    for(int i = 1; i < m_NbControlPoints; ++i)
        m_RealizationsPerRandomVariable["Delta#" + std::to_string(i)] = 1;
    
    
    for(size_t i = 0; i <  m_NbIndependentSources*(m_ManifoldDimension - 1); ++i)
        m_RealizationsPerRandomVariable["Beta#" + std::to_string(i)] = 1;
        
    m_RealizationsPerRandomVariable["Ksi"] = m_NumberOfSubjects;
    m_RealizationsPerRandomVariable["Tau"] = m_NumberOfSubjects;
    
    for(int i = 0; i < m_NbIndependentSources; ++i)
        m_RealizationsPerRandomVariable["S#" + std::to_string(i)] = m_NumberOfSubjects;
    
    auto R = SimulateRealizations();
    
    /// Update the model
    ComputeDeltas(R);
    ComputeThicknesses(R);
    ComputeOrthonormalBasis();
    ComputeAMatrix(R);
    ComputeSpaceShifts(R);
    ComputeBlock();
    
    /// Simulate the data
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(DS.GetMinimumNumberOfObservations(), DS.GetMaximumNumberOfObservations());
    std::uniform_real_distribution<double> ObsDistrib(60, 95);
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    Data D;
    double RealNoise = 0.0;
    
    for(int i = 0; i < m_NumberOfSubjects; ++i)
    {
        IndividualData ID;
        std::vector<double> TimePoints;
        for(int j = 0; j < Uni(RNG); ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));
        }
        std::sort(TimePoints.begin(), TimePoints.end());
                
        VectorType T(TimePoints.size());
        int y = 0;
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it, ++y) 
        {
            T(y) = *it;
        }
        
        m_SubjectTimePoints.push_back(T);
        
        int ObsNb = 0;
        for(auto it = TimePoints.begin(); it != TimePoints.end(); ++it, ++ObsNb)
        {
            VectorType Scores(m_ManifoldDimension);
            VectorType ParallelCurve = ComputeParallelCurve(i, ObsNb);
            for(    auto IterCurve = ParallelCurve.begin(), IterScore = Scores.begin() 
                    ; IterCurve != ParallelCurve.end() && IterScore!= Scores.end()
                    ; ++IterCurve, ++IterScore)
            {
                double Noise = NoiseDistrib(RNG);
                RealNoise += Noise;
                *IterScore = *IterCurve + Noise;
            }
            ID.push_back(std::pair<VectorType, double>(Scores, *it));
        }
        D.push_back(ID);
    }
    
    std::vector<VectorType> IndividualObservationDate;
    double SumObservations = 0.0, K = 0.0;
    unsigned int Indiv = 0;
    for(auto it = D.begin(); it != D.end(); ++it, ++Indiv)
    {
        VectorType IndividualObs(it->size());
        K += it->size();
        int i = 0;
        for(auto it2 = it->begin(); it2 != it->end(); ++it2, ++i)
        {
            SumObservations += it2->first.squared_magnitude();
            IndividualObs(i) = it2->second;
        }
        IndividualObservationDate.push_back(IndividualObs);
    }
    m_IndividualObservationDate = IndividualObservationDate;
    m_SumObservations = SumObservations;
    m_NbTotalOfObservations = K;
    
    std::cout << "Real Likelihood = " << ComputeLogLikelihood(D) << std::endl;
    
    return D;
}


std::vector<AbstractModel::SamplerBlock>
MeshworkModel
::GetSamplerBlocks() 
const
{

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Outputs
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


void
MeshworkModel
::InitializeFakeRandomVariables() 
{
    
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
MeshworkModel
::ComputeSubjectTimePoint(const Realizations &R, const int SubjectNumber) 
{
    
}

void 
MeshworkModel
::ComputeDeltas(const Realizations &R) 
{
    
}


void 
MeshworkModel
::ComputeThicknesses(const Realizations &R) 
{
    
}

void 
MeshworkModel
::ComputeOrthonormalBasis() 
{
    
}

void 
MeshworkModel
::ComputeAMatrix(const Realizations &R) 
{
    
}

void 
MeshworkModel
::ComputeSpaceShifts(const Realizations &R) 
{
    
}

void
MeshworkModel
::ComputeBlock() 
{
    
}


AbstractModel::VectorType
MeshworkModel
::ComputeParallelCurve(int SubjectNumber, int ObservationNumber) 
{
    
}





