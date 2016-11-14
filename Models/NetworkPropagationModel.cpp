#include "NetworkPropagationModel.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


NetworkPropagationModel
::NetworkPropagationModel(const unsigned int NbIndependentComponents,
                          std::shared_ptr<AbstractManifold> &M,
                          std::shared_ptr<MatrixType> &KernelMatrix,
                          std::shared_ptr<MatrixType> &InterpolationMatrix) 
{
    // TODO : check the input data
    //TestAssert::WarningEquality_Object(InterpolationMatrix.is_square(), true);
    //TestAssert::WarningEquality_Object(InvertKernelMatrix.is_square(), true);
    
    /// Get the data
    m_NbIndependentComponents = NbIndependentComponents;
    m_Manifold = M;
    
    
    m_InvertKernelMatrix = inverse(*KernelMatrix);
    m_InterpolationMatrix = *InterpolationMatrix;
    
    /// Open the output file
    m_OutputParameters.open("Parameters.txt", std::ofstream::out | std::ofstream::trunc);
    
    /// Unitialize the vector Interpolation Coefficient
    m_InterpolationCoefficients.set_size(m_InvertKernelMatrix.columns());
}

NetworkPropagationModel
::~NetworkPropagationModel() 
{
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void 
NetworkPropagationModel
::Initialize()
{
    /// Initialization
    m_PopulationRandomVariables.clear();
    m_IndividualRandomVariables.clear();
    
    auto P0 = std::make_shared< GaussianRandomVariable >(0.35, 0.000001);
    auto T0 = std::make_shared< GaussianRandomVariable >(70.0, 0.000025);
    auto V0 = std::make_shared< GaussianRandomVariable >(0.06, 0.00000001);
    auto Ksi = std::make_shared< GaussianRandomVariable >(0.0, 0.16);
    auto Tau = std::make_shared< GaussianRandomVariable >(0.0, 36.0);
    m_Noise = std::make_shared< GaussianRandomVariable >(0.0, 0.0001);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0) );
    m_PopulationRandomVariables.insert( RandomVariable("T0", T0) );
    m_PopulationRandomVariables.insert( RandomVariable("V0", V0) );
    m_IndividualRandomVariables.insert( RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert( RandomVariable("Tau", Tau));
        
    /// Initial Beta coefficient
    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        auto Beta = std::make_shared< GaussianRandomVariable> ((double)i/2.0 - 0.5 , 0.000001);
        m_PopulationRandomVariables.insert(RandomVariable("Beta#" + std::to_string(i), Beta));
    }
    
    /// Initial Space shift coefficient
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    { 
        auto S = std::make_shared< LaplaceRandomVariable >(0.0, 1.0/2.0);
        m_IndividualRandomVariables.insert(RandomVariable("S#" + std::to_string(i), S));
    }
    
    /// Delta - temporal translation
    for(int i = 0; i < m_InvertKernelMatrix.size() - 1; ++i)
    {
        auto Delta = std::make_shared< GaussianRandomVariable > ((double)i/10, 0.0001);
        m_PopulationRandomVariables.insert( RandomVariable("Delta#" + std::to_string(i), Delta));
    }
    
    /// Tests
    TestAssert::WarningInequality_GreaterThan(1.0, P0->GetMean(), 0.0, "NetworkPropagationModel>Initialize : wrong P0");
    TestAssert::WarningInequality_GreaterThan(m_Noise->GetVariance(), 0.0, "NetworkPropagationModel>Initialize : wrong noise");
}

void 
NetworkPropagationModel
::UpdateParameters(const std::shared_ptr<MultiRealizations> &R, const std::vector<std::string> Names) 
{
    std::string Name = Names[0];
    Name = Name.substr(0, Name.find_first_of("#"));
    
    if(Name == "Delta")
    {
        ComputeInterpolationCoefficients(R);
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }
    else if(Name == "P0" or Name == "T0" or Name == "V0")
    {   
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }
    else if(Name == "Beta")
    {
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }
    else if(Name == "S")
    {
        ComputeSpaceShifts(R);
    }
    else if(Name == "Ksi" or Name == "Tau" or Name == "Individual")
    {
        // Nothing to do
    }
    else if(Name == "All")
    {
        ComputeInterpolationCoefficients(R);
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    } 
    else
    {
        std::cout << "Should be name in NetworkPropagationModel > Update Parameters  - " << Name << std::endl;
        ComputeOrthonormalBasis(R);
        ComputeAMatrix(R);
        ComputeSpaceShifts(R);
    }
}


NetworkPropagationModel::SufficientStatisticsVector
NetworkPropagationModel
::GetSufficientStatistics(const std::shared_ptr<MultiRealizations> &R,
                          const std::shared_ptr<Data> &D) 
{
    /////////////////////////
    /// Initialization
    /////////////////////////
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")(0);
    VectorType P0(1, R->at("P0")(0) );
    VectorType V0(1, R->at("V0")(0) );
    VectorType Delta = GetPropagationCoefficients(R);
    
    /////////////////////////
    /// Compute S0
    /////////////////////////
    int K = 0;
    double SumOfObservation = 0.0;
    for(const auto& it : *D)
    {
        K += it.size();
        for(auto it2 : it)
        {
            SumOfObservation += it2.first.squared_magnitude();
        }
    }
    VectorType S0(1, SumOfObservation);
    
    /////////////////////////
    /// Compute S1 and S2
    /////////////////////////
    VectorType S1(K), S2(K);
    
    int i = 0;
    auto IterS1 = S1.begin();
    auto IterS2 = S2.begin();
    for(auto Iter = D->begin(); Iter != D->end() && IterS1 != S1.end() && IterS2 != S2.end(); ++Iter, ++i)
    {
        /// Given a particular subject, get its attributes, then, loop over its observation
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        for(auto it : *Iter)
        {
            double TimePoint = SubjectTimePoint(it.second);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Delta);

            *IterS1 = dot_product(ParallelCurve, it.first);
            *IterS2 = ParallelCurve.squared_magnitude();
            
            ++IterS1, ++IterS2;
        }
    }

    /////////////////////////
    /// Compute S3 and S4
    /////////////////////////
    VectorType S3(R->at("Ksi").size()), S4(R->at("Tau").size());
    auto IterKsi = R->at("Ksi").begin();
    auto IterTau = R->at("Tau").begin();
    auto IterS3 = S3.begin();
    auto IterS4 = S4.begin();
    for( ; IterKsi != R->at("Ksi").end() && IterTau != R->at("Tau").end() && IterS3 != S3.end() && IterS4 != S4.end()
            ; ++IterKsi, ++IterTau, ++IterS3, ++IterS4 )
    {
        *IterS3 = *IterKsi * *IterKsi;
        *IterS4 = *IterTau * *IterTau;
    }
    
    /////////////////////////
    /// Compute S8
    /////////////////////////
    VectorType S8(m_Manifold->GetDimension() - 1);
    auto IterS8 = S8.begin();
    for(auto it = Delta.begin() + 1; it != Delta.end() && IterS8 != S8.end(); ++it, ++IterS8)
    {
        *IterS8 = *it;
    }
    
    /////////////////////////
    /// Compute S9
    /////////////////////////
    VectorType S9((m_Manifold->GetDimension() - 1) * m_NbIndependentComponents);
    int j = 0;
    for(auto it = S9.begin(); it != S9.end(); ++it, ++j)
    {
        *it = R->at("Beta#" + std::to_string(j))[0];
    }
    
    SufficientStatisticsVector S = {S0, S1, S2, S3, S4, P0, VectorType(1, T0), V0, S8, S9};
    
    /////////////////////////
    /// Tests
    /////////////////////////
    /*
    /// There are two ways to compute the logLikelihood. Check if they are equal
    std::function<double()> f1 = [=, &D, &R] () { return this->ComputeLogLikelihoodGeneric(R, D); };
    std::function<double()> f2 = [=, &D, &S] ()
    {
        double Calculation = 0;
        int K = 0;
        for(const auto& it : *D)
        {
            K += it.size();
            for(auto it2 : it)
            {
                for(auto it3 : it2.first)
                {
                    Calculation += it3*it3;
                }
            }
        }
        for(auto it : S[0])
        {
            Calculation += -2 * it;
        }
        for(auto it : S[1])
        {
            Calculation += it;
        }
        Calculation /= -2*this->m_Noise->GetVariance();
        Calculation -= K * log(sqrt( 2 * M_PI * m_Noise->GetVariance() ));  
            
        return Calculation;
    };
    
    TestAssert::WarningEquality_Function(f1, f2, "Likelihood not correct. LongitudinalModel > GetSufficientStatistics");
    */
    
    return  S;
}


void 
NetworkPropagationModel
::UpdateRandomVariables(const SufficientStatisticsVector &StochSufficientStatistics,
                        const std::shared_ptr<Data> &D) 
{
    double NumberOfSubjects = D->size();
    
    /// Update P0(mean), T0(mean) and VO(mean)
    auto AbstractP0 = m_PopulationRandomVariables.at("P0");
    auto AbstractT0 = m_PopulationRandomVariables.at("T0");
    auto AbstractV0 = m_PopulationRandomVariables.at("V0");
    
    auto P0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractP0 );
    auto T0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractT0 );
    auto V0 = std::dynamic_pointer_cast< GaussianRandomVariable >( AbstractV0 );
    
    P0->SetMean(StochSufficientStatistics[5][0]);
    T0->SetMean(StochSufficientStatistics[6][0]);
    V0->SetMean(StochSufficientStatistics[7][0]);

    /// Update Delta(k)(mean)
    int i = 0;
    for(auto it = StochSufficientStatistics[8].begin(); it != StochSufficientStatistics[8].end(); ++it, ++i)
    {
        auto AbstractDelta = m_PopulationRandomVariables.at( "Delta#" + std::to_string(i) );
        auto Delta = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractDelta ) ;
        Delta->SetMean(*it) ;
    }
    
    /// Update Beta(k)(mean)
    i = 0;
    for(auto it = StochSufficientStatistics[9].begin(); it != StochSufficientStatistics[9].end(); ++it, ++i)
    {
        auto AbstractBeta = m_PopulationRandomVariables.at( "Beta#" + std::to_string(i) );
        auto Beta = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractBeta );
        std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractBeta );
        Beta->SetMean(*it);
    }
    
    /// Update Ksi & Tau
    double VarianceKsi = 0, VarianceTau = 0;
    
    auto IterKsi = StochSufficientStatistics[3].begin();
    auto IterTau = StochSufficientStatistics[4].begin();
    for( ; IterKsi != StochSufficientStatistics[3].end() && IterTau != StochSufficientStatistics[4].end(); ++IterKsi, ++IterTau)
    {
        VarianceKsi += *IterKsi;
        VarianceTau += *IterTau;
    }
    VarianceKsi /= NumberOfSubjects;
    VarianceTau /= NumberOfSubjects;
    
    auto AbstractKsi = m_IndividualRandomVariables.at("Ksi");
    auto AbstractTau = m_IndividualRandomVariables.at("Tau");
    
    auto Ksi = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractKsi );
    auto Tau = std::dynamic_pointer_cast<GaussianRandomVariable>( AbstractTau );
    
    Ksi->SetVariance(VarianceKsi);   
    Tau->SetVariance(VarianceTau);

    /// Update Uncertainty variance
    double N = m_Manifold->GetDimension();
    double K = 0;

    /// Sum Yijk²
    double NoiseVariance = StochSufficientStatistics[0](0);
    for(auto it : *D)
    {
        K += it.size();
    }

    /// Sum -2 S1 + S2
    auto IterS1 = StochSufficientStatistics[1].begin();
    auto IterS2 = StochSufficientStatistics[2].begin();
    for( ; IterS1 != StochSufficientStatistics[1].end() && IterS2 != StochSufficientStatistics[2].end(); ++IterS1, ++IterS2)
    {
        NoiseVariance += - 2 * *IterS1 + *IterS2;
    }

    /// Divide by N*K, then take the square root
    NoiseVariance /= N*K;
    m_Noise->SetVariance(NoiseVariance);
    
    
    /// Tests
    TestAssert::WarningInequality_GreaterThan(1.0, StochSufficientStatistics[5](0), 0.0, "LongitudinalModel>UpdateRandomVariables ; wrong P0");
    TestAssert::WarningInequality_GreaterThan(m_Noise->GetVariance(), 0.0, "LongitudinalModel>UpdateRandomVariables ; wrong noise variance");
}


double 
NetworkPropagationModel
::ComputeLogLikelihood(const std::shared_ptr<MultiRealizations> &R, const std::shared_ptr<Data> &D) 
{
     /// Get the data
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")(0);
    VectorType P0(1, R->at("P0")(0) );
    VectorType V0(1, R->at("V0")(0) );
    VectorType Delta = GetPropagationCoefficients(R);
    
    /// Compute the likelihood
    double LogLikelihood = 0, K = 0;
    int i = 0;
    for(auto IterData = D->begin(); IterData != D->end(); ++IterData, ++i)
    {
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        K += IterData->size();
        
        for(auto it : *IterData)
        {
            double TimePoint = SubjectTimePoint(it.second);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Delta);
            LogLikelihood  += (it.first - ParallelCurve).squared_magnitude();
        }
    }
    
    LogLikelihood  /= -2*m_Noise->GetVariance();
    LogLikelihood -= K*log(sqrt(2 * m_Noise->GetVariance() * M_PI ));
    
    return LogLikelihood ;
}

double 
NetworkPropagationModel
::ComputeIndividualLogLikelihood(const std::shared_ptr<MultiRealizations> &R,
                                 const std::shared_ptr<Data> &D, const int SubjectNumber) 
{
    /// Initialize the individual parameters
    std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(SubjectNumber, R);
    VectorType SpaceShift = m_SpaceShifts.at("W" + std::to_string(SubjectNumber));

    /// Get the global parameters
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")(0);
    VectorType P0(1, R->at("P0")(0) );
    VectorType V0(1, R->at("V0")(0) );
    VectorType Delta = GetPropagationCoefficients(R);

    
    /// Compute the likelihood
    double LogLikelihood = 0;
    double k = D->at(SubjectNumber).size();
    int i = 0;
    for(auto IterData = D->at(SubjectNumber).begin(); IterData != D->at(SubjectNumber).end(); ++IterData, ++i)
    {
        double TimePoint = SubjectTimePoint(IterData->second);
        VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Delta);
        LogLikelihood += (IterData->first - ParallelCurve).squared_magnitude();
    }
    
    /// Test
    TestAssert::WarningInequality_GreaterThan(LogLikelihood, 0.0, "LongitudinalModel>ComputeIndividualLogLikelihood ; wrong Likelihood");

    LogLikelihood /= -2*m_Noise->GetVariance();
    LogLikelihood -= k * log(sqrt( 2 * m_Noise->GetVariance() * M_PI));
    return LogLikelihood;
}



NetworkPropagationModel::Data
NetworkPropagationModel
::SimulateData(int NumberOfSubjects, int MinObs, int MaxObs) 
{
    /// Simulate realizations
    auto R = std::make_shared<MultiRealizations>( SimulateRealizations(NumberOfSubjects) );
    
    /// Initialize model attributes
    ComputeInterpolationCoefficients(R);
    ComputeOrthonormalBasis(R);
    ComputeAMatrix(R);
    ComputeSpaceShifts(R);
    
    /// Initialize
    std::random_device RD;
    std::mt19937 RNG(RD());
    std::uniform_int_distribution<int> Uni(MinObs, MaxObs);
    std::normal_distribution<double> ObsDistrib(70.0, sqrt(3.0));
    std::normal_distribution<double> NoiseDistrib(0.0, sqrt(m_Noise->GetVariance()));
    
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")(0);
    VectorType P0(1, R->at("P0")(0) );
    VectorType V0(1, R->at("V0")(0) );
    VectorType Delta = GetPropagationCoefficients(R);
    
    Data D;
    
    /// Simulate the data
    for(int i = 0; i < NumberOfSubjects; ++i) 
    {
        IndividualData InDa;
        
        /// Generate the time points of the subjets
        std::vector<double> TimePoints;
        for(int j = 0; j < Uni(RNG); ++j)
        {
            TimePoints.push_back(ObsDistrib(RNG));
        }
        std::sort(TimePoints.begin(), TimePoints.end());
        std::function<double(double)> SubjectTimePoint = GetSubjectTimePoint(i, R);
        auto SpaceShift = m_SpaceShifts.at("W" + std::to_string(i));
        
        /// Generate observations corresponding to the time points
        for(auto it : TimePoints)
        {
            VectorType Scores(SpaceShift.size());
            
            double TimePoint = SubjectTimePoint(it);
            VectorType ParallelCurve = CastedManifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint, Delta);
            int i = 0;
            for(auto it2 = ParallelCurve.begin(); it2 != ParallelCurve.end(); ++it2, ++i)
            {
                Scores(i) = *it2 + NoiseDistrib(RNG);
            }
            InDa.push_back(std::pair<VectorType, double> (Scores, it));
        }
        
        D.push_back(InDa);
    }
    
    return D;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Output(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
NetworkPropagationModel
::ComputeOutputs() 
{
    double MeanP0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("P0"))->GetMean();
    double MeanT0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("T0"))->GetMean();
    double MeanV0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("V0"))->GetMean();
    double SigmaKsi = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Ksi"))->GetVariance();
    double SigmaTau = std::dynamic_pointer_cast<GaussianRandomVariable>(m_IndividualRandomVariables.at("Tau"))->GetVariance();
    double Beta0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Beta#0"))->GetMean();
    double Beta1 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Beta#1"))->GetMean();
    double Beta2 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Beta#2"))->GetMean();
    double Delta0 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#0"))->GetMean();
    double Delta1 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#1"))->GetMean();
    double Delta2 = std::dynamic_pointer_cast<GaussianRandomVariable>(m_PopulationRandomVariables.at("Delta#2"))->GetMean();
    double Sigma = m_Noise->GetVariance();

    m_OutputParameters << MeanP0 << ", " << MeanT0 << ", " << MeanV0 << ", ";
    m_OutputParameters << SigmaKsi << ", " << SigmaTau << ", " << Sigma << ", " << Beta1 << ", " << Beta2 << ", ";
    m_OutputParameters << Delta0 << ", " << Delta1 << ", " << Delta2 << ", " << std::endl;

    
    std::cout << "Parameters : P0: " << MeanP0 << ". T0: " << MeanT0 << ". V0: " << MeanV0;
    std::cout << ". Ksi: " << SigmaKsi << ". Tau: " << SigmaTau << ". Sigma: " << Sigma;
    std::cout << ". Beta0: " << Beta0 << ". Beta1: " << Beta1 << ". Beta2: " << Beta2 << ". Delta0: " << Delta0;
    std::cout << ". Delta1: " << Delta1 << ". Delta2: " << Delta2 << std::endl;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s)  - should not be used in production, maybe in unit function but better erased:
////////////////////////////////////////////////////////////////////////////////////////////////////

void
NetworkPropagationModel
::InitializeFakeRandomVariables() 
{
    
    /////////////////////////////
    /// Population Parameters ///
    /////////////////////////////
    auto P0 = std::make_shared<GaussianRandomVariable>(0.35, 0.000001);
    auto T0 = std::make_shared<GaussianRandomVariable>(70.0, 0.0001);
    auto V0 = std::make_shared<GaussianRandomVariable>(0.06, 0.000001);
    
    m_PopulationRandomVariables.insert( RandomVariable("P0", P0));
    m_PopulationRandomVariables.insert(RandomVariable("T0", T0));
    m_PopulationRandomVariables.insert(RandomVariable("V0", V0));

    m_Noise = std::make_shared< GaussianRandomVariable >(0.0, 0.0001);

    for(int i = 0; i < m_NbIndependentComponents*(m_Manifold->GetDimension()-1) ; ++i)
    {
        auto Beta = std::make_shared< GaussianRandomVariable> ((double)i/5.0, 0.0001);
        std::string name = "Beta#" + std::to_string(i);
        m_PopulationRandomVariables.insert(RandomVariable(name, Beta));
    }
    
    for(int i = 0; i < m_Manifold->GetDimension() - 1 ; ++i)
    {
        auto Delta = std::make_shared< GaussianRandomVariable >(0.5 + (double)i, 0.0001);
        m_PopulationRandomVariables.insert(RandomVariable("Delta#" + std::to_string(i), Delta));
    }


    /////////////////////////////
    /// Individual Parameters ///
    /////////////////////////////
    auto Ksi = std::make_shared< GaussianRandomVariable >(0.0, 0.5*0.5);
    auto Tau = std::make_shared< GaussianRandomVariable >(0.0, 5.0*5.0);
    
    m_IndividualRandomVariables.insert(RandomVariable("Ksi", Ksi));
    m_IndividualRandomVariables.insert(RandomVariable("Tau", Tau));
    
    for(int i = 0; i<m_NbIndependentComponents ; ++i)
    {
        auto S = std::make_shared< LaplaceRandomVariable >(0.0, 1.0/2.0);
        m_IndividualRandomVariables.insert(RandomVariable("S#" + std::to_string(i), S));
    }
    
    
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

NetworkPropagationModel::VectorType
NetworkPropagationModel
::GetPropagationCoefficients(const std::shared_ptr<MultiRealizations> &R) 
{
    // TO DO : Maybe do not calculate the meshpoints that belong to the control point : they are equal to delta(i)
    VectorType Prop = m_InterpolationMatrix * m_InterpolationCoefficients ;
    return Prop;
    
}

std::function<double(double)>  
NetworkPropagationModel
::GetSubjectTimePoint(const int SubjectNumber, const std::shared_ptr<MultiRealizations>& R)
{
    double AccFactor = exp(R->at("Ksi")(SubjectNumber));
    double TimeShift = R->at("Tau")(SubjectNumber);
    double T0 = R->at("T0")(0);
    
    return [AccFactor, TimeShift, T0](double t) { return AccFactor * (t - TimeShift - T0) + T0; };
}

void
NetworkPropagationModel
::ComputeOrthonormalBasis(const std::shared_ptr<MultiRealizations>& R) {

    /////////////////////////
    /// Get the vectors P0 and V0
    /////////////////////////
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(
            m_Manifold);
    double T0 = R->at("T0")(0);
    VectorType P0(1, R->at("P0")(0));
    VectorType V0(1, R->at("V0")(0));
    VectorType Delta = GetPropagationCoefficients(R);


    /// Compute the transformation to do the Householder reflection in a Euclidean space
    VectorType U = CastedManifold->GetVelocityTransformToEuclideanSpace(P0, T0, V0, Delta);

    /// Compute the initial pivot vector U
    double Norm = U.magnitude();
    U(0) += copysign(1, -U(0)) * Norm;
    double NormU = U.squared_magnitude();

    /// Compute Q = I(N) - 2U . Ut / (Ut . U)
    std::vector<VectorType> Q;
    for (auto it = U.begin(); it != U.end(); ++it) {
        VectorType Coordinate(U.size());
        int i = 0;
        for (auto it2 = U.begin(); it2 != U.end(); ++it2, ++i) {
            Coordinate(i) = -2 * *it * *it2 / NormU;
        }
        Q.push_back(Coordinate);
    }

    for (int i = 0; i < Q.size(); ++i) {
        Q[i](i) += 1;
    }



    /// TESTS
    /*
    /// Test colinearity between the first vector and Q[0]
    
    
    /// Orthogonality between the basis vectors and the velocity
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    for(auto it = Q.begin() + 1; it != Q.end(); ++it)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(*it, InitialVelocity, InitialPosition);
        };
        std::function<double()> f2 = []() { return 0;};
        TestAssert::WarningEquality_Function(f1, f2, "Basis vector B not orthogonal to the velocity. Longotudinal Model > ComputeOrthonormalBasis");
    }   
    
    /// END TEST
    */
    /// Drop the first vector which is colinear to gamma_derivative(t0)

    Q.erase(Q.begin());
    m_OrthogonalBasis = Q;
}

void
NetworkPropagationModel
::ComputeAMatrix(const std::shared_ptr<MultiRealizations>& R)
{
    MatrixType AMatrix(m_Manifold->GetDimension(), m_NbIndependentComponents);
    
    for(int i = 0; i < m_NbIndependentComponents ; ++i)
    {
        VectorType Beta(m_Manifold->GetDimension() - 1);
        for(int j = 0; j < m_Manifold->GetDimension() - 1 ; ++j)
        {
            std::string Number = std::to_string(int(j + i*(m_Manifold->GetDimension() - 1)));
            Beta(j) = R->at( "Beta#" + Number)(0);
        }
        
        VectorType V = LinearCombination(Beta, m_OrthogonalBasis)   ;     
        AMatrix.set_column(i, V);
    }
    
    m_AMatrix = AMatrix;
    
    /*
    /// TESTS
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    
    /// Check if the basis is orthogonal to the velocity
    for(int i = 0; i < m_AMatrix.columns(); ++i)
    {
        auto it = m_AMatrix.get_column(i);
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        std::function<double()> f2 = []() { return 0;};
        TestAssert::WarningEquality_Function(f1, f2, "Basis vector B not orthogonal to the velocity. LongitudinalModel > ComputeAMatrix");
    }
    
    /// Check if the matrix column is orthogonal to the velocity
    for(auto it : m_OrthogonalBasis)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        std::function<double()> f2 = []() { return 0;};
        TestAssert::WarningEquality_Function(f1, f2, "A column not orthogonal to the velocity. LongitudinalModel > ComputeAMatrix");
    }
    /// END TESTS
     */
}

void
NetworkPropagationModel
::ComputeSpaceShifts(const std::shared_ptr<MultiRealizations>& R)
{
    std::map< std::string, VectorType> SpaceShifts;
    int NumberOfSubjects = (int)R->at("Tau").size();
    if(NumberOfSubjects != R->at("Ksi").size())
    {
        throw std::invalid_argument("Not the same number of realization in Tau and in Ksi. Whereas its equal to the number of subjects");
    }

    for(int i = 0; i < NumberOfSubjects; ++i)
    {
        VectorType Si(m_NbIndependentComponents);
        for(int j = 0; j < m_NbIndependentComponents; ++j)
        {
            Si(j) = R->at("S#" + std::to_string(j))(i);
        }

        VectorType V = m_AMatrix * Si;
        std::pair< std::string, VectorType > SpaceShift("W"+std::to_string(i),  V);
        SpaceShifts.insert( SpaceShift);
    }

    m_SpaceShifts = SpaceShifts;

    /*
    ///////////////////////////////
    //// Debugging : Unit tests ///
    ///////////////////////////////
    
    
    /// Get Data
    std::shared_ptr<PropagationManifold> CastedManifold = std::dynamic_pointer_cast<PropagationManifold>(m_Manifold);
    double T0 = R->at("T0")[0];
    auto InitialPosition = GetInitialPosition(R);
    auto InitialVelocity = GetInitialVelocity(R);
    auto Delta = GetPropagationCoefficients(R);
    std::function<double()> NullFunction = []() { return 0;};
    
    /// DEBUG 1 : The orthonormal basis is orthogonal to the velocity
    for(auto it : m_AMatrix)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        
        TestAssert::WarningEquality_Function(f1, NullFunction, "Basis vector B not orthogonal to the velocity. LongitudinalModel > ComputeSpaceShifts");
    }
    
    /// DEBUG 2 : Check if the matrix column is orthogonal to the velocity
    for(auto it : m_OrthogonalBasis)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it, InitialVelocity, InitialPosition);
        };
        TestAssert::WarningEquality_Function(f1, NullFunction, "A column not orthogonal to the velocity. LongitudinalModel > ComputeSpaceShifts");
    }
    
    /// DEBUG 3 : Check if the space shifts are orthogonal to the velocity
    for(auto it : m_SpaceShifts)
    {
        std::function<double()> f1 = [=, &InitialPosition, &InitialVelocity, &it] ()
        {
            return this->m_Manifold->ComputeScalarProduct(it.second, InitialVelocity, InitialPosition);
        };
        TestAssert::WarningEquality_Function(f1, NullFunction, "Space Shifts not orthogonal to the velocity. Longitudinal Model > ComputeSpaceShifts");
    }
    
    /// DEBUG 4 : <diff(g(t)) , ParallelTransport>|g(t) = 0
    /// DEBUG 5 : <ParallelTransport, ParallelTransport>g(t) = Constant
    auto ParallelTransport = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, m_SpaceShifts["W0"], 60.0, Delta);
    auto position = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, 60.0, Delta);
    double Constant = m_Manifold->ComputeScalarProduct(ParallelTransport, ParallelTransport, position);
    std::function<double()> ConstantFunction = [Constant] () { return Constant; };
    
    
    for(double t = 61; t < 80; t += 0.5)
    {
        
        std::function<double()> f1 = [&CastedManifold, &InitialPosition, &InitialVelocity, &T0, &Delta, &SpaceShifts, t] () 
        {
            auto a = CastedManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, t, Delta);
            auto b = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, SpaceShifts["W0"], t, Delta);
            auto pos = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
            return CastedManifold->ComputeScalarProduct(a, b, pos);
        };
    
        std::function<double()> f2= [&SpaceShifts, &CastedManifold, &InitialPosition, &InitialVelocity, &T0, &Delta, t] () 
        {
            auto a = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, SpaceShifts["W0"], t, Delta);
            auto pos = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
            return CastedManifold->ComputeScalarProduct(a, a, pos);
        };
    
        
        
        
        TestAssert::WarningEquality_Function(f1, NullFunction, "Parallel transport not orthogonal to velocity at " + std::to_string(t) + ". LongitudinalModel > ComputeSpaceShifts");
        TestAssert::WarningEquality_Function(f2, ConstantFunction, "Norm of parallel transport not constant over time. LongiudinalModel > ComputeSpaceShifts");
    }
    




    /// DEBUG 6 : <ParallelTransport, ParallelTransport>g(t) = Constant
    
    for(double t = 62; t < 75; t += 0.5)
    {
        std::vector<double> Geodesic = CastedManifold->ComputeGeodesic(InitialPosition, T0, InitialVelocity, t, Delta);
        std::vector<double> GeodesicDerivative = CastedManifold->ComputeGeodesicDerivative(InitialPosition, T0, InitialVelocity, t, Delta);
        std::vector<double> ParallelTransport = CastedManifold->ComputeParallelTransport(InitialPosition, T0, InitialVelocity, m_SpaceShifts["W0"], t, Delta);
        double ScalarProd = m_Manifold->ComputeScalarProduct(ParallelTransport, GeodesicDerivative, Geodesic);
        std::cout << "<diff(g(t)), ParallelTransport(W0)>|g(t)> at " << t << " : " << ScalarProd << std::endl;
        double Norm = m_Manifold->ComputeScalarProduct(ParallelTransport, ParallelTransport, Geodesic);
        std::cout << "Norm at " << t << " : " << Norm << std::endl;
    }
    




    ///////////////////////////////
    ////End Unit tests ///
    ///////////////////////////////
    */

}

void 
NetworkPropagationModel
::ComputeInterpolationCoefficients(const std::shared_ptr<MultiRealizations>& R) 
{
    
    VectorType Delta(m_InvertKernelMatrix.columns(), 0.0);
    
    int i = 0;
    for(auto it = Delta.begin() + 1; it != Delta.end(); ++it, ++i)
    {
        *it = R->at("Delta#" + std::to_string(i))(0);
    }
    
    m_InterpolationCoefficients = m_InvertKernelMatrix * Delta;
}










