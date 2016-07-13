#include "LongitudinalModel.h"
#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


LongitudinalModel
::LongitudinalModel()
{

    m_Data = new Data;

}


LongitudinalModel
::~LongitudinalModel()
{
    delete m_Data;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Getter(s) and Setter(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<double>
LongitudinalModel
::GetAlgorithmParameters()
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // In the parameters to return:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<double> AlgorithmParameters;

    // Add the initial position mean
    AlgorithmParameters.push_back(m_P0->GetMean());

    // Add the initial time mean
    AlgorithmParameters.push_back(m_T0->GetMean());

    // Add the initial velocity mean
    AlgorithmParameters.push_back(m_V0->GetMean());

    // Add the pre acceleration factor variance
    AlgorithmParameters.push_back(m_PreAccelerationFactor[0]->GetMean());

    // Add the time shift variance
    AlgorithmParameters.push_back(m_TimeShift[0]->GetMean());

    // Add the uncertainty variance
    AlgorithmParameters.push_back(*m_UncertaintyVariance);

    return AlgorithmParameters;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
LongitudinalModel
::Initialize()
{

    // Initialize initial position
    double P0Mean = 0.5;
    double P0Variance = 0.01;
    m_P0 = std::make_shared< GaussianRandomVariable >(P0Mean, P0Variance);

    if(m_P0->GetCurrentState() < 0.0 or m_P0->GetCurrentState() > 1.0)
    {
        throw std::invalid_argument(" P0 should be initialize within ]0,1[");
    }



    // Initialize initial time
    double T0Mean = 2.0;
    double T0Variance = 0.01;
    m_T0 = std::make_shared< GaussianRandomVariable >(T0Mean, T0Variance);


    // Initialize initial velocity
    double V0Mean = 2.0;
    double V0Variance = 0.01;
    m_V0 = std::make_shared< GaussianRandomVariable >(V0Mean, V0Variance);


    // Initialize initial propagation coefficient
    double DeltaVariance = 0.01;
    m_PropagationCoefficient.clear();
    for(int i = 0; i<m_Manifold->GetNumberOfDimension() ; ++i)
    {
        double DeltaMean = (double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        m_PropagationCoefficient.push_back(Delta);
    }

    // Initialize the propagation coefficient in the manifold
    m_Manifold->SetPropagationCoefficient(m_PropagationCoefficient);

    // Initialize initial A Matrix coefficient
    m_AMatrixCoefficient.clear();
    double BetaVariance = 0.01;
    int ArraySize = (m_Manifold->GetNumberOfDimension() - 1)*m_Manifold->GetNumberOfIndependentComponents();
    for(int i = 0; i<ArraySize ; ++i)
    {
        double BetaMean = (double)i;
        auto Beta = std::make_shared< GaussianRandomVariable >(BetaMean, BetaVariance);
        m_AMatrixCoefficient.push_back(Beta);
    }

    // Initialize the first Orthonormal basis (B1, ..., B(n-1)Ns)
    ComputeOrthonormalBasis();

    // Initialize the uncertainty variance
    m_UncertaintyVariance = std::make_shared<double>(0.1);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Subject Specific parameters :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set the number of number of subjects
    m_NumberOfSubjects = (int)m_Data->size();

    // Initialize the pre acceleration factor
    double KsiVariance = 0.02;
    double KsiMean = 0.0;  // Always equal to 0
    m_PreAccelerationFactor.clear();
    for(int i=0; i<m_NumberOfSubjects ; ++i)
    {
        auto Ksi = std::make_shared<GaussianRandomVariable>(KsiMean, KsiVariance);
        m_PreAccelerationFactor.push_back(Ksi);
    }

    // Initialize time shifts
    m_TimeShift.clear();
    double TauVariance = 0.02;
    double TauMean = 0.0; // always equal to 0
    for(int i=0; i<m_NumberOfSubjects ; ++i)
    {
        auto Tau = std::make_shared< GaussianRandomVariable>(TauMean, TauVariance);
        m_TimeShift.push_back(Tau);
    }

    // Initialize space shift coefficient Sij
    m_SpaceShiftCoefficient.clear();
    double Location = 0.0;
    double Scale = 1.0/2.0;

    for(int i = 0 ; i<m_NumberOfSubjects ; ++i)
    {
        std::vector<std::shared_ptr<LaplaceRandomVariable>> SpaceShiftCoefficient;
        for(int  j = 0; j < m_Manifold->GetNumberOfIndependentComponents() ; ++j )
        {
            auto Coordinate = std::make_shared<LaplaceRandomVariable>(Location, Scale);
            SpaceShiftCoefficient.push_back(Coordinate);
        }
        m_SpaceShiftCoefficient.push_back(SpaceShiftCoefficient);
    }

    // Compute the first likelihood
    m_Likelihood = ComputeLikelihood();


}

void
LongitudinalModel
::Update()
{
    ComputeOrthonormalBasis();
    ComputeAMatrix();
    ComputeSpaceShifts();

}

double
LongitudinalModel
::ComputeLikelihood()
{
    double T0 = m_T0->GetCurrentState();
    double P0 = m_P0->GetCurrentState();
    double V0 = m_V0->GetCurrentState();

    double Likelihood = 0.0;

    // Compute Likelihood
    int SubjectNumber = 0;
    for(Data::iterator it = m_Data->begin(); it != m_Data->end() ; ++it)
    {
        double AccFactor = exp(m_PreAccelerationFactor[SubjectNumber]->GetCurrentState());
        double TimeShift = m_TimeShift[SubjectNumber]->GetCurrentState();
        std::vector<double> SpaceShift = m_SpaceShift[SubjectNumber];

        for(std::vector<std::pair<double, std::vector<double>>>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
        {
            double TimePoint = AccFactor * (it2->first - T0 - TimeShift) + T0;
            std::vector<double> ParallelCurve = m_Manifold->ComputeParallelCurve(P0, T0, V0, SpaceShift, TimePoint);
            std::vector<double> Observation = it2->second;

            Likelihood += ObservationDifferenceNorm(Observation, ParallelCurve);
            //std::cout << Likelihood << " - ";
        }

        SubjectNumber +=1;
    }


    if(isnan(Likelihood)) std::cout << "Nan likelihood" << std::endl;
    Likelihood = exp( - 1.0/(2.0 * *m_UncertaintyVariance) * Likelihood);
    return Likelihood;

}

std::vector<double>
LongitudinalModel
::ComputeParallelTransport(int i, double ObservationTimePoint)
{
    double P0 = m_P0->GetCurrentState();
    double T0 = m_T0->GetCurrentState();
    double V0 = m_V0->GetCurrentState();
    double AccFactor = exp(m_PreAccelerationFactor[i]->GetCurrentState());
    double TimeShift = m_TimeShift[i]->GetCurrentState();
    double TimePoint = AccFactor*(ObservationTimePoint - T0 - TimeShift) + T0;

    return m_Manifold->ComputeParallelCurve(P0, T0, V0, m_SpaceShift[i], TimePoint);
}

std::vector<RandomVariableToSample>
LongitudinalModel
::GetRandomVariableToSample()
{
    std::vector<RandomVariableToSample> RVToSample;


    // Sample P0 and its proposition distribution
    double P0CandidateVariance = 0.01;
    auto P0Candidate = std::make_shared<GaussianRandomVariable>(m_P0->GetCurrentState(), P0CandidateVariance);
    RandomVariableToSample P0 (m_P0, P0Candidate);
    RVToSample.push_back(P0);


    // Sample T0 and its proposition distribution
    double T0CandidateVariance = 0.01;
    auto T0Candidate = std::make_shared<GaussianRandomVariable>(m_T0->GetCurrentState(), T0CandidateVariance);
    RandomVariableToSample T0 (m_T0, T0Candidate);
    RVToSample.push_back(T0);


    // Sample V0 and its proposition distribution
    double V0CandidateVariance = 0.01;
    auto V0Candidate = std::make_shared<GaussianRandomVariable>(m_V0->GetCurrentState(), V0CandidateVariance);
    RandomVariableToSample V0 (m_V0, V0Candidate);
    RVToSample.push_back(V0);


    // Sample propagation coefficients and their proposition distribution
    double DeltaCandidateVariance = 0.01;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_PropagationCoefficient.begin() ; it != m_PropagationCoefficient.end() ; ++it)
    {
        auto DeltaCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), DeltaCandidateVariance);
        RandomVariableToSample Delta (*it, DeltaCandidate);
        RVToSample.push_back(Delta);
    }


    // Sample the A Matrix coefficient
    double BetaCandidateVariance = 0.01;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_AMatrixCoefficient.begin() ; it != m_AMatrixCoefficient.end() ; ++it)
    {
        auto BetaCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), BetaCandidateVariance);
        RandomVariableToSample Beta (*it, BetaCandidate);
        RVToSample.push_back(Beta);
    }


    // Sample the pre-acceleration factors
    double KsiCandidateVariance = 0.01;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_PreAccelerationFactor.begin(); it != m_PreAccelerationFactor.end() ; ++it)
    {
        auto KsiCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), KsiCandidateVariance);
        RandomVariableToSample Ksi (*it, KsiCandidate);
        RVToSample.push_back(Ksi);
    }


    // Sample the time shift
    double TauCandidateVariance = 0.01;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_TimeShift.begin() ; it != m_TimeShift.end() ; ++it)
    {
        auto TauCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), TauCandidateVariance);
        RandomVariableToSample Tau (*it, TauCandidate);
        RVToSample.push_back(Tau);
    }


    // Sample the Space shift coefficient
    double SCandidateVariance = 0.01;
    for(std::vector<std::vector<std::shared_ptr< LaplaceRandomVariable>>>::iterator it = m_SpaceShiftCoefficient.begin() ; it != m_SpaceShiftCoefficient.end() ; ++it)
    {
        for(std::vector<std::shared_ptr<LaplaceRandomVariable>>::iterator it2 = it->begin() ; it2 != it->end() ; ++it2)
        {
            auto SCandidate = std::make_shared<GaussianRandomVariable>((*it2)->GetCurrentState(), SCandidateVariance);
            RandomVariableToSample S (*it2, SCandidate);
            RVToSample.push_back(S);
        }
    }

    return RVToSample;

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

Data
LongitudinalModel
::SimulateSpecificData(std::vector<std::vector<double>> TimePoint)
{


    //////////////////////////////////////////////////////////////////////////
    // Model Parameters
    //////////////////////////////////////////////////////////////////////////

    int NumberOfSubjects = (int)TimePoint.size();
    m_NumberOfSubjects= NumberOfSubjects;

    m_P0 = std::make_shared<GaussianRandomVariable>(0.5, 0.01);
    m_T0 = std::make_shared<GaussianRandomVariable>(1.0, 0.01);
    m_V0 = std::make_shared<GaussianRandomVariable>(1.0, 0.01);

    m_PropagationCoefficient.clear();
    for(int i = 0; i< m_Manifold->GetNumberOfDimension(); ++i)
    {
        m_PropagationCoefficient.push_back(std::make_shared<GaussianRandomVariable>((double)i, 0.01));
    }

    m_Manifold->SetPropagationCoefficient(m_PropagationCoefficient);

    m_AMatrixCoefficient.clear();
    for(int i = 0; i<(m_Manifold->GetNumberOfDimension()-1)*m_Manifold->GetNumberOfIndependentComponents() ; ++i)
    {
        m_AMatrixCoefficient.push_back(std::make_shared<GaussianRandomVariable>((double)i, 0.01));
    }

    m_PreAccelerationFactor.clear();
    m_TimeShift.clear();
    for(int i = 0; i<NumberOfSubjects; ++i)
    {
        m_PreAccelerationFactor.push_back(std::make_shared<GaussianRandomVariable>(0.0, 0.01));
        m_TimeShift.push_back(std::make_shared<GaussianRandomVariable>(0.0, 0.01));
    }

    m_UncertaintyVariance = std::make_shared<double>(0.1);

    m_SpaceShiftCoefficient.clear();
    for(int i = 0; i<NumberOfSubjects; ++i)
    {
        std::vector<std::shared_ptr< LaplaceRandomVariable >> Coord;
        for(int j = 0; j<m_Manifold->GetNumberOfDimension(); ++j)
        {
            Coord.push_back(std::make_shared<LaplaceRandomVariable>(0.0, 1.0/2.0));
        }
        m_SpaceShiftCoefficient.push_back(Coord);
    }

    ComputeOrthonormalBasis();
    ComputeAMatrix();
    ComputeSpaceShifts();

    //////////////////////////////////////////////////////////////////////////
    // Simulate Data
    //////////////////////////////////////////////////////////////////////////

    Data D;

    std::random_device RD;
    std::default_random_engine Generator(RD());
    std::normal_distribution<double> Distribution(0, *m_UncertaintyVariance.get());

    int i = 0;
    for(std::vector<std::vector<double>>::iterator it = TimePoint.begin() ; it != TimePoint.end() ; ++it)
    {
        std::vector<std::pair<double, std::vector<double>>> SubjectData;
        for(std::vector<double>::iterator it2 = it->begin() ; it2 != it->end(); ++it2)
        {
            std::vector<double> Observation;

            std::vector<double> Eta = ComputeParallelTransport(i,  *it2);
            for(std::vector<double>::iterator it3 = Eta.begin() ; it3 != Eta.end(); ++it3)
            {
                double Gen = Distribution(Generator);
                Observation.push_back(*it3 + Gen);
            }
            SubjectData.push_back(std::pair<double, std::vector<double>> (*it2, Observation));

        }
        D.push_back(SubjectData);
        i+= 1;
    }


    return D;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected attribute(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
LongitudinalModel
::ComputeOrthonormalBasis()
{
    double P0 = m_P0->GetCurrentState();
    double T0 = m_T0->GetCurrentState();
    double V0 = m_V0->GetCurrentState();
    std::vector<double> GammaDerivative0 = m_Manifold->ComputeGeodesicDerivative(P0, T0, V0, T0);

    bool null = true;
    for(std::vector<double>::iterator it = GammaDerivative0.begin(); it != GammaDerivative0.end() ; ++it)
    {
        if(*it != 0)
        {
            null = false;
            break;
        }
    }
    if(null)
    {
        throw std::invalid_argument(" Gamma derivative equals zero : impossible to calculate the orthogonal basis");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// HouseHolder reflection / transformation to get the first basis
    /// https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    // GammaDerivative0, calculated earlier, is the base/ground vector to use the Householder reflection
    std::vector<std::vector<double>> OrthogonalBasis;

    // NormDerivative is the norm of the ground vector
    double NormDerivative = m_Manifold->ComputeMetric(GammaDerivative0, GammaDerivative0, m_Manifold->ComputeGeodesic(P0, T0, V0, T0));
    NormDerivative = sqrt(NormDerivative*NormDerivative);

    // Change the value of the first coordinate of the ground vector
    GammaDerivative0[0] -= copysign(NormDerivative, GammaDerivative0[0]);

    // Normalize the ground vector
    double GroundVectorNorm = m_Manifold->ComputeMetric(GammaDerivative0, GammaDerivative0, m_Manifold->ComputeGeodesic(P0, T0, V0, T0));
    for(std::vector<double>::iterator it = GammaDerivative0.begin() ; it != GammaDerivative0.end() ; ++it)
    {
        *it = *it/GroundVectorNorm;
    }

    // Compute the Orthogonal basis
    for(std::vector<double>::iterator it = GammaDerivative0.begin() ; it != GammaDerivative0.end(); ++it)
    {
        std::vector<double> BasisCoordinate;
        for(std::vector<double>::iterator it2 = GammaDerivative0.begin() ; it2 != GammaDerivative0.end() ; ++it2)
        {
            BasisCoordinate.push_back(- 2 * *it * *it2);
        }
        OrthogonalBasis.push_back(BasisCoordinate);
    }

    // Add the identity
    for(int i = 0; i<OrthogonalBasis.size() ; ++i)
    {
        OrthogonalBasis[i][i] += 1;
    }

    // Drop the first vector, which is colinear to Gamma0Derivative
    OrthogonalBasis.erase(OrthogonalBasis.begin());


    m_OrthonormalBasis = OrthogonalBasis;

}

void
LongitudinalModel
::ComputeAMatrix()
{

    m_AMatrix.clear();

    int Dim = m_Manifold->GetNumberOfDimension();
    int Indep = m_Manifold->GetNumberOfIndependentComponents();

    for(int j = 0; j<Dim ; ++j)
    {
        for(int i = 0 ; i<Indep ; ++i)
        {
            double val = 0;
            for(int k = 0; k<Dim-1; ++k)
            {
                val += m_AMatrixCoefficient[k+i*(Dim-1)]->GetCurrentState() * m_OrthonormalBasis[k][j];
            }
            m_AMatrix.push_back(val);
        }
    }


}

void
LongitudinalModel
::ComputeSpaceShifts()
{
    m_SpaceShift.clear();

    int Indep = m_Manifold->GetNumberOfIndependentComponents();

    for(int i = 0; i<m_NumberOfSubjects ; ++i)
    {
        std::vector< std::shared_ptr<LaplaceRandomVariable> > SpaceShiftCoefficient = m_SpaceShiftCoefficient[i];
        std::vector<double> IndividualSpaceShift;
        for(int j = 0; j<m_Manifold->GetNumberOfDimension() ; ++j)
        {
            double val = 0.0;
            for(int k = 0; k<Indep ; ++k)
            {
                val += m_AMatrix[k+j*Indep] * SpaceShiftCoefficient[k]->GetCurrentState();
            }
            IndividualSpaceShift.push_back(val);
        }
        m_SpaceShift.push_back(IndividualSpaceShift);
    }

}

double
LongitudinalModel
::ObservationDifferenceNorm(std::vector<double> Observation, std::vector<double> ParallelCurve)
{
    typedef std::vector<double>::iterator doubIter;

    double norm = 0.0;
    for(std::pair<doubIter, doubIter> i(Observation.begin(), ParallelCurve.begin()) ;
            i.first != Observation.end() && i.second != ParallelCurve.end();
            ++i.first, ++i.second)
    {
        //if(isnan(*i.first)) std::cout << " Observation " << *i.second << std::endl;
        norm += (*i.first - *i.second) * (*i.first - *i.second) ;
    }

    return norm;
}

std::vector<std::vector< double >>
LongitudinalModel
::InitializeSufficientStochasticStatistics()
{

    int K = 0;
    for(Data::iterator it = m_Data->begin() ; it != m_Data->end() ; ++it)
    {
        K += it->size();
    }

    std::vector<double> S1(K);
    std::vector<double> S2(K);
    std::vector<double> S3(m_Data->size());
    std::vector<double> S4(m_Data->size());
    std::vector<double> S5(1);
    std::vector<double> S6(1);
    std::vector<double> S7(1);
    std::vector<double> S8(m_Manifold->GetNumberOfDimension());
    std::vector<double> S9((m_Manifold->GetNumberOfDimension() - 1)*m_Manifold->GetNumberOfIndependentComponents());

    std::vector<std::vector<double>> S;
    S.push_back(S1);
    S.push_back(S2);
    S.push_back(S3);
    S.push_back(S4);
    S.push_back(S5);
    S.push_back(S6);
    S.push_back(S7);
    S.push_back(S8);
    S.push_back(S9);

    return S;

}


std::vector<std::vector<double>>
LongitudinalModel
::ComputeSufficientStatistics()
{
    std::vector<std::vector<double>> SufficientStatistics;

    // Compute S1 and S2
    std::vector<double> S1;
    std::vector<double> S2;

    int i = 0;
    for(Data::iterator it = m_Data->begin() ; it != m_Data->end() ; ++it)
    {
        for(std::vector<std::pair<double, std::vector<double>>>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
        {
            std::vector<double> ParallelTransport = ComputeParallelTransport(i, it2->first);
            double SumS1 = 0;
            double SumS2 = 0;
            int j = 0;
            for(std::vector<double>::iterator it3 = ParallelTransport.begin() ; it3 != ParallelTransport.end() ; ++it3)
            {
                SumS1 += *it3 * it2->second[j];
                SumS2 += *it3 * *it3;
                j += 1;
            }
            S1.push_back(SumS1);
            S2.push_back(SumS2);
        }
        i += 1;
    }
    SufficientStatistics.push_back(S1);
    SufficientStatistics.push_back(S2);


    // Compute S3
    std::vector<double> S3;

    std::vector<std::shared_ptr< GaussianRandomVariable >>::iterator it;
    for(it = m_PreAccelerationFactor.begin() ; it != m_PreAccelerationFactor.end() ; ++it)
    {
        double vect = it->get()->GetCurrentState() * it->get()->GetCurrentState();
        S3.push_back(vect);
    }
    SufficientStatistics.push_back(S3);


    // Compute S4
    std::vector<double> S4;

    std::vector<std::shared_ptr< GaussianRandomVariable >>::iterator it2;
    for(it2 = m_TimeShift.begin() ; it2 != m_TimeShift.end() ; ++it2)
    {
        double vect = it2->get()->GetCurrentState() * it2->get()->GetCurrentState();
        S4.push_back(vect);
    }
    SufficientStatistics.push_back(S4);


    // Compute S5, S6 and S7
    std::vector<double> S5, S6, S7;
    S5.push_back(m_P0->GetCurrentState());
    S6.push_back(m_T0->GetCurrentState());
    S7.push_back(m_V0->GetCurrentState());

    SufficientStatistics.push_back(S5);
    SufficientStatistics.push_back(S6);
    SufficientStatistics.push_back(S7);


    //Compute S8
    std::vector<double> S8;
    std::vector<std::shared_ptr< GaussianRandomVariable >>::iterator it3;
    for(it3 = m_PropagationCoefficient.begin() ; it3 != m_PropagationCoefficient.end() ; ++it3)
    {
        S8.push_back(it3->get()->GetCurrentState());
    }
    SufficientStatistics.push_back(S8);

    // Compute S9
    std::vector<double> S9;
    std::vector<std::shared_ptr< GaussianRandomVariable >>::iterator it4;
    for(it4 = m_AMatrixCoefficient.begin() ; it4 != m_AMatrixCoefficient.end() ; ++it4)
    {
        S9.push_back( it4->get()->GetCurrentState() );
    }
    SufficientStatistics.push_back(S9);

    return SufficientStatistics;
}

void
LongitudinalModel
::ComputeMaximizationStep(std::vector<std::vector<double>> S)
{
    // Update P0_mean, T0_Mean, V0_Mean
    m_P0->SetMean(S[4][0]);
    m_T0->SetMean(S[5][0]);
    m_V0->SetMean(S[6][0]);

    // Update delta(k)_mean
    int i = 0;
    for(std::vector<double>::iterator it = S[7].begin() ; it != S[7].end() ; ++it)
    {
        m_PropagationCoefficient[i]->SetMean(*it);
        i += 1;
    }

    // Update beta(k)_mean
    int j=0;
    for(std::vector<double>::iterator it = S[8].begin() ; it != S[8].end() ; ++it)
    {
        m_AMatrixCoefficient[j]->SetMean(*it);
        j += 1;
    }

    // Update Ksi_Variance
    double Value = 0;
    for(std::vector<double>::iterator it = S[2].begin() ; it != S[2].end() ; ++it)
    {
        Value += *it;
    }
    Value = Value/(double)S[2].size();
    for(std::vector<std::shared_ptr< GaussianRandomVariable >>::iterator it = m_PreAccelerationFactor.begin() ; it != m_PreAccelerationFactor.end() ; ++it)
    {
        it->get()->SetVariance(Value);
    }

    // Update tau_Variance
    Value = 0;
    for(std::vector<double>::iterator it = S[3].begin() ; it != S[3].end() ; ++it)
    {
        Value += *it;
    }
    Value = Value/(double)S[3].size();
    for(std::vector<std::shared_ptr< GaussianRandomVariable >>::iterator it = m_TimeShift.begin() ; it != m_TimeShift.end() ; ++it)
    {
        it->get()->SetVariance(Value);
    }

    // Update the uncertainty variance
    Value = 0;
    int K = 0;
    for(Data::iterator it = m_Data->begin() ; it != m_Data->end() ; ++it)
    {
        for(std::vector<std::pair<double, std::vector<double>>>::iterator it2 = it->begin() ; it2 != it->end() ; ++it2)
        {
            K += 1;
            for(std::vector<double>::iterator it3 = it2->second.begin() ; it3 != it2->second.end() ; ++it3)
            {
                // ACTUALLY THIS VALUE CAN BE COMPUTED ONLY ONCE
                Value += *it3;
            }
        }
    }

    for(std::vector<double>::iterator it = S[0].begin() ; it != S[0].end() ; ++it)
    {
        Value -= 2* *it;
    }

    for(std::vector<double>::iterator it = S[1].begin() ; it != S[1].end() ; ++it)
    {
        Value += *it;
    }

    Value = Value / (K * m_Data->size());
    Value = sqrt(Value);

    *m_UncertaintyVariance = Value;

}

