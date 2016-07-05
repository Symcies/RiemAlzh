#include "LongitudinalModel.h"

LongitudinalModel
::LongitudinalModel()
{
    // Initialize initial position
    double P0Mean = 1.0;
    double P0Variance = 1.0;
    m_P0 = std::make_shared< GaussianRandomVariable >(P0Mean, P0Variance);


    // Initialize initial time
    double T0Mean = 1.0;
    double T0Variance = 1.0;
    m_T0 = std::make_shared< GaussianRandomVariable >(T0Mean, T0Variance);


    // Initialize initial velocity
    double V0Mean = 1.0;
    double V0Variance = 1.0;
    m_V0 = std::make_shared< GaussianRandomVariable >(V0Mean, V0Variance);


    // Initialize initial propagation coefficient
    double DeltaVariance = 1.0;
    for(int i = 0; i<m_Manifold->GetNumberOfDimension() ; ++i)
    {
        double DeltaMean = (double)i;
        auto Delta = std::make_shared< GaussianRandomVariable >(DeltaMean, DeltaVariance);
        m_PropagationCoefficient.push_back(Delta);
    }


    // Initialize initial A Matrix coefficient
    double BetaVariance = 1.0;
    int ArraySize = (m_Manifold->GetNumberOfDimension() - 1)*m_Manifold->GetNumberOfIndependentComponents();
    for(int i = 0; i<ArraySize ; ++i)
    {
        double BetaMean = (double)i;
        auto Beta = std::make_shared< GaussianRandomVariable >(BetaMean, BetaVariance);
        m_AMatrixCoefficient.push_back(Beta);
    }


    // Initialize the first Orthonormal basis (B1, ..., B(n-1)Ns)
    InitializeOrthonormalBasis();

}

void
LongitudinalModel
::UpdateSubjectSpecificParameters(Data *D)
{

    // Set the number of subjects
    m_NumberOfSubjects = (int)D->size();

    // Initiate the data
    m_Data = D;

    // Initialize Pre acceleration factor
    double KsiVariance = 1.0;
    double KsiMean = 0.0;  // Always equal to 0
    for(int i=0; i<m_NumberOfSubjects ; ++i)
    {
        auto Ksi = std::make_shared<GaussianRandomVariable>(KsiMean, KsiVariance);
        m_PreAccelerationFactor.push_back(Ksi);
    }


    // Initialize time shifts
    double TauVariance = 1.0;
    double TauMean = 0.0; // always equal to 0
    for(int i=0; i<m_NumberOfSubjects ; ++i)
    {
        auto Tau = std::make_shared< GaussianRandomVariable>(TauMean, TauVariance);
        m_TimeShift.push_back(Tau);
    }


    // Initialize the uncertainty
    m_UncertaintyVariance = std::make_shared<double>(1.0);


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
            std::vector<double> Biomarkers = it2->second;

            Likelihood += ObservationDifferenceNorm(ParallelCurve, Biomarkers);
        }

        SubjectNumber +=1;
    }

    return exp( - 1.0/(2.0 * *m_UncertaintyVariance) * Likelihood);

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
    double P0CandidateVariance = 1.0;
    auto P0Candidate = std::make_shared<GaussianRandomVariable>(m_P0->GetCurrentState(), P0CandidateVariance);
    RandomVariableToSample P0 (m_P0, P0Candidate);
    RVToSample.push_back(P0);


    // Sample T0 and its proposition distribution
    double T0CandidateVariance = 1.0;
    auto T0Candidate = std::make_shared<GaussianRandomVariable>(m_T0->GetCurrentState(), T0CandidateVariance);
    RandomVariableToSample T0 (m_T0, T0Candidate);
    RVToSample.push_back(T0);


    // Sample V0 and its proposition distribution
    double V0CandidateVariance = 1.0;
    auto V0Candidate = std::make_shared<GaussianRandomVariable>(m_V0->GetCurrentState(), V0CandidateVariance);
    RandomVariableToSample V0 (m_V0, V0Candidate);
    RVToSample.push_back(V0);


    // Sample propagation coefficients and their proposition distribution
    double DeltaCandidateVariance = 1.0;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_PropagationCoefficient.begin() ; it != m_PropagationCoefficient.end() ; ++it)
    {
        auto DeltaCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), DeltaCandidateVariance);
        RandomVariableToSample Delta (*it, DeltaCandidate);
        RVToSample.push_back(Delta);
    }


    // Sample the A Matrix coefficient
    double BetaCandidateVariance = 1.0;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_AMatrixCoefficient.begin() ; it != m_AMatrixCoefficient.end() ; ++it)
    {
        auto BetaCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), BetaCandidateVariance);
        RandomVariableToSample Beta (*it, BetaCandidate);
        RVToSample.push_back(Beta);
    }


    // Sample the pre-acceleration factors
    double KsiCandidateVariance = 1.0;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_PreAccelerationFactor.begin(); it != m_PreAccelerationFactor.end() ; ++it)
    {
        auto KsiCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), KsiCandidateVariance);
        RandomVariableToSample Ksi (*it, KsiCandidate);
        RVToSample.push_back(Ksi);
    }


    // Sample the time shift
    double TauCandidateVariance = 1.0;
    for(std::vector<std::shared_ptr<GaussianRandomVariable>>::iterator it = m_TimeShift.begin() ; it != m_TimeShift.end() ; ++it)
    {
        auto TauCandidate = std::make_shared<GaussianRandomVariable>((*it)->GetCurrentState(), TauCandidateVariance);
        RandomVariableToSample Tau (*it, TauCandidate);
        RVToSample.push_back(Tau);
    }


    // Sample the Space shift coefficient
    double SCandidateVariance = 1.0;
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
// Protected attribute(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


void
LongitudinalModel
::SetAlgorithmParametersBYVALUE()
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // In the parameters to return:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    m_AlgorithmParameters.clear();

    // Add the initial position mean
    m_AlgorithmParameters.push_back(m_P0->GetMean());

    // Add the initial time mean
    m_AlgorithmParameters.push_back(m_T0->GetMean());

    // Add the initial velocity mean
    m_AlgorithmParameters.push_back(m_V0->GetMean());

    // Add the pre acceleration factor variance
    m_AlgorithmParameters.push_back(m_PreAccelerationFactor[0]->GetMean());

    // Add the time shift variance
    m_AlgorithmParameters.push_back(m_TimeShift[0]->GetMean());

    // Add the uncertainty variance
    m_AlgorithmParameters.push_back(*m_UncertaintyVariance);
}

void
LongitudinalModel
::InitializeOrthonormalBasis()
{
    if(m_V0->GetCurrentState() == 0.0)
    {
        throw std::invalid_argument(" V0 equals zero : impossible to calculate the orthogonal basis");
    }

    std::vector<std::vector<double>> FirstBasis;

    for(int i=1; i<m_Manifold->GetNumberOfDimension() ; ++i)
    {
        std::vector<double> vector;
        for(int j = 0; j<m_Manifold->GetNumberOfDimension() ; ++j)
        {
            if(i==j)
            {
                vector.push_back(1.0);
            }
            else
            {
                vector.push_back(0.0);
            }
        }
        FirstBasis.push_back(vector);
    }

    ////////////////////////////////////////////////////////////////////////////
    ///// JUST NEED TO CHECK THAT (e2, e3, ..., eN) is a basis of Orthogonal(V0)
    ///// <=> (v0, e2, ..., eN) is a basis
    ////////////////////////////////////////////////////////////////////////////

    m_OrthonormalBasis = FirstBasis;
}

void
LongitudinalModel
::ComputeOrthonormalBasis()
{
    if(m_V0->GetCurrentState() == 0.0)
    {
        throw std::invalid_argument(" V0 equals zero : impossible to calculate the orthogonal basis");
    }

    /// Modified Gram-Schmidt Process
    /// http://cavern.uark.edu/~arnold/4353/CGSMGS.pdf
    /// http://ocw.mit.edu/courses/mathematics/18-335j-introduction-to-numerical-methods-fall-2010/lecture-notes/MIT18_335JF10_lec10a_hand.pdf

    std::vector<std::vector<double>> NewOrthogonalBasis;

    /// Copy the existing orthonormal basis and adding v0 to it at first place
    /// So that it is a basis of the N-dimensional manifold
    std::vector<std::vector<double>> OrthogonalBasis = m_OrthonormalBasis;
    std::vector<double> VectorV0;

    for(int i = 0; i<m_Manifold->GetNumberOfDimension(); ++i) VectorV0.push_back(m_V0->GetCurrentState());

    OrthogonalBasis.insert(OrthogonalBasis.begin(), VectorV0);

    for(int k=0; k<m_Manifold->GetNumberOfDimension(); ++k)
    {
        std::vector<double> w = OrthogonalBasis[k];
        double r = 0;
        for(int j=0; k-1; ++j)
        {
            for(int i=0; i<m_Manifold->GetNumberOfDimension(); ++i)
            {
                r += NewOrthogonalBasis[j][i]*w[i];
            }

            for(int i=0; i<m_Manifold->GetNumberOfDimension(); ++i)
            {
                w[i] -= r*NewOrthogonalBasis[j][i];
            }
        }

        std::vector<double> NewBasisVector;
        double norm = 0;
        for(int j = 0; j<m_Manifold->GetNumberOfDimension(); ++j)
        {
            norm += w[j]*w[j];
        }
        norm = sqrt(norm);
        for(int j = 0; j<m_Manifold->GetNumberOfDimension(); ++j)
        {
            NewBasisVector.push_back(w[j]/norm);
        }
        NewOrthogonalBasis.push_back(NewBasisVector);
    }

    /// Erase the V0 vector
    NewOrthogonalBasis.erase(NewOrthogonalBasis.begin());
    m_OrthonormalBasis = NewOrthogonalBasis;



}

void
LongitudinalModel
::ComputeAMatrix()
{

    m_AMatrix.clear();

    int Dimension = m_Manifold->GetNumberOfDimension();

    for(int i = 0; i<Dimension; ++i)
    {
        for(int j = 0; j<m_Manifold->GetNumberOfIndependentComponents() ; ++j)
        {
            double val = 0;
            for(int k = 0; k<Dimension-1; ++k)
            {
                val += m_AMatrixCoefficient[k+j*Dimension]->GetCurrentState() * m_OrthonormalBasis[k][i];
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

    for(int i = 0; i<m_NumberOfSubjects ; ++i)
    {
        std::vector< std::shared_ptr<LaplaceRandomVariable> > SpaceShiftCoefficient = m_SpaceShiftCoefficient[i];
        std::vector<double> IndividualSpaceShift;
        for(int j = 0; j<m_Manifold->GetNumberOfDimension(); ++j)
        {
            double val = 0.0;
            for(int k = 0; k<m_Manifold->GetNumberOfIndependentComponents(); ++k)
            {
                val += m_AMatrix[i + j*m_Manifold->GetNumberOfDimension()] * SpaceShiftCoefficient[j]->GetCurrentState();
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

