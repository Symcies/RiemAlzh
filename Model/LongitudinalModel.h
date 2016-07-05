#ifndef _LongitudinalModel_h
#define _LongitudinalModel_h

#include <stdexcept>            // throw
#include <memory>               // unique_ptr and smart_ptr
#include <cmath>

#include <AbstractManifold.h>
#include <AbstractRandomVariable.h>
#include <GaussianRandomVariable.h>
#include <MultivariateLogisticManifold.h>
#include <LaplaceRandomVariable.h>

typedef std::vector<std::vector<std::pair<double, std::vector<double>>>> Data;
typedef std::pair<std::shared_ptr<AbstractRandomVariable>, std::shared_ptr<AbstractRandomVariable> > RandomVariableToSample;
typedef std::vector<double> AlgorithmParameters;

class LongitudinalModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    LongitudinalModel();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getter(s) and Setter(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Get the parameters for the algorithm
    AlgorithmParameters GetAlgorithmParametersBYVALUE() { return m_AlgorithmParameters; }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Update the parameters that are subjects specific : can be done only when the algorithm has initialize its longitudinal model
    void UpdateSubjectSpecificParameters(Data *D);

    // Compute the likelihood
    double ComputeLikelihood();

    // Compute the parallel transport of subject i at observation j
    std::vector<double> ComputeParallelTransport(int i, double ObservationTimePoint);

    // Get the random variable to sample
    std::vector<RandomVariableToSample> GetRandomVariableToSample();


    // Initialize the sufficient statistics
    std::vector<std::vector< double >> InitializeSufficientStochasticStatistics();

    // Compute the sufficient statistics
    std::vector<std::vector<double>> ComputeSufficientStatistics();

    // Update the parameters thanks to the MCMC SAEM maximization step
    void ComputeMaximizationStep(std::vector<std::vector<double>> StochasticSufficientStatistics);




protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Number of subjects
    int  m_NumberOfSubjects;

    // Subject related data
    Data *m_Data;

    // Initial Position p0
    std::shared_ptr< GaussianRandomVariable > m_P0;

    // Initial Time t0
    std::shared_ptr< GaussianRandomVariable > m_T0;

    // Initial Velocity v0
    std::shared_ptr< GaussianRandomVariable > m_V0;

    // Propagation coefficient (delta)k among the biomarkers
    std::vector<std::shared_ptr< GaussianRandomVariable >> m_PropagationCoefficient;

    // Coefficient (beta)k of the A matrix
    std::vector<std::shared_ptr< GaussianRandomVariable >> m_AMatrixCoefficient;

    // Pre acceleration coefficient (ksi)i of the subjects
    std::vector<std::shared_ptr< GaussianRandomVariable >> m_PreAccelerationFactor;

    // Time shift (tau)i of the subjects
    std::vector<std::shared_ptr< GaussianRandomVariable >> m_TimeShift;

    // Uncertainties
    std::shared_ptr<double> m_UncertaintyVariance;

    // Vectors of the space shift coefficient of the subjects
    std::vector<std::vector< std::shared_ptr< LaplaceRandomVariable > >> m_SpaceShiftCoefficient;

    // Space shift vector of each subject
    std::vector<std::vector< double >> m_SpaceShift;

    // Manifold used for the simulation
    std::shared_ptr< MultivariateLogisticManifold > m_Manifold;

    // Orthonormal Basis ( B1, ..., B(N-1)Ns ) used to construct A
    std::vector<std::vector<double>> m_OrthonormalBasis;

    // A Matrix to construct the time shifts: N rows, Ns column
    std::vector<double> m_AMatrix;

    // Parameters to send to the algorithm class
    // FIGURE OUT IF IT WOULD BE BETTER TO PASS AS VARIABLES NOT TO UPDATE IT AT EACH STEP
    AlgorithmParameters m_AlgorithmParameters;



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set the parameters to send back to the algorithm
    void SetAlgorithmParametersBYVALUE();

    // Initialize the first orthonormal basis
    void InitializeOrthonormalBasis();

    // Calculate the orthonomal basis thanks to Gram Schmidt
    void ComputeOrthonormalBasis();

    // Compute the A Matrix
    void ComputeAMatrix();

    // Compute the space shifts
    void ComputeSpaceShifts();

    // Compute the square of the difference between the observation yij and its corresponding parallel curve eta(...)
    double ObservationDifferenceNorm(std::vector<double> Observation, std::vector<double> ParallelCurve);

};


#endif //_LongitudinalModel_h
