#ifndef _Algorithm_h
#define _Algorithm_h


#include "../Models/AbstractModel.h"
#include "../Samplers/AbstractSampler.h"

class Algorithm {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    typedef std::vector< std::vector< std::pair< std::vector<double>, double> > > Data;
    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, std::vector<double>> Realizations;
    typedef std::vector< std::vector< double >> SufficientStatisticsVector;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    Algorithm();
    ~Algorithm();



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void SetModel(const std::shared_ptr<AbstractModel>& M ) { m_Model = M; }

    inline void SetSampler(const std::shared_ptr<AbstractSampler>& S) { m_Sampler = S; }


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the MCMC SAEM algorithm
    void ComputeMCMCSAEM(const std::shared_ptr<Data>& D);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) of the MCMC SAEM:
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initialize the stochastic approximation
    void InitializeStochasticSufficientStatistics(const SufficientStatisticsVector& S);

    // Initialize the realization of the model (and related manifold) random variables
    void InitializeRealization(unsigned int NbIndividuals);

    // Initialize the manifold parameters
    void InitializeModelParameters(std::shared_ptr<Realizations>& R);

    // Compute the simulation step : Gibbs Sampling
    void ComputeSimulationStep(const std::shared_ptr<Data>& D);

    // Compute the stochastic coefficient 
    void ComputeStochasticApproximation(double iteration, std::vector< std::vector< double >> SufficientStatistics);

    // Compute the decreasing step size of the approximation step
    double DecreasingStepSize(double Iteration, double NoMemoryTime);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Abstract Model
    std::shared_ptr<AbstractModel> m_Model;

    /// Realisation of the random variables of the model
    std::shared_ptr<Realizations> m_Realization;

    /// Abstract Sampler - for the MCMC SAEM 
    std::shared_ptr<AbstractSampler> m_Sampler;

    /// Stochastic Sufficient Statistics used in the stochastic approximation step
    std::vector< std::vector< double >> m_StochasticSufficientStatistics;


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) of eventual other algorithm:
    ////////////////////////////////////////////////////////////////////////////////////////////////////


};


#endif //_Algorithm_h
