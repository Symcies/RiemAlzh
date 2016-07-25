#ifndef _AbstractModel_h
#define _AbstractModel_h

#include <memory>

#include "../RandomVariables/AbstractRandomVariable.h"
#include "../RandomVariables/LaplaceRandomVariable.h"
#include "../Manifolds/AbstractManifold.h"

class AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    typedef std::vector< std::vector< std::pair< std::vector<double>, unsigned int> > > Data;
    typedef std::vector< std::pair< std::vector<double>, unsigned int> > IndividualData;
    typedef std::map< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariableMap;
    typedef std::pair< std::string, std::shared_ptr< AbstractRandomVariable >> RandomVariable;
    typedef std::map<std::string, double> Realizations;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void SetManifold(std::shared_ptr< AbstractManifold> M) { m_Manifold = M; }

    inline RandomVariableMap GetPopulationrandomVariables() { return m_PopulationRandomVariables; }

    inline RandomVariableMap GetIndividualRandomVariables() { return m_IndividualRandomVariables; }

    inline std::vector< std::vector< double >> GetSufficientStatistics() { return m_SufficientStatistics; }

    inline void SetSufficientStatistics(std::vector< std::vector< double >> NewStatistics) { m_SufficientStatistics = NewStatistics; }


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the random variables : Population-wide and subject-specifid
    virtual void InitializeRandomVariables() = 0;

    /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatistics GetSufficientStatistics(const Realizations& R, const Data& D) = 0;

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const std::vector< std::vector< double >>& SufficientStatistics, const Data& D) = 0;
   
protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Riemanian manifold
    std::shared_ptr< AbstractManifold > m_Manifold;

    /// Random variables shared among the population
    RandomVariableMap m_PopulationRandomVariables;

    /// Random variables that are subject-specific
    RandomVariableMap m_IndividualRandomVariables;


};


#endif //_AbstractModel_h
