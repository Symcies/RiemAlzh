#ifndef _LongitudinalModel_h
#define _LongitudinalModel_h


#include "AbstractModel.h"
#include "../Utilities/MatrixFunctions.h"

class LongitudinalModel : public AbstractModel {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    LongitudinalModel(unsigned int NbIndependentComponents);
    ~LongitudinalModel();


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize the random variables : Population-wide and subject-specific
    virtual void InitializeRandomVariables();

     /// Update the sufficient statistics according to the model variables / parameters 
    virtual SufficientStatisticsVector GetSufficientStatistics(const Realizations& R, const Data& D);

    /// Update the fixed effects thanks to the approximation step of the algorithm
    virtual void UpdateRandomVariables(const SufficientStatisticsVector& SufficientStatistics, const Data& D);

    /// Compute the likelihood of the model
    virtual double ComputeLikelihood(const Realizations& R, const Data& D);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Initialize Population random variables
    void InitializePopulationRandomVariables();

    /// Initialize Individual random variables
    void InitializeIndividualRandomVariables();

    /// Compute Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns>
    void ComputeOrthonormalBasis( const Realizations& R); // TODO : Use a library to do it faster

    /// Compute the A Matrix used to get the space shifts
    void ComputeAMatrix( const Realizations& R); // TODO : Use a library to do it faster

    // Compute the space shifts
    void ComputeSpaceShifts(const Realizations& R, const Data& D); // TODO : Use a library to do it faster


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Number of independent components
    unsigned int m_NbIndependentComponents;

    /// Noise model
    std::shared_ptr< GaussianRandomVariable > m_Noise;

    /// Orthonormal Basis vec<B1, ..., B(N-1)> where Bi is vec<Ns> (Basis orthogonal to gamma0_deriv(T0)
    std::vector< std::vector< double >> m_OrthogonalBasis; // TODO : Initialize somewhere

    /// A Matrix vec<A1, ..., A(N)> where Ai is vec<Ns> (Ai is a column)
    std::vector< std::vector< double >> m_AMatrix; // TODO : Initialize somewhere

    /// Space shifts w(i) of the model
    std::map< std::string, std::vector< double >> m_SpaceShifts; // TODO : Initialize somewhere


};


#endif //_LongitudinalModel_h
