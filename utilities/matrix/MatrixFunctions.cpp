#include "MatrixFunctions.h"




LinearAlgebra<ScalarType>::VectorType
LinearCombination(LinearAlgebra<ScalarType>::VectorType Coefficients, std::vector<LinearAlgebra<ScalarType>::VectorType> Vectors)
{
    typedef LinearAlgebra<ScalarType>::VectorType VectorType;
    
    if(Coefficients.size() != Vectors.size())
    {
        std::cout << " The coefficients and vectors do not have the same size. How to do a linear combination?";
    }
    
    VectorType ReturnVector(Vectors[0].size(), 0.0);
    
    auto IterVect = Vectors.begin();
    auto IterCoef = Coefficients.begin();
    for(    ; IterCoef != Coefficients.end() && IterVect != Vectors.end(); ++IterCoef, ++IterVect)
    {
        ReturnVector += *IterCoef * *IterVect;
    }
    
    return ReturnVector;

}
