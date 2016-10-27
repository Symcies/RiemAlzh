#include "MatrixFunctions.h"


double
ComputeEuclideanScalarProduct(std::vector<double> A, std::vector<double> B)
{
	if(A.size() != B.size())
	{
        std::cout << "The matrixes do not have the same size. How to compute their euclidian scalar product?" << std::endl;
	}

	double ScalarProduct = 0;
    auto u = A.begin();
    auto v = B.begin();
	for( ; u != A.end() && v != B.end() ; ++u, ++v)
	{
		ScalarProduct += *u * *v;
	}

	return ScalarProduct;
}



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

double
NormOfVectorDifference(LinearAlgebra<ScalarType>::VectorType U, LinearAlgebra<ScalarType>::VectorType V)
{

    if(U.size() != U.size())
    {
        std::cout << " The vectors do not have the same size. How to do a difference?";
    }
    
    LinearAlgebra<ScalarType>::VectorType R = U - V;
    double Result = R.squared_magnitude();
         
    return Result;
}


// TO BE DELETED ABSOLUTELY !
LinearAlgebra<ScalarType>::VectorType
ConvertToVectorType(std::vector<double> VectorSTL)
{
    LinearAlgebra<ScalarType>::VectorType V(VectorSTL.size());
    
    int i = 0;
    for(auto it = VectorSTL.begin(); it != VectorSTL.end(); ++it, ++i)
    {
        V(i) = *it;
    }
    
    return V;
    
}