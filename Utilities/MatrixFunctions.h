
#ifndef _MatrixFunctions_h
#define _MatrixFunctions_h

typedef double ScalarType;

#include <vector>
#include <iostream>
#include <cmath>
#include "../LinearAlgebra/LinearAlgebra.h"

double ComputeEuclideanScalarProduct(std::vector<double>, std::vector<double>);

LinearAlgebra<ScalarType>::VectorType LinearCombination(LinearAlgebra<ScalarType>::VectorType Coefficients,  std::vector<LinearAlgebra<ScalarType>::VectorType> Vectors);

/// Computes transpose(U-V).(U-V)
double NormOfVectorDifference(LinearAlgebra<ScalarType>::VectorType U, LinearAlgebra<ScalarType>::VectorType V);


LinearAlgebra<ScalarType>::VectorType
ConvertToVectorType(std::vector<double> VectorSTL);

#endif //_MatrixFunctions_h
