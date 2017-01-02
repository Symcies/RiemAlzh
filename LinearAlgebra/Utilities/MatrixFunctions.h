
#ifndef _MatrixFunctions_h
#define _MatrixFunctions_h

typedef double ScalarType;

#include <vector>
#include <iostream>
#include <cmath>
#include "../LinearAlgebra.h"

LinearAlgebra<ScalarType>::VectorType 
LinearCombination(LinearAlgebra<ScalarType>::VectorType Coefficients,  std::vector<LinearAlgebra<ScalarType>::VectorType> Vectors);

#endif //_MatrixFunctions_h
