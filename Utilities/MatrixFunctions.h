
#ifndef _MatrixFunctions_h
#define _MatrixFunctions_h

#include <vector>
#include <iostream>
#include <cmath>

double ComputeEuclideanScalarProduct(std::vector<double>, std::vector<double>);

std::vector<double> LinearCombination(std::vector<double> Coefficients,  std::vector<std::vector<double>> Vectors);

/// Computes transpose(U-V).(U-V)
double NormOfVectorDifference(std::vector<double> U, std::vector<double> V);

#endif //_MatrixFunctions_h
