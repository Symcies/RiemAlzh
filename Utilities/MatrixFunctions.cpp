#include "MatrixFunctions.h"


double
ComputeEuclidianScalarProduct(std::vector<double> A, std::vector<double> B)
{
	if(A.size() != B.size())
	{
		cout << "The matrixes does not have the same size. How to compute their euclidian scalar product?"
	}

	double ScalarProduct = 0;
	for(std::pair<std::vector<double>, std::vector<double>> i(A.begin(), B.begin()) ;
		i.first != A.end() && i.second() != B.end() ;
		++i.first, ++ i.second)
	{
		ScalarProduct += i.first * i.second;
	}

	return ScalarProduct;
}
