#include "MatrixFunctions.h"


double
ComputeEuclideanScalarProduct(std::vector<double> A, std::vector<double> B)
{
	if(A.size() != B.size())
	{
        std::cout << "The matrixes do not have the same size. How to compute their euclidian scalar product?" << std::endl;
	}

	double ScalarProduct = 0;
	for(std::pair<std::vector<double>::iterator, std::vector<double>::iterator> i(A.begin(), B.begin()) ;
		i.first != A.end() && i.second != B.end() ;
		++i.first, ++ i.second)
	{
		ScalarProduct += *i.first * *i.second;
	}

	return ScalarProduct;
}



std::vector<double>
LinearCombination(std::vector<double> Coefficients, std::vector<std::vector<double>> Vectors)
{
    if(Coefficients.size() != Vectors.size())
    {
        std::cout << " The coefficients and vectors do not have the same size. How to do a linear combination?";
    }
    std::vector<double> ReturnVector(Vectors[0].size(), 0.0);

    typedef std::vector<double>::iterator DoubleIter;
    for(int i = 0; i < Vectors.size(); ++i)
    {
        double Coeff = Coefficients[i];
        std::vector<double> Vector = Vectors[i];

        for(std::pair<DoubleIter, DoubleIter> i(Vector.begin(), ReturnVector.begin());
            i.first != Vector.end() && i.second != ReturnVector.end();
            ++i.first, ++i.second)
        {
            *i.second += Coeff* *i.first;
        }
    }

    return ReturnVector;

}

double
NormOfVectorDifference(std::vector<double> U, std::vector<double> V)
{

    if(U.size() != U.size())
    {
        std::cout << " The vectors do not have the same size. How to do a difference?";
    }

    typedef std::vector<double>::iterator DoubleIter;
    double Diff = 0;
    for(std::pair<DoubleIter, DoubleIter> i(U.begin(), V.begin());
        i.first != U.end() && i.second != V.end();
        ++i.first, ++i.second)
    {
        Diff += (i.first-i.second) * (i.first-i.second);
    }

    return Diff;
}