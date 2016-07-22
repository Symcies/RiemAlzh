#include "LaplaceRandomVariable.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

LaplaceRandomVariable
::LaplaceRandomVariable(double Location, double Scale)
{
    m_Location = Location;
    m_Scale = Scale;
}

LaplaceRandomVariable
::~LaplaceRandomVariable()
{ }


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

double
LaplaceRandomVariable
::Sample()
{
    std::random_device RD;
    std::mt19937_64 Generator(RD());
    std::uniform_real_distribution<double> Distribution(1/2, 1/2);

    double UniDraw = Distribution(Generator);
    return m_Location - m_Scale*copysign( 1.0, UniDraw) * log( 1 - 2*fabs(UniDraw) );

}


double
LaplaceRandomVariable
::Likelihood(double X)
{
    return exp( - fabs(X - m_Location) / m_Scale) / (2*m_Scale) ;
}