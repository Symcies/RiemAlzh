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
    std::uniform_real_distribution<double> Distribution(-1.0/2.0, 1.0/2.0);

    double UniDraw = Distribution(Generator);
    return  m_Location - m_Scale*copysign( 1.0, UniDraw) * log( 1 - 2.0*fabs(UniDraw) );

}


double
LaplaceRandomVariable
::Likelihood(double X)
{
    return exp( - fabs(X - m_Location) / m_Scale) / (2.0*m_Scale) ;
}


double 
LaplaceRandomVariable
::LogLikelihood(double X) 
{
    double LogLikelihood = -log(2.0*m_Scale);
    LogLikelihood += -fabs(X - m_Location) / m_Scale;
    return LogLikelihood;
}