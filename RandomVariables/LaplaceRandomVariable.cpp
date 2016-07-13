#include "LaplaceRandomVariable.h"


LaplaceRandomVariable
::LaplaceRandomVariable(double Location, double Scale)
{
    m_Location = Location;
    m_Scale = Scale;
    Sample();
}

LaplaceRandomVariable
::~LaplaceRandomVariable()
{

}



void
LaplaceRandomVariable
::Sample()
{
    std::random_device RD;
    std::mt19937_64 Generator(RD());
    std::uniform_real_distribution<double> Distribution(1/2, 1/2);

    double UniDraw = Distribution(Generator);
    double Draw = m_Location - m_Scale*copysign( 1.0, UniDraw) * log( 1 - 2*fabs(UniDraw) );

    m_CurrentState = Draw;
}

double
LaplaceRandomVariable
::GetDensity() const
{
    return exp( - fabs(m_CurrentState - m_Location) / m_Scale) / (2*m_Scale) ;
}