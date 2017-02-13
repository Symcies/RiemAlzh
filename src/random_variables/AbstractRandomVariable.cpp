#include "AbstractRandomVariable.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


AbstractRandomVariable::VectorType
AbstractRandomVariable
::Samples(unsigned int NumberOfSamples) 
{
    VectorType Samples(NumberOfSamples);
    ScalarType * s = Samples.memptr();
    
    for(size_t i = 0; i < NumberOfSamples; ++i)
        s[i] = Sample();
    
    return Samples;
}