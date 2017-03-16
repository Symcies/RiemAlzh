#include "AbstractRandomVariable.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


AbstractRandomVariable::VectorType AbstractRandomVariable::Samples(unsigned int samples_num)
{
    VectorType samples(samples_num);
    ScalarType * s = samples.memptr();

    for(size_t i = 0; i < samples_num; ++i)
        s[i] = Sample();

    return samples;
}
