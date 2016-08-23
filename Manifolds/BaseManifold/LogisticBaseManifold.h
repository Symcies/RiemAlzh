#ifndef _LogisticBaseManifold_h
#define _LogisticBaseManifold_h

#include <math.h>

#include "AbstractBaseManifold.h"

class LogisticBaseManifold : public AbstractBaseManifold {
public:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    LogisticBaseManifold();
    ~LogisticBaseManifold();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compute the geodesic
    virtual double ComputeGeodesic(double P0, double T0, double V0, double TimePoint);

    /// Compute the geodesic derivative
    virtual double ComputeGeodesicDerivative(double P0, double T0, double V0, double TimePoint);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
};


#endif //_LogisticBaseManifold_h
