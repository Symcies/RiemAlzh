#ifndef _AlgorithmSettings_h
#define _AlgorithmSettings_h

#include <string>
#include <iostream>

#include "tinyxml2.h"

namespace io {

class AlgorithmSettings {

public:

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    AlgorithmSettings(char *XMLFile);

    ~AlgorithmSettings();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get the maximum number of iterations
    unsigned int GetMaximumNumberOfIterations() { return m_MaximumNumberOfIterations; }

    /// Get the number of burn-in iterations
    unsigned int GetNumberOfBurnInIterations() { return m_NumberOfBurnInIterations; }

    /// Get the number of iterations to wait till newt output display
    unsigned int GetCounterToDisplayOutputs() { return m_CounterToDisplayOutputs; }

    /// Get the number of iteration to wait till next data saving
    unsigned int GetCounterToSaveData() { return m_CounterToSaveData; }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Iteration attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Maximum number of iteration of the Algorithm
    unsigned int m_MaximumNumberOfIterations;

    /// Number of Burn-In Iterations
    unsigned int m_NumberOfBurnInIterations;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Display attribute(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Number of iterations to wait till newt output display
    unsigned int m_CounterToDisplayOutputs;

    /// Number of iteration to wait till next data saving
    unsigned int m_CounterToSaveData;


    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Methods(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /// Convert the max number of iteration from the xml file

};


}
#endif //_AlgorithmSettings_h

