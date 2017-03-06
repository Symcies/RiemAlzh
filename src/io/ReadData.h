#ifndef ReadData_h
#define ReadData_h

typedef double ScalarType;

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <ios>
#include <src/observations/Observations.h>

#include "Observations.h"
#include "DataSettings.h"
#include "LinearAlgebra.h"

namespace io {

class ReadData {
public:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// typedef :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    typedef typename LinearAlgebra<ScalarType>::MatrixType MatrixType;
    typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;

    typedef std::vector<std::vector<std::pair<LinearAlgebra<ScalarType>::VectorType, double> > > Data;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Constructor(s) / Destructor :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Encapsulation method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Other method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////
  
    /// Open observations
    static Observations ReadObservations(DataSettings& DS);
    
    /// Open network propagation files
    static MatrixType OpenKernel(std::string FilePath);


protected:
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Attribute(s)
    ////////////////////////////////////////////////////////////////////////////////////////////////////


};

}

#endif //ReadData_h
