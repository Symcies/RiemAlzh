#pragma once

typedef double ScalarType;

#include <iostream>
#include <fstream>
#include <string>

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
    ReadData();
    ~ReadData();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Method(s) :
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Open observations
    static Observations ReadObservations(DataSettings& ds);

    /// Open network propagation files
    static MatrixType OpenKernel(std::string file_path);

private:
    /// Copy constructor, private to prevent copy
    ReadData(const ReadData&);

    /// Assignment operator, private to prevent copy
    ReadData& operator=(const ReadData&);

};

} //end io namespace
