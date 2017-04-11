#pragma once

typedef double ScalarType;

#include <iostream>
#include <fstream>
#include <string>

#include "InputsAssert.h"
#include "LinearAlgebra.h"
#include "Observations.h"
#include "RealDataSettings.h"

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
  static Observations ReadObservations(const RealDataSettings& ds);

  /// Open network propagation files
  static MatrixType OpenKernel(std::string file_path);

private:
  static VectorType ExtractObservation(std::ifstream& f_stream, unsigned int dimension);
  static IndividualObservations CreateIndividualObs(const RealDataSettings& ds, ReadData::VectorType& time_points,
    std::vector<ReadData::VectorType>& landmarks, std::vector<ReadData::VectorType>& cognitive_scores);
  static int CountLineNumber(std::ifstream& file);

  /// Copy constructor, private to prevent copy
  ReadData(const ReadData&);

  /// Assignment operator, private to prevent copy
  ReadData& operator=(const ReadData&);

};

} //end io namespace
