#pragma once

typedef double ScalarType;

#include <string>
#include <unordered_map>
#include <tuple>
#include <iostream>

#include "LinearAlgebra.h"
#include "tinyxml2.h"

#include "InputsAssert.h"

namespace io {

class ModelSettings {

public:
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
  typedef std::unordered_map<std::string, std::pair<std::vector<double>, ScalarType>> InitialRVParameters;
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  ModelSettings(std::string xml_file);

  ~ModelSettings();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  std::string GetType() { return type_; }

  unsigned int GetIndependentSourcesNumber() const { return independent_sources_nb_; }

  std::string GetInvertKernelPath() const { return invert_kernel_matrix_path_; }

  std::string GetInterpolationKernelPath() const { return interpolation_matrix_path_; }
  
  InitialRVParameters GetInitialRandomVariables() const { return init_random_variables_; }

  InitialRVParameters GetSecondRandomVariables() const { return second_random_variables_; }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

private:

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// Type of the model
  std::string type_;

  /// Name of the output file
  std::string output_file_name_;
  
  /// Number of sources
  unsigned int independent_sources_nb_;

  /// Path to the kernel matrix
  std::string invert_kernel_matrix_path_;

  /// Path to the interpolation matrix
  std::string interpolation_matrix_path_;
  
  /// Initial model random variables
  InitialRVParameters init_random_variables_;
  InitialRVParameters second_random_variables_;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Model specific methods(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// Path to kernel Kxd

  /// Path to kernel invKd

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// Methods(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /// Load the model variables
  void LoadInitialRandomVariables(const tinyxml2::XMLElement *settings);
  
  /// Load the random varianle parameters
  std::vector<double> LoadRVParameters(const tinyxml2::XMLElement* parameters);
  
  /// Chooses the right model
  void LoadModel(const tinyxml2::XMLElement *settings);

  /// Load the fast network model
  void LoadFastNetwork(const tinyxml2::XMLElement *settings);

  /// Load the meshwork model
  void LoadMeshworkModel(const tinyxml2::XMLElement *settings);

  /// Load the network model
  void LoadNetworkModel(const tinyxml2::XMLElement *settings);

  /// Load the multivariate
  void LoadMultivariate(const tinyxml2::XMLElement *settings);

  /// Load the univariate
  void LoadUnivariate(const tinyxml2::XMLElement *settings);
  
  /// Load the gaussian model
  void LoadGaussian(const tinyxml2::XMLElement *settings);

  ///Print the type of the model
  void PrintModelInfo();

  /// Copy constructor, private to prevent copy
  ModelSettings(const ModelSettings&);

  /// Assignment operator, private to prevent copy
  ModelSettings& operator=(const ModelSettings&);

};

}
