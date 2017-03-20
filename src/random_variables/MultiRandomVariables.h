#pragma once

#include <memory>
#include <unordered_map>

#include "Realizations.h"
#include "AbstractRandomVariable.h"
#include "GaussianRandomVariable.h"

typedef double ScalarType;


class MultiRandomVariables {

public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef typename LinearAlgebra<ScalarType>::VectorType VectorType;
  typedef typename std::unordered_map<int, std::shared_ptr<AbstractRandomVariable>> IntRandomVariableHash;
  typedef typename std::unordered_map<std::string, int> StringIntHash;
  typedef typename std::unordered_map<std::string, ScalarType > StringScalarHash;
  typedef typename std::unordered_map<int, ScalarType > IntScalarHash;
  typedef typename std::unordered_map<int, int > IntIntHash;
  typedef typename std::unordered_map<int, std::string> IntStringHash;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  MultiRandomVariables();
  ~MultiRandomVariables();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  std::unique_ptr<AbstractRandomVariable> GetRandomVariable(std::string name) const;

  std::unique_ptr<AbstractRandomVariable> GetRandomVariable(int key) const;

  void Clear();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////



  void AddRandomVariable(std::string name, std::string type, const std::vector<double>& params);

  /// Update a random variable based on its name
  void UpdateRandomVariable(std::string name, IntScalarHash params);
  void UpdateRandomVariable(std::string name, StringScalarHash params);

  /// Update a random variable based on its key
  void UpdateRandomVariable(int key, IntScalarHash params);
  void UpdateRandomVariable(int key, StringScalarHash params);

  /// Simulate realizations
  Realizations SimulateRealizations(StringIntHash num_of_real_per_rand_var);

private:



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Hash table of the random variables
  IntRandomVariableHash rand_var_;

  /// Convert the names into the int key
  StringIntHash string_to_int_key_;

  /// Convert the int into their names
  IntStringHash int_to_string_key_;

  /// Convert the variable key into their type
  IntStringHash rand_var_string_type_key_;
  IntIntHash    rand_var_int_type_key_;

  /// Int key counter
  int key_count_ = 0;

};
