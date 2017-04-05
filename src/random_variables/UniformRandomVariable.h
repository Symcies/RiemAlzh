#pragma once

#include "AbstractRandomVariable.h"

class UniformRandomVariable : public AbstractRandomVariable {
 public:
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  UniformRandomVariable(ScalarType min, ScalarType max);
  virtual ~UniformRandomVariable();

  UniformRandomVariable(const UniformRandomVariable&);
  virtual UniformRandomVariable& operator=(const UniformRandomVariable&);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Getter(s) and Setter(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////

  virtual ScalarType GetParameter(std::string param_name) const;
  virtual ScalarType GetParameter(int param_key) const;

  inline double GetMin() const { return min_; }

  inline double GetMax() const { return max_; }

  inline void SetMin(double min) { min_ = min; };

  inline void SetMax(double max) {max_ = max; };

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Draw a sample
  virtual double Sample();

  /// Compute the likelihood given a current state
  virtual double Likelihood(double x);

  /// Compute the loglikelihood given a current state
  virtual double LogLikelihood(double x);

  /// Update the random variable parameters
  virtual void Update(StringScalarHash params);
  virtual void Update(IntScalarHash    params);

private:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// min of the uniform distribution
  double min_;

  /// max of the uniform distribution
  double max_;

};
