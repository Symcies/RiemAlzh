#pragma once

#include "ModelSettings.h"

#include "AbstractModel.h"
#include "UnivariateModel.h"
#include "MultivariateModel.h"

class ModelFactory {
 public:
  
  static std::shared_ptr<AbstractModel> NewModel(io::ModelSettings& settings);
};


