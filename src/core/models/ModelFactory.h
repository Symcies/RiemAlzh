#pragma once

#include "ModelSettings.h"

#include "AbstractModel.h"
#include "UnivariateModel.h"
#include "MultivariateModel.h"
#include "FastNetworkModel.h"
#include "NetworkModel.h"
#include "MeshworkModel.h"
#include "GaussianModel.h"
#include "GaussianMixtureModel.h"

class ModelFactory {
 public:
  
  static std::shared_ptr<AbstractModel> NewModel(io::ModelSettings& settings);
};


