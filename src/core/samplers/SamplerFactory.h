#pragma once 

#include "AbstractSampler.h"
#include "BlockedGibbsSampler.h"
#include "MHwGSampler.h"

class SamplerFactory {
 public:
  
  static std::shared_ptr<AbstractSampler> NewSampler(std::string type);
};
