#include "SamplerFactory.h"


std::shared_ptr<AbstractSampler> SamplerFactory::NewSampler(std::string type) {
  
  if(type == "MHwG")        return std::make_shared<MHwGSampler>();
  else if (type == "Block") return std::make_shared<BlockedGibbsSampler>();
  
}