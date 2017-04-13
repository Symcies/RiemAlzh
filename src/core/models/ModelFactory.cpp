#include "ModelFactory.h"


std::shared_ptr<AbstractModel> ModelFactory::NewModel(io::ModelSettings &settings) {
  
  if (settings.GetType() == "Univariate") 
    return std::make_shared<UnivariateModel>(settings);
  else if (settings.GetType() == "Multivariate")
    return std::make_shared<MultivariateModel>(settings);
  else 
    std::cerr << "The model " << settings.GetType() << " does not belong to the possible models" << std::endl;
  
}