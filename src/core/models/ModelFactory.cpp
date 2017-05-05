#include "ModelFactory.h"


std::shared_ptr<AbstractModel> ModelFactory::NewModel(io::ModelSettings &settings) {
  
  if (settings.GetType() == "Univariate") 
    return std::make_shared<UnivariateModel>(settings);
  else if (settings.GetType() == "Multivariate")
    return std::make_shared<MultivariateModel>(settings);
  else if (settings.GetType() == "FastNetwork")
    return std::make_shared<FastNetworkModel>(settings);
  else if (settings.GetType() == "Meshwork")
    return std::make_shared<MeshworkModel>(settings);
  else if (settings.GetType() == "Gaussian")
    return std::make_shared<GaussianModel>(settings);
  else if(settings.GetType() == "GaussianMixture")
    return std::make_shared<GaussianMixtureModel>(settings);
  else 
    std::cerr << "The model " << settings.GetType() << " does not belong to the possible models" << std::endl;
  
}