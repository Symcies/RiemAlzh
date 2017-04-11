#include "Builder.h"

 
std::shared_ptr<io::DataSettings> Builder::BuilderDataSettings(std::string xml_file) {
  if (xml_file == ""){
    std::cerr << "Problem in DataSettings. xml_file name null." << std::endl;
    //TODO: define exit behavior
    
  }
  
  tinyxml2::XMLDocument file;
  file.LoadFile(xml_file.c_str());
  
  auto settings = file.FirstChildElement("data-settings");
  std::string real_data = settings->FirstChildElement("data-type")->GetText();
  
  if (real_data == "true") {
    return std::make_shared<io::RealDataSettings>(xml_file);
  } else {
    return std::make_shared<io::SimulatedDataSettings>(xml_file);
  }
}