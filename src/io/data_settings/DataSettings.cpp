#include "DataSettings.h"

namespace io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

DataSettings::DataSettings(const char *xml_file) {

  if (xml_file == ""){
    std::cerr << "Problem in DataSettings. xml_file name null." << std::endl;
    //TODO: define exit behavior
    return;
  }

  tinyxml2::XMLDocument file;
  file.LoadFile(xml_file);

  //TODO: check if auto is good practice or can generate crash
  auto settings = file.FirstChildElement("data-settings");

  std::string real_data = settings->FirstChildElement("data-type")->GetText();
  if (real_data == "true") {
    is_data_real_ = true;
  }
  else if (real_data == "false") {
    is_data_real_ = false;
  }
  else {
    //TODO: better error message
    std::cerr << "Problem in DataSettings (loading data)" << std::endl;
  }


}

DataSettings::~DataSettings() {

}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////



} //end namespace
