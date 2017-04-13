#pragma once

#include "DataSettings.h"
#include "RealDataSettings.h"
#include "SimulatedDataSettings.h"


class Builder {
 public:
  
  static std::shared_ptr<io::DataSettings> BuilderDataSettings(std::string xml_file);
};

