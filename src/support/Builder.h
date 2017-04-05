#pragma once

#include "DataSettings.h"
#include "RealDataSettings.h"
#include "SimulatedDataSettings.h"


class Builder {
 public:
  
  static std::shared_ptr<io::DataSettings> BuilderDataSettings(const char *xml_file);
};

