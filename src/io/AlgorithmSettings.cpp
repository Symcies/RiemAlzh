#include "AlgorithmSettings.h"

namespace io {


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////


AlgorithmSettings::AlgorithmSettings(char *xml_file) {
    tinyxml2::XMLDocument parameters;
    parameters.LoadFile(xml_file);

    auto settings = parameters.FirstChildElement("algorithm-settings");

    max_num_iter_ = atoi(
            settings->FirstChildElement("max-iterations")->GetText());
    num_burnin_iter_ = atoi(settings->FirstChildElement("burn-in")->GetText());
    output_iter_ = atoi(
            settings->FirstChildElement("step-size-to-display")->GetText());
    data_save_iter_ = atoi(settings->FirstChildElement("step-size-to-save")->GetText());
}


AlgorithmSettings::~AlgorithmSettings() {

}

////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////



}
