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
    num_burn_iter_ = atoi(settings->FirstChildElement("burn-in")->GetText());
    counter_to_next_output_display_ = atoi(
            settings->FirstChildElement("step-size-to-display")->GetText());
    counter_to_next_data_save_ = atoi(settings->FirstChildElement("step-size-to-save")->GetText());
}


AlgorithmSettings::~AlgorithmSettings() {

}

////////////////////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////////////////////



}
