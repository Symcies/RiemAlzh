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
    LoadRealDataSettings(settings->FirstChildElement("real-data"));
    InitArgsOfOtherType(true);
  }
  else if (real_data == "false") {
    LoadSimulatedDataSettings(settings->FirstChildElement("simulated-data"));
    InitArgsOfOtherType(false);
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

void DataSettings::LoadRealDataSettings(const tinyxml2::XMLElement* settings) {

  is_data_real_ = true;
  data_path_ = settings->FirstChildElement("folder-path")->GetText();

  group_path_ = data_path_ + settings->FirstChildElement("group-file")->GetText();
  timepoints_path_ = data_path_ + settings->FirstChildElement("timepoints-file")->GetText();

  // TODO : For the path to data, path to timepoints and path to group, need to check if it is ok!!!
  /// Read the Cognitive scores
  const tinyxml2::XMLElement* data = settings->FirstChildElement("observations");
  const tinyxml2::XMLElement* cog_scores = data->FirstChildElement("cognitive-scores");
  const tinyxml2::XMLElement* landmarks  = data->FirstChildElement("landmarks");

  LoadRealCognitiveScores(cog_scores);
  LoadRealLandmarks(landmarks);



}


void DataSettings::LoadRealCognitiveScores(const tinyxml2::XMLElement *settings) {

  std::string Presence = settings->FirstChildElement("presence")->GetText();

  if(Presence == "yes")
  {
    are_cog_scores_present_ = true;
    cog_scores_path_ = data_path_ + settings->FirstChildElement("path-to-data")->GetText();

    if (!std::ifstream(cog_scores_path_))
    {
      std::cerr << "The file located at " << cog_scores_path_ << " does not exist" << std::endl;
    }

    cog_scores_dim_ = atoi(settings->FirstChildElement("dimension")->GetText());


    std::cout << "The model is reading cognitive scores located at " << cog_scores_path_ << std::endl;
  }
  else
  {
    are_cog_scores_present_ = false;
  }
}

void DataSettings::LoadRealLandmarks(const tinyxml2::XMLElement *settings)
{
  std::string Presence = settings->FirstChildElement("presence")->GetText();

  if(Presence == "yes")
  {
    are_landmarks_present_ = true;
    landmarks_path = data_path_ + settings->FirstChildElement("path-to-data")->GetText();

    if(!std::ifstream(landmarks_path))
    {
      std::cerr << "The file located at " << landmarks_path << " does not exist" << std::endl;
    }

    landmarks_dim_ = atoi(settings->FirstChildElement("dimension")->GetText());

    std::cout << "The model is reading landmarks located at " << landmarks_path << std::endl;
  }
  else
  {
    are_landmarks_present_ = false;
    landmarks_dim_ = -1;
  }
}


void DataSettings::LoadSimulatedDataSettings(const tinyxml2::XMLElement* settings) {

  is_data_real_ = false;

  subjects_total_num_  = atoi(settings->FirstChildElement("number-of-individuals")->GetText());
  min_observation_num_ = atoi(settings->FirstChildElement("min-number-of-observations")->GetText());
  max_observation_num_ = atoi(settings->FirstChildElement("max-number-of-observations")->GetText());

  std::cout << "The model is simulating between ";
  std::cout << min_observation_num_ << " and " << max_observation_num_;
  std::cout << " observations for " << subjects_total_num_ << " subjects" << std::endl;

}

void DataSettings::InitArgsOfOtherType(bool real){
  if (real){
    subjects_total_num_ = -1;
    min_observation_num_ = -1;
    max_observation_num_ = -1;
  }
  else {
    are_cog_scores_present_ = false;
    are_landmarks_present_ = false;
    cog_scores_dim_ = -1;
    landmarks_dim_ = -1;
  }
}

} //end namespace
