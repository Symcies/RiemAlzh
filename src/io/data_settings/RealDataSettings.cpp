#include "RealDataSettings.h"


namespace  io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

RealDataSettings::RealDataSettings(const char *xml_file) : DataSettings(xml_file) {
  
  tinyxml2::XMLDocument file;
  file.LoadFile(xml_file);
  auto settings = file.FirstChildElement("data-settings")->FirstChildElement("real-data");
  data_path_ = settings->FirstChildElement("folder-path")->GetText();

  group_path_ = data_path_ + settings->FirstChildElement("group-file")->GetText();
  timepoints_path_ = data_path_ + settings->FirstChildElement("timepoints-file")->GetText();

  // TODO : For the path to data, path to timepoints and path to group, need to check if it is ok!!!
  /// Read the Cognitive scores
  const tinyxml2::XMLElement* data = settings->FirstChildElement("observations");
  const tinyxml2::XMLElement* cog_scores = data->FirstChildElement("cognitive-scores");
  const tinyxml2::XMLElement* landmarks  = data->FirstChildElement("landmarks");

  LoadCognitiveScores(cog_scores);
  LoadLandmarks(landmarks);
  
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void RealDataSettings::LoadLandmarks(const tinyxml2::XMLElement *settings) {
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
    landmarks_dim_ = 0;
  }
}
    
void RealDataSettings::LoadCognitiveScores(const tinyxml2::XMLElement *settings) { 
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


}




