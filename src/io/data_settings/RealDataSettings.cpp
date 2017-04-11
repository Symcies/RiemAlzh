#include "RealDataSettings.h"


namespace  io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

RealDataSettings::RealDataSettings(std::string xml_file) : DataSettings(xml_file) {

  if (InputsAssert::IsFileCorrect(xml_file, true)){

    tinyxml2::XMLDocument file;
    file.LoadFile(xml_file.c_str());
    auto settings = file.FirstChildElement("data-settings")->FirstChildElement("real-data");

    /// Extract paths
    std::string data_path = settings->FirstChildElement("folder-path")->GetText();
    if (InputsAssert::IsFileCorrect(&data_path[0], false)){
      data_path_ = data_path;
    }
    std::string group_path = data_path_ + settings->FirstChildElement("group-file")->GetText();
    if (InputsAssert::IsFileCorrect(&group_path[0], false)){
      group_path_ = group_path;
    }
    std::string timepoints_path = data_path_ + settings->FirstChildElement("timepoints-file")->GetText();
    if (InputsAssert::IsFileCorrect(&timepoints_path[0], false)){
      timepoints_path_ = timepoints_path;
    }

    /// Read the Cognitive scores
    const tinyxml2::XMLElement *data = settings->FirstChildElement("observations");
    const tinyxml2::XMLElement *cog_scores = data->FirstChildElement("cognitive-scores");
    const tinyxml2::XMLElement *landmarks = data->FirstChildElement("landmarks");

    LoadCognitiveScores(cog_scores);
    LoadLandmarks(landmarks);
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////

void RealDataSettings::LoadLandmarks(const tinyxml2::XMLElement *settings) {
  std::string presence = settings->FirstChildElement("presence")->GetText();

  if(InputsAssert::ToLowerCase(presence) == "yes")
  {
    are_landmarks_present_ = true;

    std::string landmarks_path = data_path_ + settings->FirstChildElement("path-to-data")->GetText();
    if (InputsAssert::IsFileCorrect(landmarks_path, false)){
      landmarks_path_ = landmarks_path;
    }

    landmarks_dim_ = std::stoi(settings->FirstChildElement("dimension")->GetText());

    std::cout << "The model is reading landmarks located at " << landmarks_path << std::endl;
  }
  else if (InputsAssert::ToLowerCase(presence) == "no")
  {
    are_landmarks_present_ = false;
    landmarks_dim_ = 0;
    std::cout << "No landmarks." << std::endl;
  }
  else {
    std::cerr << "Presence must be 'yes' or 'no' in the DataSettings file. It was " << presence <<std::endl;
  }
}
    
void RealDataSettings::LoadCognitiveScores(const tinyxml2::XMLElement *settings) {
  std::string presence = settings->FirstChildElement("presence")->GetText();
  if(InputsAssert::ToLowerCase(presence) == "yes")
  {
    are_cog_scores_present_ = true;

    std::string cog_scores_path = data_path_ + settings->FirstChildElement("path-to-data")->GetText();
    if (InputsAssert::IsFileCorrect(cog_scores_path, false)){
      cog_scores_path_ = cog_scores_path;
    }

    cog_scores_dim_ = std::stoi(settings->FirstChildElement("dimension")->GetText());


    std::cout << "The model is reading cognitive scores located at " << cog_scores_path_ << std::endl;
  }
  else if (InputsAssert::ToLowerCase(presence) == "no")
  {
    are_cog_scores_present_ = false;
    cog_scores_dim_ = 0;
    std::cout << "No cog scores." << std::endl;

  }
  else {
    std::cerr << "Presence must be 'yes' or 'no' in the DataSettings file. It was " << presence <<std::endl;
  }
}


}




