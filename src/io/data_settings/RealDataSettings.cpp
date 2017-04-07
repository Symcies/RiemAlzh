#include "RealDataSettings.h"


namespace  io {


////////////////////////////////////////////////////////////////////////////////////////////////
/// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////

RealDataSettings::RealDataSettings(const char *xml_file) : DataSettings(xml_file) {

  if (InputsAssert::IsFilePathCorrect(xml_file) && InputsAssert::IsXMLValid(xml_file)){

    tinyxml2::XMLDocument file;
    file.LoadFile(xml_file);
    auto settings = file.FirstChildElement("data-settings")->FirstChildElement("real-data");

    std::string data_path = settings->FirstChildElement("folder-path")->GetText();
    if (InputsAssert::IsFilePathCorrect(&data_path[0])){
      data_path_ = data_path;
    }
    std::string group_path = data_path_ + settings->FirstChildElement("group-file")->GetText();
    if (InputsAssert::IsFilePathCorrect(&group_path[0])){
      group_path_ = group_path;
    }
    std::string timepoints_path = data_path_ + settings->FirstChildElement("timepoints-file")->GetText();
    if (InputsAssert::IsFilePathCorrect(&timepoints_path[0])){
      timepoints_path_ = timepoints_path;
    }

    // TODO : For the path to data, path to timepoints and path to group, need to check if it is ok!!!
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
    if (InputsAssert::IsFilePathCorrect(landmarks_path)){
      landmarks_path_ = landmarks_path;
    }

    landmarks_dim_ = atoi(settings->FirstChildElement("dimension")->GetText());

    std::cout << "The model is reading landmarks located at " << landmarks_path << std::endl;
  }
  else if (InputsAssert::ToLowerCase(presence) == "no")
  {
    are_landmarks_present_ = false;
    landmarks_dim_ = 0;
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
    if (InputsAssert::IsFilePathCorrect(cog_scores_path)){
      cog_scores_path_ = cog_scores_path;
    }

    cog_scores_dim_ = atoi(settings->FirstChildElement("dimension")->GetText());


    std::cout << "The model is reading cognitive scores located at " << cog_scores_path_ << std::endl;
  }
  else if (InputsAssert::ToLowerCase(presence) == "no")
  {
    are_cog_scores_present_ = false;
    landmarks_dim_ = 0;
  }
  else {
    std::cerr << "Presence must be 'yes' or 'no' in the DataSettings file. It was " << presence <<std::endl;
  }
}


}




