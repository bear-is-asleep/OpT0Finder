#ifndef __FMPARAMS_UTILS_CXX__
#define __FMPARAMS_UTILS_CXX__

#include "FMParamsUtils.h"
#include <sstream>
#include <fstream>
namespace flashmatch {

  std::string FMConfigFile2String(std::string fname)
  {
    std::ifstream filestrm(fname.c_str());
    std::string   contents;
    std::string   line;

    while(std::getline(filestrm, line)) {

      if(line.empty()) continue;
      
      std::stringstream linestrm(line);
      std::string       valid_line;
      
      std::getline(linestrm, valid_line, '#');
      
      if(valid_line.empty()) continue;
      
      contents += " " + valid_line;
    }
    filestrm.close();    
    return contents;
  }

  FMParams CreateFMParamsFromFile(std::string fname,std::string cfg_name)
  {
    FMParams res(cfg_name,FMConfigFile2String(fname));
    return res;
  }

}

#endif
