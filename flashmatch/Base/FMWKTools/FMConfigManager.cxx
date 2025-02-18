#ifndef __FLASHMATCHBASE_FMCONFIGMANAGER_CXX__
#define __FLASHMATCHBASE_FMCONFIGMANAGER_CXX__

#include "FMConfigManager.h"
#include "FMParamsUtils.h"
namespace flashmatch {

  FMConfigManager* FMConfigManager::_me = nullptr;

  void FMConfigManager::AddConfigFile(const std::string cfg_file)
  {
    if(_cfg_files.find(cfg_file)!=_cfg_files.end()) {
      std::cerr << "Duplicate addition of config fiel: " << cfg_file << std::endl;
      throw std::exception();
    }

    _cfg.add_pset(CreateFMParamsFromFile(cfg_file));
  }

  void FMConfigManager::AddConfigString(const std::string cfg_str)
  {
    FMParams p;
    p.add(cfg_str);
    _cfg.add_pset(p);
  }
  
  const FMParams& FMConfigManager::GetConfig(const std::string cfg)
  {
    return _cfg.get_pset(cfg);
  }

}

#endif
