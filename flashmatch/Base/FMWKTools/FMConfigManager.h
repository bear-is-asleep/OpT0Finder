/**
 * \file FMConfigManager.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class flashmatch::FMConfigManager
 *
 * @author drinkingkazu
 */

/** \addtogroup core_Base

    @{*/
#ifndef __FLASHMATCHBASE_FMCONFIGMANAGER_H__
#define __FLASHMATCHBASE_FMCONFIGMANAGER_H__

#include <iostream>
#include "FMParams.h"
#include <set>

namespace flashmatch {
  /**
     \class FMConfigManager
     \brief Utility class to register a set of configurations
     Provides also a shared instance through which registered configurations can be shared beyond a single owner.\n
     Using flashmatch::FMParams, the uniqueness of configuration parameters is guaranteed (no worry to "overwrite")\n
  */
  class FMConfigManager {
    
  public:
    
    /// Default constructor
    FMConfigManager() {}
         
    /// Default destructor
    ~FMConfigManager(){}
    /// Shared static reference getter
    static const FMConfigManager& get() 
    {
      if(!_me) _me = new FMConfigManager;
      return *_me;
    }
    /// Adder of configuration from a file
    void AddConfigFile(const std::string cfg_file);
    /// Adder of configuration from parsed string
    void AddConfigString(const std::string cfg_str);
    /// Configuration retrieval method
    const flashmatch::FMParams& GetConfig(const std::string cfg);

  private:

    static FMConfigManager* _me;
    std::set<std::string> _cfg_files;
    flashmatch::FMParams _cfg;
    
  };
}
#endif
/** @} */ // end of doxygen group 

