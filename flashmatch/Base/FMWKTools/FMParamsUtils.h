/**
 * \file FMParamsUtils.h
 *
 * \ingroup FMParams
 * 
 * \brief Utility functions in Base/FMParams
 *
 * @author Kazu - Nevis 2015
 */

/** \addtogroup FMParams

    @{*/

#ifndef __FMPARAMS_UTILS_H__
#define __FMPARAMS_UTILS_H__

#include "FMParams.h"

namespace flashmatch {

  /// Given a configuration string, format to create flashmatch::FMParams
  //std::string FormatFMParamsString(std::string fname);
  /// Given a configuration file (full path), read & parse contents to create flashmatch::FMParams
  std::string FMConfigFile2String(std::string fname);
  /// Given a configuration file (full path), create and return flashmatch::FMParams
  FMParams CreateFMParamsFromFile(std::string fname,std::string cfg_name="cfg");

}

#endif
/** @} */ // end of doxygen group
