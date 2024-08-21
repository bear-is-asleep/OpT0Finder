#ifndef OPT0FINDER_TPCFLASHFILTER_CXX
#define OPT0FINDER_TPCFLASHFILTER_CXX

#include "TPCFlashFilter.h"
#include <map>
#include <numeric>
//#include <functional>
namespace flashmatch {

  static TPCFlashFilterFactory __global_TPCFlashFilterFactory__;

  TPCFlashFilter::TPCFlashFilter(const std::string name)
    : BaseFlashFilter(name)
  {}

  void TPCFlashFilter::_Configure_(const Config_t &pset)
  {

  }

  IDArray_t TPCFlashFilter::Filter(const FlashArray_t& flash_v)
  {
    
  }


}

#endif