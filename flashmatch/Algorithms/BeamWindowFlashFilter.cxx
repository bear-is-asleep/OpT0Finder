#ifndef OPT0FINDER_BEAMWINDOWFLASHFILTER_CXX
#define OPT0FINDER_BEAMWINDOWFLASHFILTER_CXX

#include "BeamWindowFlashFilter.h"
#include "TimeRange.h"
#include <map>
#include <numeric>
//#include <functional>
namespace flashmatch {

  static BeamWindowFlashFilterFactory __global_BeamWindowFlashFilterFactory__;

  BeamWindowFlashFilter::BeamWindowFlashFilter(const std::string name)
    : BaseFlashFilter(name)
    , _window_low(0.0)     // Default time [ns] as a lower threshold
    , _window_high(1600)   // Default time [ns] as an upper threshold
    , _window_tol(200)     // Default time [ns] as a padding around thresholds
    , _npe_threshold(10)   // Default # p.e. as a threshold
  {}

  void BeamWindowFlashFilter::_Configure_(const Config_t &pset)
  {
    std::vector<double> window = pset.get<std::vector<double>>("BeamWindow");
    _window_low = window[0];
    _window_high = window[1];
    _window_tol = pset.get<double>("WindowTolerance");
    _npe_threshold = pset.get<double>("NPEThreshold"  );
  }

  IDArray_t BeamWindowFlashFilter::Filter(const FlashArray_t& flash_v)
  {
    // Prepare a return flashmatch::IDArray_t object
    IDArray_t result;
    
    // Loop over flash array
    for(size_t index=0; index<flash_v.size(); ++index) {

      auto const& flash = flash_v[index]; // Retrieve this flash
      
      double npe = std::accumulate(flash.pe_v.begin(),flash.pe_v.end(),0.0); // Sum p.e.

      if(npe < _npe_threshold || flash.time < (_window_low - _window_tol) || flash.time > (_window_high + _window_tol)) continue; // Ignore if outside beam window

      result.push_back(index);

    }

    return result;
  }


}

#endif
