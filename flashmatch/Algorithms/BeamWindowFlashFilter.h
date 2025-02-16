/**
 * \file BeamWindowFlashFilter.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class BeamWindowFlashFilter
 *
 * @author Justin Mueller
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_BEAMWINDOWFLASHFILTER_H
#define OPT0FINDER_BEAMWINDOWFLASHFILTER_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 0
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashFilterFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashFilterFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

namespace flashmatch {
  /**
     \class BeamWindowFlashFilter
     Type flashmatch::BaseFlashFilter algorithm class. It applies a very simple  \n
     filter on the flash time using a beam window.
  */
  class BeamWindowFlashFilter : public BaseFlashFilter{
    
  public:
    
    /// Default constructor
    BeamWindowFlashFilter(const std::string name="BeamWindowFlashFilter");
    
    /// Default destructor
    ~BeamWindowFlashFilter(){}

    /// Implementation of a virtual function
    IDArray_t Filter(const FlashArray_t&);

  protected:

    void _Configure_(const Config_t &pset);

  private:

    double _window_low;    ///< threshold [ns]: to ignore any flash with time below this value.
    double _window_high;   ///< threshold [ns]: to ignore any flash with time above this value.
    double _window_tol;    ///< time [ns]: used as a padding on either of above limits.
    double _npe_threshold; ///< threshold [p.e.]: to ignore any flash below this value
    
  };

  /**
     \class flashmatch::BeamWindowFlashFilterFactory
  */
  class BeamWindowFlashFilterFactory : public FlashFilterFactoryBase {
  public:
    /// ctor
    BeamWindowFlashFilterFactory() { FlashFilterFactory::get().add_factory("BeamWindowFlashFilter",this); }
    /// dtor
    ~BeamWindowFlashFilterFactory() {}
    /// creation method
    BaseFlashFilter* create(const std::string instance_name) { return new BeamWindowFlashFilter(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 

