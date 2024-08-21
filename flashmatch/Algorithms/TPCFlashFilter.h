/**
 * \file TPCFlashFilter.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class TPCFlashFilter
 *
 * @author kazuhiro
 */

/** \addtogroup Algorithms

    @{*/
#ifndef OPT0FINDER_TPCFLASHFILTER_H
#define OPT0FINDER_TPCFLASHFILTER_H

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
     \class TPCFlashFilter
     Type flashmatch::BaseFlashFilter algorithm class. It filters out flashes that are not within the \n
     specified TPC volume.
  */
  class TPCFlashFilter : public BaseFlashFilter{
    
  public:
    
    /// Default constructor
    TPCFlashFilter(const std::string name="TPCFlashFilter");
    
    /// Default destructor
    ~TPCFlashFilter(){}

    /// Implementation of a virtual function
    IDArray_t Filter(const FlashArray_t&);

  protected:

    void _Configure_(const Config_t &pset);

  private:
    
  };

  /**
     \class flashmatch::TPCFlashFilterFactory
  */
  class TPCFlashFilterFactory : public FlashFilterFactoryBase {
  public:
    /// ctor
    TPCFlashFilterFactory() { FlashFilterFactory::get().add_factory("TPCFlashFilter",this); }
    /// dtor
    ~TPCFlashFilterFactory() {}
    /// creation method
    BaseFlashFilter* create(const std::string instance_name) { return new TPCFlashFilter(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 