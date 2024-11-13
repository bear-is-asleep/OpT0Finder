#ifndef BASEFLASHHYPOTHESIS_CXX
#define BASEFLASHHYPOTHESIS_CXX

#include "BaseFlashHypothesis.h"
#include "OpT0FinderException.h"
namespace flashmatch {

  BaseFlashHypothesis::BaseFlashHypothesis(std::string name)
    : BaseAlgorithm(kFlashHypothesis,name)
    , _channel_mask(DetectorSpecs::GetME().NOpDets(),false)
    , _uncoated_pmt_list(DetectorSpecs::GetME().NOpDets(),false)
  {}

  void BaseFlashHypothesis::_Configure_(const Config_t &pset)
  {
    _threshold_proximity = pset.get<double>("ProximityThreshold",5.0);
    _segment_size = pset.get<double>("SegmentSize",0.5);

    _global_qe = pset.get<double>("GlobalQE");
    _global_qe_refl = pset.get<double>("GlobalQERefl", -1);

    _qe_v.clear();
    _qe_v = pset.get<std::vector<double> >("CCVCorrection",_qe_v);
    if(_qe_v.empty()) _qe_v.resize(DetectorSpecs::GetME().NOpDets(),1.0);
    if(_qe_v.size() != DetectorSpecs::GetME().NOpDets()) {
      FLASH_CRITICAL() << "CCVCorrection factor array has size " << _qe_v.size()
      << " != number of opdet (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
      throw OpT0FinderException();
    }

    //Debug statements
    for(size_t i=0; i<_channel_mask.size(); ++i) {
      FLASH_DEBUG() << "Channel " << i << " is masked to " << _channel_mask[i] << std::endl;
    }

  }

  //TODO: Implement this mask for a given list of channel ids or conditions
  void BaseFlashHypothesis::InitializeMask(Flash_t &flash) const{
    flash.pds_mask_v.clear();
    flash.pds_mask_v.resize(DetectorSpecs::GetME().NOpDets(), 0); // 1 means skip
  }

  Flash_t BaseFlashHypothesis::GetEstimate(const QCluster_t& tpc) const
  {
    Flash_t res;
    //res.pe_v.resize(OpDetXArray().size());

    FillEstimate(tpc,res);
    return res;
  }

}
#endif
