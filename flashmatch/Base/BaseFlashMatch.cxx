#ifndef OPT0FINDER_BASEFLASHMATCH_CXX
#define OPT0FINDER_BASEFLASHMATCH_CXX

#include "BaseFlashMatch.h"

namespace flashmatch {

  Flash_t BaseFlashMatch::GetEstimate(const QCluster_t& tpc) const
  {
    return _flash_hypothesis->GetEstimate(tpc);
  }

  void BaseFlashMatch::FillEstimate(const QCluster_t& tpc, Flash_t& opdet) const
  {
    _flash_hypothesis->FillEstimate(tpc,opdet);
  }

  void BaseFlashMatch::SetFlashHypothesis(flashmatch::BaseFlashHypothesis* alg)
  {
    _flash_hypothesis = alg;
  }

  void BaseFlashMatch::SetTPCCryo(int tpc, int cryo) {
    _tpc = tpc;
    _cryo = cryo;

    auto const& bbox = DetectorSpecs::GetME().ActiveVolume(_tpc, _cryo);
    _vol_xmax = bbox.Max()[0];
    _vol_xmin = bbox.Min()[0];
  }

  void BaseFlashMatch::InitializeMask(QCluster_t& trk) const{
    trk.tpc_mask_v.clear();
    trk.tpc_mask_v.resize(DetectorSpecs::GetME().NOpDets(), 0); // 1 means skip
  }

}

#endif
