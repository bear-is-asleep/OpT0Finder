#ifndef OPT0FINDER_TIMEYZCOMPATMATCH_CXX
#define OPT0FINDER_TIMEYZCOMPATMATCH_CXX

#include "TimeYZCompatMatch.h"
#include <cmath>
#include <sstream>

namespace flashmatch {

  static TimeYZCompatMatchFactory __global_TimeYZCompatMatchFactory__;

  TimeYZCompatMatch::TimeYZCompatMatch(const std::string name)
    : BaseProhibitAlgo(name)
  {}

  void TimeYZCompatMatch::_Configure_(const Config_t &pset)
  {
    _time_shift  = pset.get<double>("BeamTimeShift",0);
    _time_window = pset.get<double>("TouchingTrackWindow",5);
    _yz_distance = pset.get<double>("YZDistance", 100);
  }

  bool TimeYZCompatMatch::MatchCompatible(const QCluster_t& clus, const Flash_t& flash)
  {
    double qy_center(0), qz_center(0), totalq(0);
    for(auto& qpt : clus)
    {
      qy_center += qpt.q * qpt.y;
      qz_center += qpt.q * qpt.z;
      totalq += qpt.q;
    }
    qy_center /= (totalq * clus.size());
    qz_center /= (totalq * clus.size());
    bool near_yz(std::sqrt(std::pow(qy_center - flash.y, 2) + std::pow(qz_center - flash.z, 2)) < _yz_distance);

    if(clus.empty()) return false; 

    // get time of flash
    auto flash_time = flash.time - _time_shift;

    // get time of cluster by looking at the range of x-positions
    double clus_x_min = kINVALID_DOUBLE;
    double clus_x_max = -1 * kINVALID_DOUBLE;
    for (auto const& pt : clus){
      if (pt.x > clus_x_max) { clus_x_max = pt.x; }
      if (pt.x < clus_x_min) { clus_x_min = pt.x; }
    }

    // Detector boundary info
    double xmax = DetectorSpecs::GetME().ActiveVolume().Max()[0];
    double xmin = DetectorSpecs::GetME().ActiveVolume().Min()[0];

    // Assume this is the right flash, check if a trajectory is within the active volume
    // One assumes tpc0 and the other tpc1 
    double reco_x_tpc0 = clus_x_min - flash_time * DetectorSpecs::GetME().DriftVelocity();
    double reco_x_tpc1 = clus_x_max + flash_time * DetectorSpecs::GetME().DriftVelocity();

    double distance_window = _time_window * DetectorSpecs::GetME().DriftVelocity();

    bool incompatible_tpc0 = (xmin - reco_x_tpc0) > distance_window;
    bool incompatible_tpc1 = (reco_x_tpc1 - xmax) > distance_window;

    FLASH_INFO() << "Inspecting..." << std::endl
      << "Detector X span : " << xmin << " => " << xmax << std::endl
      << "TPC pts X span  : " << clus_x_min << " => " << clus_x_max << " ... " << clus.size() << " points" << std::endl
      << "Flash time      : " << flash_time << " (shifted by " << _time_shift << ")" << std::endl
      << "Hypothesis X pos: " << reco_x_tpc0 << " => " << reco_x_tpc1 << std::endl
      << "From TPC-0 edge : " << (xmin - reco_x_tpc0) << " ... incompatible? " << incompatible_tpc0 << std::endl
      << "From TPC-1 edge : " << (reco_x_tpc1 - xmax) << " ... incompatible? " << incompatible_tpc1 << std::endl;

    if ( (incompatible_tpc0 && incompatible_tpc1) || !near_yz ) return false;
    return true;

  }


}
#endif
