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
      FLASH_DEBUG() << "Channel " << i << " is masked." << std::endl;
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

  int BaseFlashHypothesis::InspectTouchingEdges(const QCluster_t& trk) const
    {
        // Return code: 0=no touch, 1=start touching, 2=end touching, 3= both touching
        int result = 0;

        // Check if the start/end is near the edge
        auto const& start = trk.front();
        auto const& end   = trk.back();
        if (((start.x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_proximity) ||
            ((start.y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_proximity) ||
            ((start.z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - start.x) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - start.y) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - start.z) < _threshold_proximity)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            result += 1;
        }
        if (((end.x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_proximity) ||
            ((end.y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_proximity) ||
            ((end.z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - end.x) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - end.y) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - end.z) < _threshold_proximity)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            result += 2;
        }
        /*
        // OLD CODE looking at ANY POINT close to the edge
        double min_x = kINVALID_DOUBLE; double max_x = -kINVALID_DOUBLE;
        min_idx = max_idx = 0;
        for (size_t pt_index = 0; pt_index < trk.size(); ++pt_index) {
            if (trk[pt_index].x < min_x) {
                min_x = trk[pt_index].x;
                min_idx = pt_index;
            }
            if (trk[pt_index].x > max_x) {
                max_x = trk[pt_index].x;
                max_idx = pt_index;
            }
        }
        if (((trk[min_idx].x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_proximity) ||
            ((trk[min_idx].y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_proximity) ||
            ((trk[min_idx].z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - trk[max_idx].x) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - trk[max_idx].y) < _threshold_proximity) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - trk[max_idx].z) < _threshold_proximity)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            return true;
        }
        */
        return result;

    }

    QCluster_t BaseFlashHypothesis::ComputeExtension(const geoalgo::Vector& A, const geoalgo::Vector& B) const 
    {
        QCluster_t extension;

        // direction to extend should be A=>B ... so B should be closer to the edge
        geoalgo::Vector pt = B;
        geoalgo::Vector AB = B - A;
        auto length = AB.Length();
        if(length < 1.e-9) {
            FLASH_WARNING() << "Extension calculation halted as the direction estimate is unreliable for a short segment (" 
            << length << " cm)" << std::endl;
            return extension;
        }
        AB.Normalize();
        const auto& box0 = DetectorSpecs::GetME().PhotonLibraryVolume();
        const auto& box1 = DetectorSpecs::GetME().ActiveVolume();

        double dist_inside = 0;
        extension.reserve(100);
        pt += (AB * _segment_size/2.);
        while(box0.Contain(pt) && dist_inside < _threshold_proximity) {
          if(!box1.Contain(pt)) {
            QPoint_t qpt;
            qpt.x = pt[0];
            qpt.y = pt[1];
            qpt.z = pt[2];
            //std::cout << qpt.x << " " << qpt.y << " " << qpt.z << std::endl;
            qpt.q = _segment_size * DetectorSpecs::GetME().LightYield() * DetectorSpecs::GetME().MIPdEdx();
            extension.emplace_back(qpt);
          }
          else{ dist_inside += _segment_size; }

          pt += (AB * _segment_size);
        }

        return extension;
    }

    QCluster_t BaseFlashHypothesis::TrackExtension(const QCluster_t& in_trk, const int touch) const
    {
        auto start_touch = bool(touch & 0x1);
        auto end_touch   = bool(touch & 0x2);

        QCluster_t trk;
        if(start_touch) {
            geoalgo::Vector B(in_trk[0].x, in_trk[0].y, in_trk[0].z);
            geoalgo::Vector A(3);
            // search for a canddiate point to define a direction (and avoid using the same point)
            for(size_t i=1; i<in_trk.size(); ++i) {
                A[0] = in_trk[i].x;
                A[1] = in_trk[i].y;
                A[2] = in_trk[i].z;
                if(A.Dist(B) > _segment_size) break;
            }
            trk = this->ComputeExtension(A,B);
        }
        trk.reserve(trk.size() + in_trk.size());
        for(auto const& pt : in_trk) trk.push_back(pt);
        if(end_touch) {
            geoalgo::Vector B(in_trk[in_trk.size()-1].x, in_trk[in_trk.size()-1].y, in_trk[in_trk.size()-1].z);
            geoalgo::Vector A = B;
            // search for a canddiate point to define a direction (and avoid using the same point)
            for(int i=in_trk.size()-2; i>=0; --i) {
                A[0] = in_trk[i].x;
                A[1] = in_trk[i].y;
                A[2] = in_trk[i].z;
                if(A.Dist(B) > _segment_size) break;
            }
            auto end_trk = this->ComputeExtension(A,B);
            trk.reserve(trk.size()+end_trk.size());
            for(auto const& pt : end_trk)
                trk.push_back(pt);
        }

        /*

        QCluster_t trk = in_trk;
        // Extend the track
        // Compute coordinates of final point first
        geoalgo::Vector A(trk[max_idx].x, trk[max_idx].y, trk[max_idx].z);
        geoalgo::Vector B(trk[min_idx].x, trk[min_idx].y, trk[min_idx].z);
        geoalgo::Vector AB = B - A;

        double x_C = DetectorSpecs::GetME().PhotonLibraryVolume().Min()[0];
        double lengthBC = (x_C - B[0])/(B[0] - A[0]) * AB.Length();
        geoalgo::Vector C = B + AB / AB.Length() * lengthBC;

                // Add to _var_trk the part betwen boundary and C
                //_custom_algo->MakeQCluster(C, B, _var_trk, -1);
        geoalgo::Vector unit = (C - B).Dir();

        QPoint_t q_pt;
        geoalgo::Vector current = B;
        int num_pts = int(lengthBC / _segment_size);
        trk.reserve(trk.size() + num_pts);
        for (int i = 0; i < num_pts+1; i++) {
            double current_segment_size = (i < num_pts ? _segment_size : (lengthBC - _segment_size*num_pts));
            current = current + unit * current_segment_size/2.0;
            q_pt.x = current[0];
            q_pt.y = current[1];
            q_pt.z = current[2];
            q_pt.q = current_segment_size * DetectorSpecs::GetME().LightYield() * DetectorSpecs::GetME().MIPdEdx();
            if (trk.front().x < trk.back().x) {
                trk.insert(trk.begin(), q_pt);
            }
            else {
                trk.emplace_back(q_pt);
            }
            current = current + unit * current_segment_size/2.0;
            //std::cout << "Adding point " << current  << " " << i << " " << num_pts << std::endl;
        }
        //std::cout << " done " << trk.size() << std::endl;
        //std::cout << trk.front().x << " " << trk.back().x << std::endl;
        */
        return trk;
    }



}
#endif
