#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX

#include "PhotonLibHypothesis.h"
#include <assert.h>

#ifndef USING_LARSOFT
#define USING_LARSOFT 0
#endif

#if USING_LARSOFT == 0
#include <omp.h>
#define NUM_THREADS 12
#endif

//using namespace std::chrono;

namespace flashmatch {

    static PhotonLibHypothesisFactory __global_PhotonLibHypothesisFactory__;

    PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
    {
        #if USING_LARSOFT == 0
        omp_set_num_threads(NUM_THREADS); 
        #endif
    }

    void PhotonLibHypothesis::_Configure_(const Config_t &pset)
    {
        // Call the base class's _Configure_ method
        BaseFlashHypothesis::_Configure_(pset);

        _reco_pe_calib = pset.get<double>("RecoPECalibFactor",1.0);
        _extend_tracks = pset.get<bool>("ExtendTracks", false);
        /*
        if(_extend_tracks) {
            FLASH_CRITICAL() << "ExtendTracks currently has a problem and cannot be enabled!" << std::endl
            << "If you would like to debug, here is useful info..." << std::endl
            << "1. InspectTouchingEdges function checks if a track is touching either of 2 PMT planes (within 1 cryostat)" << std::endl
            << "2. But TrackExtension function only extends to the negative x direction (and it does so disregard of the fact if that end is touching or not)" << std::endl
            << "This is completely wrong (a bug). Hence we disabled at the moment. TrackExtension may not be critical for flash matching performance" <<std::endl;
            throw OpT0FinderException();
        } 
        */
        _threshold_track_len = pset.get<double>("ExtensionTrackLengthMaxThreshold", 20.0);

    }


    void PhotonLibHypothesis::FillEstimate(const QCluster_t& tpc_trk, Flash_t &flash) const
    {
        size_t n_pmt = DetectorSpecs::GetME().NOpDets();//n_pmt returns 0 now, needs to be fixed
        if(flash.pe_v.empty()) flash.pe_v.resize(n_pmt);
        if(flash.pe_err_v.empty()) flash.pe_err_v.resize(n_pmt);
        if(flash.pe_true_v.empty()) flash.pe_true_v.resize(n_pmt);

        assert(flash.pe_v.size()     == n_pmt);
        assert(flash.pe_true_v.size() == n_pmt);
        assert(flash.pe_err_v.size() == n_pmt);

        for (auto& v : flash.pe_v      ) {v = 0;}
        for (auto& v : flash.pe_err_v  ) {v = 0;}
        for (auto& v : flash.pe_true_v ) {v = 0;}

        double track_length = tpc_trk.front().dist(tpc_trk.back());
        int touch = this->InspectTouchingEdges(tpc_trk);
        bool extend_tracks = (_extend_tracks && touch && track_length < _threshold_track_len);
        if(_extend_tracks) {
            FLASH_DEBUG() << "Extend? " << extend_tracks << " ... track length " << track_length << " touch-or-not " << touch << std::endl;
        }


        if(!extend_tracks) 
            this->BuildHypothesis(tpc_trk,flash);
        else {
            auto trk = this->TrackExtension(tpc_trk,touch);
            this->BuildHypothesis(trk,flash);
        }

    }

    void PhotonLibHypothesis::BuildHypothesis(const QCluster_t& trk, Flash_t &flash) const
    {
        //bool speak=false;
        //std::cout<<trk[0].x<<" "<<trk[0].y<<" "<<trk[0].z<<std::endl;
        //if(std::fabs(trk[0].x + 211.635)<1.0) speak=true;
        size_t n_pmt = DetectorSpecs::GetME().NOpDets();
	#if USING_LARSOFT == 0
        #pragma omp parallel
        #endif
        {
            size_t num_pts  = trk.size();
            size_t start_pt = 0;

            #if USING_LARSOFT == 0
            size_t thread_id = omp_get_thread_num();
            size_t num_threads = omp_get_num_threads();
            if(num_threads == 1 || num_threads>num_pts) {
                if(thread_id > 0) {
                    start_pt = 0;
                }
            }else{
                num_pts = trk.size() / (num_threads-1);
                start_pt = num_pts * thread_id;
                if(thread_id+1 == num_threads) num_pts = (trk.size() % (num_threads-1));
                //sleep(thread_id);
                //std::cout<<"Start " <<start_pt << " num pts " << num_pts << " / " << trk.size() << std::endl;
            }
            #endif

            auto const& vox_def = DetectorSpecs::GetME().GetVoxelDef();
            std::vector<double> local_pe_v(n_pmt,0.);
            std::vector<double> local_pe_refl_v(n_pmt,0.);
            double pos[3];
            size_t nproc0,nproc1; // for debug cout
            nproc0 = nproc1 = 0;  // for debug cout

            for(size_t ipt = start_pt; ipt < start_pt+num_pts; ++ipt) {

                auto const& pt = trk[ipt];
                pos[0] = pt.x;
                pos[1] = pt.y;
                pos[2] = pt.z;
                nproc0 += 1; // for debug cout

                int vox_id = vox_def.GetVoxelID(pos);
                if (vox_id < 0) continue;

		auto const& lib_data = DetectorSpecs::GetME().GetLibraryEntries(vox_id);

		#if USING_LARSOFT == 1
		auto const& lib_data_refl = DetectorSpecs::GetME().GetLibraryEntries(vox_id,true);
		#endif

                nproc1 += 1; // for debug cout
                double qsum = 0.;
                double vsum = 0.;
                for(size_t ipmt=0; ipmt < n_pmt; ++ipmt) {

                    if(_channel_mask[ipmt]) continue;

                    if(!_uncoated_pmt_list[ipmt])
		      local_pe_v[ipmt] += pt.q * lib_data[ipmt];
            qsum += pt.q * lib_data[ipmt];
            vsum += lib_data[ipmt];
		    #if USING_LARSOFT == 1
		    if(_global_qe_refl > 0.)
		      local_pe_refl_v[ipmt] += pt.q * lib_data_refl[ipmt];
                    #endif
                    //local_pe_v[ipmt] += pt.q * vis_pmt[ipmt];
                    //std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q0 << "," << q1 << "]" << std::endl;
                }
                //if(speak) std::cout << pt.x << " " << pt.y << " " << pt.z << " ... " << pt.q << ", " << qsum << ", " << vsum << std::endl;
            }

            #if USING_LARSOFT == 0
            #pragma omp critical
            #endif
            {
                double qsum = 0.; // for debug cout
                for(size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
                    double q0 = (local_pe_v[ipmt] * _global_qe * _reco_pe_calib / _qe_v[ipmt]);
                    double q1 = (local_pe_refl_v[ipmt] * _global_qe_refl * _reco_pe_calib / _qe_v[ipmt]);
                    double q = q0 + q1;
                    flash.pe_v[ipmt] +=  q;
                    qsum += q; // for debug cout
                }

                // for debug cout
                //double qsum2 = 0.;
                //for(auto const& v: flash.pe_v) qsum2 += v;
                //std::cout<<"Thread ID " << thread_id << " ... " << start_pt << " => " << start_pt+num_pts << " ... " << nproc0 << "/" << nproc1 << " qsum " << qsum << " total sum " << qsum2 << std::endl<<std::flush;
                //std::cout<<qsum<<std::endl;
            }

        }
        // for debug cout
        //double charge = 0.;
        //for(auto const& v : flash.pe_v) { charge+=v; }
        //if(charge>0)
        //    std::cout << std::endl << "Track size " << trk.size() << " points, total pe " << charge << std::endl << std::flush;
        //sleep(3);
        //throw OpT0FinderException();
        return;
    }

}
#endif