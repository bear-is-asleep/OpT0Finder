
#ifndef SEMIANALYTICALMODEL_H
#define SEMIANALYTICALMODEL_H

// SemiAnalyticalModel
//  - fast optical simulation using semi-analytical model
//  - contains functions to calculate the number of direct and reflected photons
//  incident
//    each photo-detector, along with the necessary ultility functions (geometry
//    calculations etc.)
//  - full description of model: Eur. Phys. J. C 81, 349 (2021)

// Nov 2021 by P. Green

#include "TVector3.h"

//#include "boost/math/policies/policy.hpp"

#include <map>
#include <vector>

//#define __STDCPP_WANT_MATH_SPEC_FUNCS__


#ifndef USING_LARSOFT
#define USING_LARSOFT 0
#endif

#if USING_LARSOFT == 0
//Geometry
#include "flashmatch/GeoAlgo/GeoVector.h"
#include "flashmatch/GeoAlgo/GeoAABox.h"
#include "flashmatch/GeoAlgo/GeoAlgoConstants.h"
//Config
//#include "flashmatch/Base/FMWKTools/ConfigManager.h"
#include "flashmatch/Base/FMWKInterface.h"
//#include "flashmatch/Base/FMWKTools/PSetUtils.h"
#include "flashmatch/Base/OpT0FinderException.h"
#include "flashmatch/Base/BaseFlashHypothesis.h"
#include "flashmatch/Base/FlashHypothesisFactory.h"
//Utils 
//#include "counter.h"
namespace flashmatch {
  /// Configuration object
  using Config_t = flashmatch::PSet;
}
#else
//Config
#include "fhiclcpp/ParameterSet.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashHypothesisFactory.h"
namespace flashmatch {
  /// Configuration object
  using Config_t = fhicl::ParameterSet;
}
#endif



// Define a new policy *not* internally promoting RealType to double:
// typedef boost::math::policies::policy<boost::math::policies::promote_double<false>>
//  noLDoublePromote;

namespace flashmatch {
  class SemiAnalyticalModel : public BaseFlashHypothesis {

  public:
    // Default constructor
    SemiAnalyticalModel(const std::string name="SemiAnalyticalModel");

    // Default destructor
    ~SemiAnalyticalModel() {};

    void FillEstimate(const QCluster_t&, Flash_t&) const;

    // direct / VUV light
    void detectedDirectVisibilities(std::vector<double>& DetectedVisibilities,
                                    geoalgo::Point_t const& ScintPoint) const;

    // reflected / visible light
    void detectedReflectedVisibilities(std::vector<double>& ReflDetectedVisibilities,
                                       geoalgo::Point_t const& ScintPoint,
                                       bool AnodeMode = false) const;

    // print the configuration
    void Print() const;
    void printVUVParameters() const;
    void printVISParameters() const;

    // structure for optical detector information
    struct OpticalDetector {
      double h; // height
      double w; // width
      geoalgo::Point_t center;
      int type;
      int orientation;
    };
  protected:

    void _Configure_(const Config_t& pset);

  private:
    double VUVAbsorptionLength() const;

    // structure for rectangular solid angle calculation
    struct Dims {
      double h, w; // height, width
    };


    // direct light photo-detector visibility calculation
    double VUVVisibility(geoalgo::Point_t const& ScintPoint, OpticalDetector const& opDet) const;

    // reflected light photo-detector visibility calculation
    double VISVisibility(geoalgo::Point_t const& ScintPoint,
                         OpticalDetector const& opDet,
                         const double cathode_visibility,
                         geoalgo::Point_t const& hotspot,
                         bool AnodeMode = false) const;

    // Gaisser-Hillas
    double Gaisser_Hillas(const double x, const double* par) const;

    // solid angle calculations
    // rectangular aperture
    double Rectangle_SolidAngle(const double a, const double b, const double d) const;
    double Rectangle_SolidAngle(Dims const& o,
                                geoalgo::Vector_t const& v,
                                const int OpDetOrientation) const;
    // circular aperture
    double Disk_SolidAngle(const double d, const double h, const double b) const;
    // dome aperture calculation
    double Omega_Dome_Model(const double distance, const double theta) const;

    // TODO: replace with geometry service
    bool isOpDetInSameTPC(geoalgo::Point_t const& ScintPoint, geoalgo::Point_t const& OpDetPoint) const;
    std::vector<OpticalDetector> opticalDetectors() const;

    // geometry properties
    //geo::GeometryCore const& fGeom;
    //const larg4::ISTPC fISTPC;
    int fNTPC;
    geoalgo::AABox fActiveVolume;
    TVector3 fcathode_centre, fanode_centre;
    double fplane_depth, fanode_plane_depth;

    // photodetector geometry properties
    size_t fNOpDets;
    std::vector<OpticalDetector> fOpDetector;
    double fradius;
    Dims fcathode_plane;
    Dims fanode_plane;

    // optical detector properties inherited from FWKMInterface
    int fspherical_type;
    int fspherical_orientation;
    std::vector<int> fspherical_ids;
    int frectengular_type;
    int frectengular_orientation;
    double frectengular_height;
    double frectengular_width;
    std::vector<int> frectengular_ids;

    // For VUV semi-analytic hits
    double fdelta_angulo_vuv;
    // flat PDs
    bool fIsFlatPDCorr;
    std::vector<std::vector<double>> fGHvuvpars_flat;
    std::vector<double> fborder_corr_angulo_flat;
    std::vector<std::vector<double>> fborder_corr_flat;
    // lateral PDs
    bool fIsFlatPDCorrLat;
    std::vector<std::vector<double>> fGHvuvpars_flat_lateral;
    std::vector<double> fborder_corr_angulo_flat_lateral;
    std::vector<std::vector<double>> fborder_corr_flat_lateral;

    // dome PDs
    bool fIsDomePDCorr;
    std::vector<std::vector<double>> fGHvuvpars_dome;
    std::vector<double> fborder_corr_angulo_dome;
    std::vector<std::vector<double>> fborder_corr_dome;
    // Field cage scaling
    bool fApplyFieldCageTransparency;
    double fFieldCageTransparencyLateral;
    double fFieldCageTransparencyCathode;

    // For VIS semi-analytic hits
    bool fDoReflectedLight;
    bool fIncludeAnodeReflections;
    // correction parameters for VIS Nhits estimation
    double fdelta_angulo_vis;
    double fAnodeReflectivity;
    // flat PDs
    std::vector<double> fvis_distances_x_flat;
    std::vector<double> fvis_distances_r_flat;
    std::vector<std::vector<std::vector<double>>> fvispars_flat;
    // lateral PDs
    std::vector<double> fvis_distances_x_flat_lateral;
    std::vector<double> fvis_distances_r_flat_lateral;
    std::vector<std::vector<std::vector<double>>> fvispars_flat_lateral;
    // dome PDs
    std::vector<double> fvis_distances_x_dome;
    std::vector<double> fvis_distances_r_dome;
    std::vector<std::vector<std::vector<double>>> fvispars_dome;

    // absorption length
    bool fUseXeAbsorption;
    double fvuv_absorption_length;
    std::map<double, double> abs_length_spectrum;

    std::vector<double> _qe_refl_v;

    //std::unique_ptr<flashmatch::SemiAnalyticalModel> _semi_model; ///< The semi-analytical model
  };

  /**
    \class flashmatch::SemiAnalyticalModelFactory
  */
  class SemiAnalyticalModelFactory : public FlashHypothesisFactoryBase {
  public:
    /// ctor
    SemiAnalyticalModelFactory() { FlashHypothesisFactory::get().add_factory("SemiAnalyticalModel",this); }
    /// dtor
    ~SemiAnalyticalModelFactory() {}
    /// creation method
    BaseFlashHypothesis* create(const std::string instance_name) { return new SemiAnalyticalModel(instance_name); }
  };

} // namespace

#endif