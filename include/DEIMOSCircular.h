#ifndef SHAPELENS_DEIMOSCIRCULAR_H
#define SHAPELENS_DEIMOSCIRCULAR_H

#include "DEIMOS.h"
#include "History.h"

namespace shapelens {
  /// Class for moment-based lensing and deconvolution with an 
  /// circular weight function.
  /// This is a simplification to the approach of DEIMOSElliptical in that
  /// it uses a circular instead of an elliptcal weight function. 
  /// The optimization for centroid and scale are still done.
  ///
  /// \b NOTE: If Config::USE_WCS is set to \p true, the moments
  /// and all moment-related quantities are
  /// measured in WCS units.
  class DEIMOSCircular : public DEIMOS {
  public:
    /// Default constructor.
    DEIMOSCircular();
    /// Constructor from filename.
    DEIMOSCircular(std::string filename);
    /// Constructor for measuring moments up to order \p N from \p obj.
    /// \p C is the correction order for the deweighting step, 
    /// \p scale the size of the weight function that will be matched to \p obj.
    DEIMOSCircular(const Object& obj, int N, int C, data_t scale);
    /// Constructor for measuring deconvolved moments up to order 
    /// \p N from \p obj.
    /// \p psfs are PSF moment measurements made at different scales.
    /// Only scales available in \p psfs will be used for the matching 
    /// procedure, starting with the largest. The final scale is the one
    /// that maximizes S/N of the measurement.\n
    /// \p C is the correction order for the deweighting step. 
    DEIMOSCircular (const Object& obj, const DEIMOS::PSFMultiScale& psf, int N, int C);
    /// Save to a file.
    void save(std::string filename) const;
    /// Correction order.
    int C;
    /// Matching scale of the weight function.
    data_t matching_scale;
    /// Actual width of the weight function [WCS units].
    data_t scale;
    /// Pixel scale [WCS unit/pixel].
    data_t scale_factor;
    /// Centroid of weight function.
    Point<data_t> centroid;
    /// S/N of moment measurement.
    /// This map stores the S/N of the last measurement for any value of
    /// DEIMOSElliptical::matching_scale employed during the matching procedure.
    std::map<data_t, data_t> SN;
    /// Matching and deconvolution flags.
    /// The flags are ordered as DEIMOS processes objects, 
    /// i.e. matching (first centroid then ellipticity) followed by
    /// deconvolution (if necessary):
    /// - <tt>flags[0]</tt>: centroid determination failed
    /// - <tt>flags[1]</tt>: deconvolution results in non-sensical moments
    std::bitset<2> flags;
    /// Whether the centroid is kept fixed during matching.
    static bool FIX_CENTROID;

  protected:
    void match(const Object& obj, Moments& mo_w);
    void deweight(const Moments& mo_w);
    void computeDeweightingMatrix(const Moments& mo_w);
    void computeCovariances(const Moments& mo_w);
    data_t computeSN(const Moments& mo_w);
    tmv::Matrix<data_t> D;
    History history;
  };
} // end namespace
#endif
