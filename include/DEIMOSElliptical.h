#ifndef SHAPELENS_DEIMOSELLIPTICAL_H
#define SHAPELENS_DEIMOSELLIPTICAL_H

#include "DEIMOS.h"
#include "History.h"

namespace shapelens {
  /// Class for moment-based lensing and deconvolution with an 
  /// elliptical weight function.
  /// This approach is published in 
  /// <i>Melchior et al., MNRAS, 412, 1552 (2011)</i>.
  /// It extends DEIMOS by incorporating moment measurement with a matched
  /// weight function and an appropriate deweighting correction.
  /// It uses an iterative matching algorithm to determine the optimal
  /// size, center, and ellipticity of a weight function to measure the
  /// convolved object. The objective of the matching step is to maximum
  /// the S/N of the moment measurements.
  ///
  /// For galaxies, a DEIMOS deconvolution will be performed with
  /// a predefined PSF shape. To allow for size adjustments during the matching
  /// procedure, measurements of the PSF moments at several should be provided.
  /// If any problems occur during the matching or deconvolution steps, it
  /// is recorded in DEIMOSElliptical::flags.
  ///
  /// A typical example for obtaining deconvolved galaxy moments is given below.
  /// \p scales is an ordered list (aka \p set) of scales that should cover
  /// the range of sizes that allow an accurate measurement of the PSF moments.
  /// \code
  /// int C = 4; // weighting correction order
  /// DEIMOS::PSFMultiScale psfs;
  /// for (std::set<data_t>::iterator iter = scales.begin(); iter != scales.end(); iter++) {
  ///   DEIMOSElliptical p(psf, 2, C, *iter);
  ///   data_t flux = p.mo(0,0);
  ///   p.mo /= flux;
  ///   psfs.insert(*iter, p.mo);
  /// }
  /// DEIMOSElliptical d(obj, psfs, 2, C);
  /// \endcode
  /// 
  /// \b NOTE: If Config::USE_WCS is set to \p true, the moments mo
  /// and all moment-related quantities (foremost scale, centroid, epsilon) are
  /// measured in WCS units.
  class DEIMOSElliptical : public DEIMOS {
  public:
    /// An elliptical Gaussian weight function.
    /// The functional form of the weight function is given by
    /// \f$W(x) = \exp\Bigl(\frac{-x^{\prime 2}}{2s^s}\Bigr)\f$,
    /// where \f$x^\prime\f$ undergoes a suitable CoordinateTransformation to
    /// create the ellipticity.
    class WeightFunction : public GaussianWeightFunction {
    public:
      /// Constructor for ellipticial weight function.
      WeightFunction(data_t scale, const Point<data_t>& centroid, const std::complex<data_t>& eps);
      /// Get value of weight function at position \p P.
      virtual data_t operator()(const Point<data_t>& P) const;
    protected:
      LensingTransformation T;
      const Point<data_t>& centroid;
    };

    /// Default constructor.
    DEIMOSElliptical(int N=0, int C=0);
    /// Constructor to read in \p filename.
    /// \p filename should have been created by save().
    DEIMOSElliptical(std::string filename);
    /// Constructor for measuring moments up to order \p N from \p obj.
    /// \p C is the correction order for the deweighting step, 
    /// \p scale the size of the weight function that will be matched to \p obj.
    DEIMOSElliptical(const Object& obj, int N, int C, data_t scale);
    /// Constructor for measuring deconvolved moments up to order 
    /// \p N from \p obj.
    /// \p psfs are PSF moment measurements made at different scales.
    /// Only scales available in \p psfs will be used for the matching 
    /// procedure, starting with the largest. The final scale is the one
    /// that maximizes S/N of the measurement.\n
    /// \p C is the correction order for the deweighting step. 
    DEIMOSElliptical (const Object& obj, const DEIMOS::PSFMultiScale& psfs, int N, int C);
    /// Save moments and auxiliary information to \p filename.
    void save(std::string filename) const;
    /// Correction order.
    /// To correct for the employment of a weight function, higher-order moments
    /// are taken into account. The moment measurement thus comprises all
    /// moments up to order <tt>N+C</tt>.
    int C;
    /// Final matching scale of the weight function.
    /// Corresponds to the size of the semi-minor axis [pixel]:
    data_t matching_scale;
    /// Actual width of the weight function [WCS units].
    /// Determined during the matching procedure as
    /// \f$s = s_m/\sqrt{1 + \epsilon^2 - 2|\epsilon|}\f$,
    /// where \f$s_m\f$ denotes DEIMOSElliptical::matching_scale
    /// and \f$\epsilon\f$ the convolved ellipticity DEIMOSElliptical::eps.
    /// If WCS units are used, \f$s\f$ is rescaled accordingly. 
    data_t scale;
    /// Pixel scale [WCS unit/pixel].
    data_t scale_factor;
    /// Centroid of weight function.
    Point<data_t> centroid;
    /// Ellipticity of weight function.
    std::complex<data_t> eps;
    /// S/N of moment measurement.
    /// This map stores the S/N of the last measurement for any value of
    /// DEIMOSElliptical::matching_scale employed during the matching procedure.
    std::map<data_t, data_t> SN;
    /// Matching and deconvolution flags.
    /// The flags are ordered as DEIMOS processes objects, 
    /// i.e. matching (first centroid then ellipticity) followed by
    /// deconvolution (if necessary):
    /// - <tt>flags[0]</tt>: centroid determination failed
    /// - <tt>flags[1]</tt>: non-sensical moments or ellipiticty
    /// - <tt>flags[2]</tt>: ellipticity matching failed
    /// - <tt>flags[3]</tt>: deconvolution results in non-sensical moments
    std::bitset<4> flags;
    /// Whether the centroid is kept fixed during matching.
    static bool FIX_CENTROID;

    friend class DEIMOSForward;
  protected:
    void match(const Object& obj, Moments& mo_w);
    void deweight(const Moments& mo_w);
    void computeCovariances(const Moments& mo_w);
    data_t computeSN(const Moments& mo_w);
    data_t getEpsScale() const;
    tmv::Matrix<data_t> D;
    History history;
  };
} // end namespace
#endif
