#ifndef SHAPELENS_KSB_H
#define SHAPELENS_KSB_H

#include "Object.h"
#include "Moments.h"
#include <set>

namespace shapelens {
  /// Class for Kaiser, Squires & Broadhurst shear estimates.
  /// The implementation is a standard-KSB approach as described in
  /// <i>Kaiser et al., ApJ, 449, 460 (1995)</i> with terminology taken
  /// from <i>Hoekstra et al., ApJ, 504, 636 (1998)</i>. 
  /// No averaging over KSB::P_sh or KSB::P_sm 
  /// is employed, all galaxies are treated independently. In addition, two
  /// modifications to the original shear estimator KSB::gamma() are provided:
  /// - In KSB::gammaTr(), KSB::P_sm is approximated 
  ///   by 1/2 of its trace.
  /// - KSB::gamma1() and KSB::gammaTr1() are strictly first-order shear 
  ///   estimators of the original and the trace trick estimator.
  ///
  /// See <i>Viola et al., MNRAS, 410, 2156 (2011)</i> for details.
  ///
  /// The argumented constructors take an object and a measurement scale;
  /// alternatively a set of scales can be provided, the final measurement
  /// will then be performed with the scale that maximizes S/N.
  ///
  /// \b NOTE: If Config::USE_WCS is set to \p true, all moments measurements
  /// and thus all derived quantities (foremost KSB::chi) are
  /// measured in WCS units.
  class KSB {
  public:
    /// Constructor.
    KSB();
    /// Constructor for shear estimates from \p obj.
    /// The measurement employs a circular Gaussian weight function with
    /// a width given by \p scale.\n
    /// If KSB::FIX_CENTROID is \p false (default),
    /// the centroid will be iteratively adjusted to have zero dipole moment
    /// under the weight funtion.
    KSB(const Object& obj, data_t scale);
    /// Constructor for optimized shear estimates from \p obj.
    /// Centroid position and optimum weight function scale are determined
    /// iteratively by requiring that the dipole moment of the weighted
    /// brightness distribution vanishes and that the measurement S/N is
    /// maximized.\n
    /// If KSB::FIX_CENTROID is \p false (default),
    /// the centroid will be iteratively adjusted to have zero dipole moment
    /// under the weight funtion.
    KSB(const Object& obj, const std::set<data_t>& scales);
    /// Resulting complex ellipticity \f$\chi\f$.
    std::complex<data_t> chi;
    /// S/N of the measurement.
    data_t SN;
    /// Centroid position for the measurement.
    Point<data_t> centroid;
    /// Scale for the measurement.
    data_t scale;
    /// Shear tensor.
    tmv::Matrix<data_t> P_sh;
    /// First-order shear tensor.
    tmv::Matrix<data_t> P1;
    /// Smear tensor.
    tmv::Matrix<data_t> P_sm;

    // shear estimators without PSF correction
    // for stellar shapes
    /// Original KSB shear estimator.
    std::complex<data_t> gamma() const;
    /// First-order shear estimator. 
    std::complex<data_t> gamma1() const;
    /// Shear estimator with the \i trace-trick.
    std::complex<data_t> gammaTr() const;
    /// First-order shear estimator with the \i trace-trick.
    std::complex<data_t> gammaTr1() const;

    // shear estimators with PSF correction: for galaxies
    /// Original KSB shear estimator with PSF correction.
    std::complex<data_t> gamma(const KSB& psf) const;
    /// First-order shear estimator with PSF correction.
    std::complex<data_t> gamma1(const KSB& psf) const;
    /// Shear estimator with the \i trace-trick and PSF correction.
    std::complex<data_t> gammaTr(const KSB& psf) const;
    /// First-order shear estimator with the \i trace-trick and PSF correction.
    std::complex<data_t> gammaTr1(const KSB& psf) const;

    /// Whether the centroid is kept fixed for the measurement.
    static bool FIX_CENTROID;

  protected:
    data_t __trQ(const Moments& mo) const;
    data_t __psi(const Moments& mo) const;
    data_t __mu(const Moments& mo) const;
    data_t __nu(const Moments& mo) const;
    data_t __pi(const Moments& mo) const;
    std::complex<data_t> __p(const KSB& star) const;
    void findCentroid(const Object& obj, GaussianWeightFunction& w);
    void measureMoments(const Object& obj, GaussianWeightFunction& w);
    data_t scale_factor;
  };

} // end namespace shapelens

#endif
