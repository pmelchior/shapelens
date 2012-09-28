#ifndef SHAPELENS_LENSHELPER
#define SHAPELENS_LENSHELPER
#include <complex>
#include "Moments.h"

namespace shapelens {
  
  /// Convert \f$\epsilon\f$ to \f$\chi\f$.
  void eps2chi(const std::complex<data_t>& eps, std::complex<data_t>& chi);
  /// Convert \f$\chi\f$ to \f$\epsilon\f$.
  void chi2eps(const std::complex<data_t>& chi, std::complex<data_t>& eps);
  /// Apply shear \f$\gamma\f$ to unlensed ellipticity \f$\epsilon\f$.
  void lensEps(const std::complex<data_t>& gamma, std::complex<data_t>& eps);
  /// Get complex ellipticity \f$\epsilon\f$ from Moments.
  std::complex<data_t> epsilon(const Moments& mo);
  /// Get complex ellipticity \f$\chi\f$ from Moments.
  std::complex<data_t> chi(const Moments& mo);
  /// Get equivalent scale of Gaussian from Moments.
  /// Returns \f$\sqrt{(Q_{11} + Q_{22})/F}\f$, 
  /// i.e. the flux-normalized sum of the
  /// symmetric second-order moments \f$Q_{ii}\f$.
  data_t scale(const Moments& mo);
  /// Get first flexion distortion \f$\zeta\f$ from Moments.
  std::complex<data_t> zeta(const Moments& mo);
  /// Get second flexion distortion \f$\delta\f$ from Moments.
  std::complex<data_t> delta(const Moments& mo);
  /// Resolution factor \f$R_2\f$.
  /// From Hirata et al. (2004), eq. 8: 
  /// \f$R_2 = 1 - [s^2(p) / s^2(g^*)]\f$, where \f$s(.)\f$
  /// denotes the shapelens::scale of either the PSF or the convolved galaxy.
  data_t R2(const Moments& g_star, const Moments& p);
} // end namespace

#endif
