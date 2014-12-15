#ifndef SHAPELENS_LENSHELPER
#define SHAPELENS_LENSHELPER
#include <complex>
#include "Moments.h"
#include "MathHelper.h"

namespace shapelens {
  
  /// Convert \f$\epsilon\f$ to \f$\chi\f$.
  inline void eps2chi(const std::complex<data_t>& eps, std::complex<data_t>& chi) {
    chi = 2.*eps/(1+pow2(abs(eps)));
  }
  /// Convert \f$\chi\f$ to \f$\epsilon\f$.
  inline void chi2eps(const std::complex<data_t>& chi, std::complex<data_t>& eps) {
     eps = chi/(1. + sqrt(complex<data_t>(1 - pow2(abs(chi)))));
  }
  /// Apply shear \f$\gamma\f$ to unlensed ellipticity \f$\epsilon\f$.
  inline void lensEps(const std::complex<data_t>& gamma, std::complex<data_t>& eps) {
    eps = (eps + gamma)/(1. + conj(gamma)*eps);
  }
  /// Get tangential projection of ellipticity or shear.
  /// Assumes right-handed coordinate system for \p eps and
  /// tangential projection of spherical center coordinates.
  inline data_t epsTangential(const std::complex<data_t>& eps, const Point<data_t>& center, const Point<data_t>& center_lens) {
    data_t phi = atan2(center(1)-center_lens(1), center(0) - center_lens(0));
    return -real(eps)*cos(2*phi) - imag(eps)*sin(2*phi);
  }
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
