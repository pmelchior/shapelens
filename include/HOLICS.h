#ifndef SHAPELENS_HOLICS_H
#define SHAPELENS_HOLICS_H

#include "Object.h"
#include "Moments.h"

namespace shapelens {
  class HOLICS {
  public :
    HOLICS(const Object& obj, data_t scale);
    std::complex<data_t> zeta;
    std::complex<data_t> delta;
    std::complex<data_t> F() const;
    std::complex<data_t> G() const;
    std::complex<data_t> F(const HOLICS& psf) const;
    std::complex<data_t> G(const HOLICS& psf) const;
    
  protected:
    data_t __trQ(const Moments& mo) const;
    data_t __xi(const Moments& mo) const;
    data_t __v0(const Moments& mo) const;
    std::complex<data_t> __zeta(const Moments& mo) const;
    std::complex<data_t> __delta(const Moments& mo) const;
    std::complex<data_t> __zeta_q(const HOLICS& star) const;
    std::complex<data_t> __delta_q(const HOLICS& star) const;

    data_t sigma, 
      trQ, trQ_, trQ__, 
      xi, xi_, xi__, xi___, 
      M, M_,
      v0_, v0___,
      C_0_zeta, C_Delta_zeta, C_0_delta,
      Delta_0_L,
      P_0_zeta, P_0_D, P_Delta_D, P_Delta_zeta, P_0_delta;
  };

} // end namespace shapelens

#endif
