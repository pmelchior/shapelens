#include "../include/LensHelper.h"
#include "../include/MathHelper.h"

// All equations taken from Bartelmann & Schneider (2001)
namespace shapelens {
  
  void eps2chi(const std::complex<data_t>& eps, std::complex<data_t>& chi) {
    chi = 2.*eps/(1+pow2(abs(eps)));
  }
  void chi2eps(const std::complex<data_t>& chi, std::complex<data_t>& eps) {
    eps = chi/(1. + sqrt(complex<data_t>(1 - pow2(abs(chi)))));
  }
  void lensEps(const std::complex<data_t>& gamma, std::complex<data_t>& eps) {
    eps = (eps + gamma)/(1. + conj(gamma)*eps);
  }

  std::complex<data_t> epsilon(const Moments& mo) {
    if (mo.getOrder() >= 2) {
      if (mo(2,0) == 0. && mo(1,1) == 0. && mo(0,2) == 0.)
	return std::complex<data_t>(0,0);
      else {
	std::complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
	e/= (std::complex<data_t>(mo(2,0) + mo(0,2)) + 2.*sqrt(std::complex<data_t>(mo(0,2)*mo(2,0) - mo(1,1)*mo(1,1))));
	return e;
      }
    } else
      return std::complex<data_t>(0,0);
  }

  std::complex<data_t> chi(const Moments& mo) {
    if (mo.getOrder() >= 2) {
      if (mo(2,0) == 0. && mo(1,1) == 0. && mo(0,2) == 0.)
	return std::complex<data_t>(0,0);
      else {
	std::complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
	e/= mo(2,0) + mo(0,2);
	return e;
      }
    } else
      return std::complex<data_t>(0,0);
  }

  // HOLICS equations from OMU 2007
  complex<data_t> zeta(const Moments& mo) {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4); // eq. 25
      complex<data_t> zeta(mo(3,0) + mo(1,2),    // eq. 26
			   mo(2,1) + mo(0,3));
      zeta /= xi;
      return zeta;
    } else
      return complex<data_t>(0,0);
  }

  complex<data_t> delta(const Moments& mo) {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4);  // eq. 25
      complex<data_t> delta(mo(3,0) - 3*mo(1,2),  // eq. 27
			    3*mo(2,1) - mo(0,3));
      delta /= xi;
      return delta;
    } else
      return complex<data_t>(0,0);
  }

  data_t scale(const Moments& mo) {
    if (mo.getOrder() >= 2) {
      if ((mo(2,0)+mo(0,2))/mo(0,0) > 0)
	return sqrt((mo(2,0)+mo(0,2))/mo(0,0));
      else
	return 0;
    }
    else
      return 0;
  }
  
  data_t R2(const Moments& g, const Moments& p) {
    if (g.getOrder() >= 2) {
      data_t scale_g = scale(g);
      if (scale_g > 0)
	return 1 - pow2(scale(p)/scale(g));
      else
	return 0;
    }
    else
      return 0;
  } 

} // end namespace
