#include "../include/KSB.h"
#include "../include/WeightFunction.h"
#include "../include/LensHelper.h"

namespace shapelens {
  
  bool KSB::FIX_CENTROID = false;
  
  KSB::KSB() {}
  
  KSB::KSB(const Object& obj, data_t scale_) :
    centroid(obj.centroid), scale(scale_) {
    scale_factor = obj.grid.getScaleFactor();
    scale *= scale_factor;
    if (Config::USE_WCS)
      obj.grid.getWCS().transform(centroid);
    
    GaussianWeightFunction w(scale,centroid);
    measureMoments(obj, w);
  }
 
  KSB::KSB(const Object& obj, const std::set<data_t>& scales) :
    centroid(obj.centroid), scale(*scales.begin()) {
    scale_factor = obj.grid.getScaleFactor();
    scale *= scale_factor;
    if (Config::USE_WCS)
      obj.grid.getWCS().transform(centroid);
    
    GaussianWeightFunction w(scale,centroid);
    measureMoments(obj, w);
    std::set<data_t>::const_iterator iter = scales.begin();

    if (iter != scales.end()) {
      KSB best = *this; 
      iter++;
      for (iter; iter != scales.end(); iter++) {
	scale = *iter * scale_factor;
	w = GaussianWeightFunction (scale,centroid);
	measureMoments(obj, w);
	if (SN > best.SN) {
	  best = *this;
	}
      }
      *this = best;
    }
  }

  void KSB::measureMoments(const Object& obj, GaussianWeightFunction& w) {
    
    // find centroid under w
    if (KSB::FIX_CENTROID == false)
      findCentroid(obj, w);
  
    // various moment combinations measured with different derivatives 
    // of the weight function
    data_t trQ, trQ_, M, mu_, mu__, psi_, psi__, pi_, pi__, nu_, nu__;

    // moments measured with w
    Moments mo(obj, w, 2, &centroid);
    M = mo(0,0);
    trQ = __trQ(mo);
    chi = shapelens::chi(mo);

    // measure moments with W' (w.r.t r^2)
    w.setDerivative(-1);
    Moments mo_(obj, w, 4, &centroid);
    trQ_ = __trQ(mo_);
    psi_ = __psi(mo_);
    mu_ = __mu(mo_);
    pi_ = __pi(mo_);
    nu_ = __nu(mo_);

    // measure moments with W'' (w.r.t r^2)
    w.setDerivative(-2);
    Moments mo__(obj, w, 4, &centroid);
    mu__ = __mu(mo__);
    psi__ = __psi(mo__);
    pi__ = __pi(mo__);
    nu__ = __nu(mo__);

    // measure S/N
    SN = M/(obj.noise_rms * sqrt(M_PI) * w.getScale());

    // Compute P_sh 
    P_sh.resize(2,2);
    P_sh(0,0) = -2*real(chi)*pi_/trQ-2*real(chi)*real(chi)+2*psi_/trQ + 2.0 ;
    P_sh(0,1) = -4*real(chi)*nu_/trQ-2*imag(chi)*real(chi)+4*mu_/trQ;
    P_sh(1,0) = -2*imag(chi)*pi_/trQ-2*imag(chi)*real(chi)+4*mu_/trQ;
    P_sh(1,1) = -4*imag(chi)*nu_/trQ-2*imag(chi)*imag(chi)+8*mo_(2,2)/trQ+2.0;

    // Compute B
    tmv::Matrix<data_t> B(2,2);
    B(0,0)=psi_;
    B(0,1)=2*mu_;
    B(1,0)=B(0,1);
    B(1,1)=4*mo_(2,2);

    // Compute P_sh, strictly I order
    P1.resize(2,2);
    P1(0,0) = 2*B(0,0)/trQ + 2. ;
    P1(0,1) = 2*B(0,1)/trQ;
    P1(1,0) = P1(0,1);
    P1(1,1) = 2*B(1,1)/trQ + 2. ;
        
    // compute P_sm
    P_sm.resize(2,2);
    P_sm(0,0) = M + 2*trQ_ + psi__;
    P_sm(0,1) = 2*mu__;
    P_sm(1,0) = P_sm(0,1);
    P_sm(1,1) = M + 2*trQ_ + 4*mo__(2,2);
    P_sm /= trQ;
    
    std::complex<data_t> e_sm;
    real(e_sm) = (2*(mo_(2,0) - mo_(0,2)) + pi__)/trQ;
    imag(e_sm) = (4*mo_(1,1) + 2*nu__)/trQ;
    P_sm(0,0) -= real(chi)*real(e_sm);
    P_sm(0,1) -= real(chi)*imag(e_sm);
    P_sm(1,0) -= imag(chi)*real(e_sm);
    P_sm(1,1) -= imag(chi)*imag(e_sm);
  }

  void KSB::findCentroid(const Object& obj, GaussianWeightFunction& w) {
    data_t delta;
    do {
      Moments mo(obj, w, 1, &centroid); 
      Point<data_t> shift(mo(1,0)/mo(0,0), mo(0,1)/mo(0,0));
      delta = sqrt(shift(0)*shift(0) + shift(1)*shift(1)) / scale_factor;
      centroid += shift;
      w.setCentroid(centroid);
    } while (delta > 1e-2);
  }

  data_t KSB::__trQ(const Moments& mo) const {
    return mo(2,0) + mo(0,2);
  }
  data_t KSB::__psi(const Moments& mo) const {
    return mo(4,0) - 2*mo(2,2) + mo(0,4);
  }
  data_t KSB::__mu(const Moments& mo) const {
    return mo(3,1) - mo(1,3);
  }
  data_t KSB::__nu(const Moments& mo) const {
    return mo(3,1) + mo(1,3);
  }
  data_t KSB::__pi(const Moments& mo) const {
    return mo(4,0) - mo(0,4);
  }

  // helper function
  complex<data_t> cmulti(const tmv::GenMatrix<data_t>& M, complex<data_t> d) {
    return complex<data_t> (M(0,0)*real(d) + M(0,1)*imag(d),
			    M(1,0)*real(d) + M(1,1)*imag(d));
  }

  // eq. (22)
  complex<data_t> KSB::__p(const KSB& star) const {
    return cmulti(star.P_sm.inverse(), star.chi);
  }
  // eq. (12)
  complex<data_t> KSB::gamma() const {
    return cmulti(P_sh.inverse(), chi);
  }
  // P_sh stricktly to I order
  complex<data_t> KSB::gamma1() const {
    return cmulti(P1.inverse(), chi);
  }
  //using trace P^g
  complex<data_t> KSB::gammaTr() const{
    return (2./(P_sh(0,0)+P_sh(1,1)))*chi;
  }
  //using trace first order
  complex<data_t> KSB::gammaTr1() const{
    return (2./(P1(0,0)+P1(1,1)))*chi;
  }
  // eqs. (23) and (27, no delta)
  complex<data_t> KSB::gamma(const KSB& psf) const {
    tmv::Matrix<data_t> P_gamma = P_sh - psf.P_sh*psf.P_sm.inverse()*P_sm;
    complex<data_t> p = __p(psf);
    return cmulti(P_gamma.inverse(), chi - cmulti(P_sm,p));
  }

  complex<data_t> KSB::gamma1(const KSB& psf) const {
    tmv::Matrix<data_t> P_gamma = P1 - psf.P1*psf.P_sm.inverse()*P_sm;
    complex<data_t> p = __p(psf);
    return cmulti(P_gamma.inverse(), chi - cmulti(P_sm,p));
  }

  complex<data_t> KSB::gammaTr(const KSB& psf) const {
    tmv::Matrix<data_t> P_gamma = P_sh - psf.P_sh*psf.P_sm.inverse()*P_sm;
    complex<data_t> p = __p(psf);
    return (2./(P_gamma(0,0)+P_gamma(1,1)))*(chi - cmulti(P_sm,p));
  }
  
  complex<data_t> KSB::gammaTr1(const KSB& psf) const {
    tmv::Matrix<data_t> P_gamma = P1 - psf.P1*psf.P_sm.inverse()*P_sm;
    complex<data_t> p = __p(psf);
    return (2./(P_gamma(0,0)+P_gamma(1,1)))*(chi - cmulti(P_sm,p));
  }
} // end namespace
