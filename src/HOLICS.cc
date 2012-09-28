#include "../include/HOLICS.h"
#include "../include/WeightFunction.h"
#include "../include/MathHelper.h"

namespace shapelens {
  HOLICS::HOLICS(const Object& obj, data_t scale) {
    GaussianWeightFunction w(scale, obj.centroid);

    sigma = w.getScale();
    
    Moments mo(obj, w, 4);
    M = mo(0,0);
    trQ = __trQ(mo);
    xi = __xi(mo);
    zeta = __zeta(mo);
    delta = __delta(mo);
    
    // measure moments with W'
    w.setDerivative(1);
    Moments mo_(obj, w, 6);
    M_ = mo_(0);
    trQ_ = __trQ(mo_);
    xi_ = __xi(mo_);
    v0_ = __v0(mo_);
    
    // measure moments with W''
    w.setDerivative(2);
    Moments mo__(obj, w, 4);
    trQ__ = __trQ(mo__);
    xi__ = __xi(mo__);

    // measure moments with W'''
    w.setDerivative(3);
    Moments mo___(obj, w, 6);
    xi___ = __xi(mo___);
    v0___ = __v0(mo___);

    // form high-order quantities...
    // eq. (C27)
    C_0_zeta = 9./4 + 3./(4*sigma*sigma)*v0_/xi;
    C_Delta_zeta = -2*trQ/xi - xi_/(sigma*sigma*xi);
    C_0_delta = 3./4 + v0_/(4*sigma*sigma*xi);
    // eq. (B8)
    Delta_0_L = 1.5 * trQ/M + (3*xi_)/(4*sigma*sigma*M);
    Delta_0_L /= (1 + trQ_/(sigma*sigma*M));
    // eq. (C18)
    P_0_zeta =  1./xi *( M + 5*trQ_/(sigma*sigma) + 7*xi__/(2*pow4(sigma)) + v0___/(2*pow6(sigma)));
    P_0_D = M_/(M*sigma*sigma) + 2*trQ__/(M*pow4(sigma)) + xi___/(2*M*pow6(sigma));
    P_Delta_D = 1 + trQ_/(M*sigma*sigma);
    P_Delta_zeta = 2*trQ/xi + xi_/(xi*sigma*sigma);
    P_0_delta = 1./xi *( M + 3*trQ_/(sigma*sigma) + 3*xi__/(2*pow4(sigma)) + v0___/(6*pow6(sigma)));
  }

  // eq. (20)
  complex<data_t> HOLICS::__zeta(const Moments& mo) const {
    return complex<data_t>(mo(3,0) + mo(1,2),
			   mo(2,1) + mo(0,3))/__xi(mo);
  }
  // eq. (20)
  complex<data_t> HOLICS::__delta(const Moments& mo) const {
    return complex<data_t>(mo(3,0) - 3*mo(1,2), 
			   3*mo(2,1) - mo(0,3))/__xi(mo);
  }
  // trace of eq. (12)
  data_t HOLICS::__trQ(const Moments& mo) const {
    return mo(2,0) + mo(0,2);
  }
  // eq. (21)
  data_t HOLICS::__xi(const Moments& mo) const {
    return mo(4,0) + 2*mo(2,2) + mo(0,4);
  }
  // eq. (A6)
  data_t HOLICS::__v0(const Moments& mo) const {
    return mo(6,0) + 3*mo(2,4) + 3*mo(4,2) + mo(0,6);
  }

  
  // eq. (34)
  complex<data_t> HOLICS::F() const {
    return zeta/(9./4 + 3*v0_/(4*xi*sigma*sigma) - 
		 (2*trQ/xi + xi_/(xi*sigma*sigma))*Delta_0_L);
  }
  // eq. (35)
  complex<data_t> HOLICS::G() const {
    return delta/(3./4 + v0_/(4*xi*sigma*sigma));
  }
  // eqs. (52) & (78)
  complex<data_t> HOLICS::F(const HOLICS& s) const {
    const HOLICS& g = *this;
    complex<data_t> zeta_iso = g.zeta - 
      (g.P_0_zeta - s.P_0_D/s.P_Delta_D*s.P_Delta_zeta)*__zeta_q(s);
    return zeta_iso/((g.C_0_zeta + g.C_Delta_zeta*g.Delta_0_L) -
		     (s.C_0_zeta + s.C_Delta_zeta*s.Delta_0_L)*
		     (s.P_0_zeta - s.P_0_D*s.P_Delta_zeta/s.P_Delta_D)*
		     (g.P_0_zeta - g.P_0_D*g.P_Delta_zeta/g.P_Delta_D));
  }
  // eqs. (53) & (79)
  complex<data_t> HOLICS::G(const HOLICS& s) const {
    const HOLICS& g = *this;
    complex<data_t> delta_iso = g.delta - g.P_0_delta*__delta_q(s);
    return delta_iso/(g.C_0_delta - 
		      (s.C_0_delta/s.P_0_delta)*g.P_0_delta);
  }
  // eq. (56)
  complex<data_t> HOLICS::__zeta_q(const HOLICS& s) const {
    return s.zeta/(s.P_0_zeta - 
		   s.P_0_D/s.P_Delta_D*s.P_Delta_zeta);
  }
  // eq. (57)
  complex<data_t> HOLICS::__delta_q(const HOLICS& s) const {
    return s.delta/s.P_0_delta;
  }

} // end namespace
