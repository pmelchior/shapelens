#include "../include/DEIMOSForward.h"
#include "../include/LensHelper.h"
#include "../include/MathHelper.h"

namespace shapelens {

  bool DEIMOSForward::FIX_CENTROID = false;
  data_t DEIMOSForward::ETA_MAX = 0.2;

  /*DEIMOSForward::DEIMOSForward(const MultiExposureObject& meo_, const std::vector<DEIMOS::PSFMultiScale>& mePSFMultiScale_, int N, int C) :
    DEIMOS(N), meo(meo_), mePSFMultiScale(mePSFMultiScale_), K(meo_.size()) {
    initialize(C);
    minimize();
    }*/

  DEIMOSForward::DEIMOSForward(const MultiExposureObject& meo_, const MultiExposureObject& mepsf_, int N, int C, const std::set<data_t>& scales) :
    DEIMOS(N), meo(meo_), mepsf(mepsf_), K(meo_.size()) {
    initialize(scales, C);
    minimize(scales);
  }

  void DEIMOSForward::initialize(const std::set<data_t>& scales, int C) { 
    // initialize moments as unit flux point source
    // -> first iteration: object looks like PSF
    mo(0,0) = 1;

    Moments tmp(N);
    tmv::Matrix<data_t> P(mo.size(), mo.size(), 0);
    DEIMOSElliptical d(DEIMOS::N, C);
    centroid(0) = centroid(1) = 0;
    for (int k = 0; k < K; k++) {
      mem.push_back(tmp);
      mem_psf.push_back(tmp);
      meP.push_back(P);

      // set up deimos for each exposure
      d.matching_scale = *scales.begin();
      d.setNoiseImage(meo[k]);
    
      // set scale factors and centroids for WCS exposures
      d.scale_factor = meo[k].grid.getScaleFactor();
      d.centroid = meo[k].centroid;
      if (Config::USE_WCS)
	meo[k].grid.getWCS().transform(d.centroid);
      centroid += d.centroid;
      meD.push_back(d);
    }
    // initialize with avg. centroid
    centroid(0) /= K;
    centroid(1) /= K;
  }
  
  void DEIMOSForward::rollBack(data_t scale, const Moments& best_mo, const Point<data_t>& best_centroid, const std::vector<DEIMOSElliptical>& best_meD, const tmv::Matrix<data_t>& best_S, const MultiExposureMoments& best_mem_psf) {
    matching_scale = scale;
    for (int k = 0; k < K; k++)
      meD[k] = const_cast<DEIMOSElliptical&>(best_meD[k]);
    mo = best_mo;
    centroid = best_centroid;
    S = best_S;
    mem_psf = best_mem_psf;
    flags[1] = 0; // clear flag since previous run was good
    history << " Reverting to previous measurement." << std::endl;
  }

  void DEIMOSForward::minimize(const std::set<data_t>& scales) {
    for (std::set<data_t>::const_iterator siter = scales.begin(); siter != scales.end(); siter++) {
      matching_scale = *siter;
      std::vector<DEIMOSElliptical> best_meD = meD;
      Moments best_mo = mo;
      Point<data_t> best_centroid = centroid;
      tmv::Matrix<data_t> best_S = S;
      MultiExposureMoments best_mem_psf = mem_psf;
      minimizeAtFixedScale();
      if (flags[1] == 0) {
	std::complex<data_t> eps = shapelens::epsilon(mo);
	SN[matching_scale] = mo(0,0)/sqrt(S(0,0));
	eta[matching_scale] = computeEta();
	history << "# Moments at matching_scale = " << matching_scale << ":" << std::endl;
	history << "# " << mo << "\t" << eps << "\t" << shapelens::scale(mo) << "\t" << SN[matching_scale] << "\t" << eta[matching_scale] << std::endl;
	if (eta[matching_scale] > ETA_MAX) {
	  history << "# ETA_MAX exceeded.";
	  if (siter != scales.begin()) {
	    siter--;
	    rollBack(*siter, best_mo, best_centroid, best_meD, best_S, best_mem_psf);
	  } else {
	    history << " Stopping here." << std::endl;
	  }
	  break;
	}
	// save moments in case a later attempt fails
	// assumes that the optimum measurement is the one with the largest
	// eta < ETA_MAX, and that eta grows with matching_scale
	for (int k = 0; k < K; k++)
	  best_meD[k] = meD[k];
	best_mo = mo;
	best_centroid = centroid;
	best_S = S;
	best_mem_psf = mem_psf;
      } else {
	history << "# Minimization failed.";
	if (siter != scales.begin()) {
	  siter--;
	  rollBack(*siter, best_mo, best_centroid, best_meD, best_S, best_mem_psf);
	} else {
	  history << " Aborting." << std::endl;
	}
	break;
      }
    }
  }

  void DEIMOSForward::minimizeAtFixedScale() {
    t = 1; // iteration counter
    Moments mo_save = mo;
    history << "# Modelling object id = " << meo[0].id << " from K = " << K << " exposures at scale = " << matching_scale << std::endl;
    history << "# iter\tmoment guess" << std::endl;
    history << "# exp\tconvolved moments\tscale\teps" << std::endl;
    history << "# ->\tmeasured moments" << std::endl;
    history << "# ==>\tnew moment guess\tflags\tchi^2\tdiff" << std::endl;
    history << "# " << std::string(70, '-') << std::endl;
    
    // compute PSF moments at given scale
    for (int k = 0; k < K; k++) {
      DEIMOSElliptical dp(mepsf[k], N, meD[k].C, matching_scale);
      mem_psf[k] = dp.mo;
    }

    while (true) {
      history << "# " << t << "\t" << mo << std::endl;
      computeMomentsFromGuess();

      // compute chi^2 and best-fit moments
      data_t chi2 = 0;
      unsigned long n_pix = 0;
      tmv::Vector<data_t> diff(mo.size());
      mo.setZero();
      S.setZero();
      for (int k = 0; k < K; k++) {
	DEIMOSElliptical& d = meD[k];
	diff = mem[k];
	diff -= d.mo;
	tmv::Matrix<data_t> S_1 = d.S.inverse();
	data_t chi2_k = diff * (S_1 * (tmv::Vector<data_t>) diff);
	chi2 += chi2_k;
	unsigned long n_pix_k = meo[k].size();
	n_pix += n_pix_k;
	tmv::Matrix<data_t> X = meP[k].transpose() * S_1;
	tmv::Vector<data_t> mo_k = X * (tmv::Vector<data_t>) d.mo;
	mo += mo_k;
	tmv::Matrix<data_t> S_k = X*meP[k];
	S += S_k;
      }
      S = S.inverse();
      mo = S * mo;
      //S /= (data_t) K;
            
      // non-sensical ellipticity check
      flags[1] = flagMoments(mo);
      history << "# ==>\t" << mo << "\t" << flags;
      if (flags[1] == 0) {
	// convergence test:
	// relative changes of the moment vector better than what is
	// needed for rel. epsilon errors of 1e-3
	tmv::Vector<data_t> diff = mo - mo_save;
	data_t ddiff = diff*diff/(mo*mo);
	history << "\t" << chi2 << "\t" << ddiff << std::endl;
	if (ddiff < 1e-8)
	  break;
      
	// update centroid
	if (FIX_CENTROID == false) {
	  Point<data_t> centroid_shift;
	  centroid_shift(0) = mo(1,0)/mo(0,0);
	  centroid_shift(1) = mo(0,1)/mo(0,0);
	  // shift centroids
	  centroid += centroid_shift;
	  for (int k = 0; k < K; k++)
	    meD[k].centroid = centroid; // always shift centroid together
	}
	mo_save = mo;
      } 
      else { // ran into deep trouble here: aborting!
	history << "!" << std::endl;
	break;
      }
      // no convergence reached in due time: set flag 0
      // this is often not a bad object, but with more jitter during
      // the iterations
      if (t==20) {
	flags[0] = 1;
	break;
      }
      t++;
    }
  }

  void DEIMOSForward::computeMomentsFromGuess() {
    // 1) convolve guess with psf moments
    // 2) measure deweighted moments from each exposure
    //    with weight function shape based on guess of convolved moments
    // 3) get the convolved moment errors
    // 4) compute the contribution to chi^2
    for (int k = 0; k < K; k++) {
      convolveExposure(k);
      Moments& mo_c = mem[k];
      DEIMOSElliptical& d = meD[k];
      // set ellipticities and sizes for weight functions in each exposure
      d.matching_scale = matching_scale;
      d.eps = shapelens::epsilon(mo_c);
      d.scale = d.getEpsScale();
      history << "# " << t << "." << k << "\t" << mo_c << "\t" << d.matching_scale << "\t" << d.scale/d.scale_factor << "\t" << d.eps << std::endl;

      DEIMOSElliptical::WeightFunction w(d.scale, d.centroid, d.eps);
      Moments mo_w(meo[k], w, DEIMOS::N + d.C, &d.centroid);
      d.deweight(mo_w);
      d.computeCovariances(mo_w);
      history << "# ->\t" << d.mo << std::endl;
    }
  }

  void DEIMOSForward::convolveExposure(unsigned int k) {
    // create matrix representation of convolution eq. 9
    // a copy from DEIMOS.cc, but this one saves matrix P
    // and applies to convolved moments mem[k]
    DEIMOSElliptical& d = meD[k];
    const Moments& p = mem_psf[k];
    tmv::Matrix<data_t>& P = meP[k];
    for (int n = 0; n <= mo.getOrder(); n++) {
      for (int i = 0; i <= n; i++) {
        int  j = n-i;
	int n = mo.getIndex(i,j);  // convolved moment index
        for (int k = 0; k <= i; k++) {
          for (int l = 0; l <= j; l++) {
            int m = mo.getIndex(k,l); // unconvolved moment index
            P(n,m) = binomial(i,k) * binomial(j,l) * p(i-k,j-l);
          }
        }
      }
    }
    // convolve moments
    mem[k] = P*mo;
    d.R2 = shapelens::R2(mem[k], p);
  }
  
  data_t DEIMOSForward::computeEta() const {
    tmv::Matrix<data_t> A(mo.size(), mo.size());
    A.setZero();
    Moments mo_(mo.getOrder());
    // calculated moments in coordintate frame with
    // eps = eps1 >= 0 (semi-major axis  = x axis)
    // see Teague (1980), eq. 11
    std::complex<data_t> eps_mo = shapelens::epsilon(mo);
    data_t phi = 0.5*atan2(imag(eps_mo),real(eps_mo));
    int i = mo.getIndex(2,0);
    int j = mo.getIndex(2,0);
    A(i,j) = (1+cos(2*phi))/2;
    j = mo.getIndex(1,1);
    A(i,j) = sin(2*phi);
    j = mo.getIndex(0,2);
    A(i,j) = (1-cos(2*phi))/2;
    i = mo.getIndex(1,1);
    j = mo.getIndex(2,0);
    A(i,j) = -sin(2*phi)/2;
    j = mo.getIndex(1,1);
    A(i,j) = cos(2*phi);
    j = mo.getIndex(0,2);
    A(i,j) = sin(2*phi)/2;
    i = mo.getIndex(0,2);
    j = mo.getIndex(2,0);
    A(i,j) = (1-cos(2*phi))/2;
    j = mo.getIndex(1,1);
    A(i,j) = -sin(2*phi);
    j = mo.getIndex(0,2);
    A(i,j) = (1+cos(2*phi))/2;
    mo_ = A*mo;
    // transform moment errors, too
    tmv::Matrix<data_t> S_ = A * S * A.transpose();
    i = mo.getIndex(2,0);
    data_t dQ11_ = sqrt(S_(i,i))/mo_(i);
    i = mo.getIndex(0,2);
    data_t dQ22_ = sqrt(S_(i,i))/mo_(i);
    // relative moment error in rotated frame:
    // independent of orientation!
    return 2*dQ11_*dQ22_/(dQ11_ + dQ22_);
  }


  const MultiExposureMoments& DEIMOSForward::getPSFMoments() const {
    return mem_psf;
  }
}
