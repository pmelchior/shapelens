#include "../include/DEIMOSForward.h"
#include "../include/LensHelper.h"
#include "../include/MathHelper.h"

namespace shapelens {

  bool DEIMOSForward::FIX_CENTROID = false;

  DEIMOSForward::DEIMOSForward(const MultiExposureObject& meo_, const std::vector<DEIMOS::PSFMultiScale>& mePSFMultiScale_, int N, int C) :
    DEIMOS(N), meo(meo_), mepsf(MultiExposureObject()), mePSFMultiScale(mePSFMultiScale_), K(meo_.size()) {
    initialize(C);
    minimize();
  }

  void DEIMOSForward::initialize(int C) { 
    // initialize moments as unit flux point source
    // -> first iteration: object looks like PSF
    mo(0,0) = 1;

    Moments tmp(N);
    tmv::Matrix<data_t> P(mo.size(), mo.size(), 0);
    DEIMOS::PSFMultiScale psfs;
    DEIMOSElliptical d(DEIMOS::N, C);
    centroid(0) = centroid(1) = 0;
    for (int k = 0; k < K; k++) {
      mem.push_back(tmp);
      meP.push_back(P);

      // set up deimos for each exposure
      d.matching_scale = psfs.getMinimumScale();
      d.setNoiseImage(meo[k]);
    
      // set scale factors and centroids for WCS exposures
      d.scale_factor = meo[k].grid.getScaleFactor();
      d.centroid = meo[k].centroid;
      if (Config::USE_WCS)
	meo[k].grid.getWCS().transform(d.centroid);
      centroid += d.centroid;
      meD.push_back(d);

      meSaveScale.push_back(0);
      meTroubleScale.push_back(1e100);
    }
    // initialize with avg. centroid
    centroid(0) /= K;
    centroid(1) /= K;
  }
  
  void DEIMOSForward::minimize() { 
    data_t last_SN = 0; 
    int t = 1;
    Moments mo_save = mo;
    history << "# Modelling object id = " << meo[0].id << " from K = " << K << " exposures" << std::endl;
    history << "# iter\tmoment guess" << std::endl;
    history << "# exp\tconvolved moments\tscale\teps" << std::endl;
    history << "# ->\tmeasured moments" << std::endl;
    history << "# ==>\tnew moment guess\tSN\teps\tchi^2\teta" << std::endl;
    history << "# " << std::string(70, '-') << std::endl;
    while (true) {
      history << "# " << t << "\t" << mo << std::endl;
      computeMomentsFromGuess(t);

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
      SN = mo(0,0)/sqrt(S(0,0));

      int i = mo.getIndex(2,0);
      data_t eta = sqrt(S(i,i)) / fabs(mo(i));
      i = mo.getIndex(0,2);
      eta += sqrt(S(i,i)) / fabs(mo(i));
      eta /= 2;

      history << "# ==>\t" << mo << "\t" << SN << "\t" << shapelens::epsilon(mo) << "\t" << chi2 << "\t" << eta << std::endl;

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

      // non-sensical ellipticity check
      flags[1] = flagMoments(mo);
      flags[0] = flags[0] | flags[1]; // keep track of difficulties
      if (flags[1] == 0) {
	mo_save = mo;
	for (int k = 0; k < K; k++)
	  meSaveScale[k] = std::max(meD[k].matching_scale, meSaveScale[k]);
      }
      else {
	mo = mo_save;
	for (int k = 0; k < K; k++)
	  meTroubleScale[k] = meD[k].matching_scale;
      }

      t++;

      // convergence: chi^2 does not change by more than relative 1e-2
      // note: chi^2 might increase slightly before convergence
      // so go for absolute value of change
      /*
      if (t==10 || (last_chi2 > 0 && fabs(last_chi2 - chi2) < 1e-2*last_chi2))
	break;
      last_chi2 = chi2;
      */

      // convergence of SN more stable than chi^2
      if (t==10 || (last_SN > 0 && fabs(last_SN - SN) < 1e-3*last_SN))
	break;
      last_SN = SN;
    }
  }

  void DEIMOSForward::computeMomentsFromGuess(int t) {
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
      d.eps = shapelens::epsilon(mo_c);
      d.matching_scale = getWeightFunctionScale(k);
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
    DEIMOSElliptical& d = meD[k];
    DEIMOS::PSFMultiScale& psfs = mePSFMultiScale[k];
    data_t scale = psfs.getScaleClosestTo(d.matching_scale);

    // when PSF images are available:
    // if this if too different, create a new PSF moment measurement
    if (mepsf.size() == K) {
      if (fabs(scale-d.matching_scale) > 0.05*d.matching_scale) {
	scale = d.matching_scale;
	DEIMOSElliptical psf(mepsf[k], DEIMOS::N, d.C, scale);
	psfs.insert(scale, psf.mo);
      }
    }

    // create matrix representation of convolution eq. 9
    // a copy from DEIMOS.cc, but this one saves matrix P
    // and applies to convolved moments mem[k]
    Moments& p = psfs[scale];
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
  }
  
  data_t DEIMOSForward::getWeightFunctionScale(unsigned int k) const {
    const DEIMOS::PSFMultiScale& psfs = mePSFMultiScale[k];
    const DEIMOSElliptical& d = meD[k];
    // simply use sqrt(trQ/F) as size: equivalent width of Gaussian 
    const Moments& m = mem[k];
    data_t s = sqrt((m(0,2) + m(2,0))/m(0,0));
    
    // Noise correction, expect noise variance to be set
    if (S(0,0) > 0) {
      const tmv::Matrix<data_t>& P = meP[k];
      tmv::Matrix<data_t> S_c = P*S*P.transpose();
      // error on sqrt(trQ), ignoring flux uncertainty
      data_t sigma_2 = 0;
      unsigned int i = m.getIndex(2,0);
      sigma_2 += S_c(i,i);
      i = m.getIndex(0,2);
      sigma_2 += S_c(i,i);
      data_t d_s = 0.822179*sqrt(sigma_2)/m(0,0);
      s -= d_s;
    }

    if (Config::USE_WCS)
      s /= d.scale_factor; // correct for WCS rescaling 

    // if we have the PSF image, we can use any scale
    if (mepsf.size() == K) {
      if (s >= meTroubleScale[k])
	s = (meSaveScale[k] + meTroubleScale[k])/2;
    }
    else { // if not we have to use the ones given to us
      s = psfs.getScaleClosestTo(s);
      if (s >= meTroubleScale[k]) {
	DEIMOS::PSFMultiScale::const_iterator iter = psfs.find(meTroubleScale[k]);
	s = iter->first;
	if (iter != psfs.begin()) {
	  iter--;
	  s = iter->first;
	}
      }
    }
    return s;
  }
}
