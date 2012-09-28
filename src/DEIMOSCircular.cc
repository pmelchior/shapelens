#include "../include/DEIMOSCircular.h"
#include "../include/LensHelper.h"
#include "../include/FITS.h"

namespace shapelens {

  bool DEIMOSCircular::FIX_CENTROID = false;

  DEIMOSCircular::DEIMOSCircular() : DEIMOS(), scale(0), scale_factor(1) {}

  DEIMOSCircular::DEIMOSCircular(const Object& obj, int N_, int C_, data_t matching_scale_) : 
    DEIMOS(N_), C(C_), matching_scale(matching_scale_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2, 0) {
     
    // measure moments within optimized weighting function
    Moments mo_w(N+C);
    scale_factor = obj.grid.getScaleFactor();
    match(obj, mo_w);

    // compute noise statistics
    setNoiseImage(obj);
    computeCovariances(mo_w);
    SN[matching_scale] = computeSN(mo_w);
 }

  DEIMOSCircular::DEIMOSCircular (const Object& obj, const DEIMOSCircular::PSFMultiScale& psfs, int N_, int C_) :
    DEIMOS(N_), C(C_), matching_scale(psfs.getMaximumScale()),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2, 0) {
    
    // measure moments within optimized weighting function
    Moments mo_w(N+C);
    scale_factor = obj.grid.getScaleFactor();
    setNoiseImage(obj);
    match(obj, mo_w);

    // initial matching failed, try with smaller scale
    while (flags.any() && matching_scale > psfs.getMinimumScale()) {
      matching_scale = psfs.getScaleSmallerThan(matching_scale);
      history << "# Matching failed (" << flags << "), restarting with s = " << matching_scale << std::endl; 
      match(obj, mo_w);
    }

    // matching successfull, try with lower scale and see whether S/N improves
    if (flags.none()) {
      while (matching_scale > psfs.getMinimumScale()) {
	data_t SN_current = SN[matching_scale];
	data_t matching_scale_current = matching_scale;
	Point<data_t> centroid_current = centroid;
	Moments mo_current = mo;
	matching_scale = psfs.getScaleSmallerThan(matching_scale);
	history << "# Trying smaller scale s = " << matching_scale << " to find optimal S/N" << std::endl;
	match(obj, mo_w);
	// if S/N goes down, reset to best S/N case
	if (SN[matching_scale] < SN_current) {
	  matching_scale = matching_scale_current;
	  history << "# Reverting to scale s = " << matching_scale << std::endl;
	  mo = mo_current;
	  scale = matching_scale * scale_factor;
	  centroid = centroid_current;
	  break;
	}
      }
      history << "# Deweighted moments:\t" << mo << std::endl;

      // deconvolve and check moments
      Moments tmp = mo;
      deconvolve(psfs.getAtScale(matching_scale));
      history << "# Deconvolved moments:\t" << mo << std::endl;
      flags[1] = flagMoments(mo);
      while (flags.any() && matching_scale > psfs.getMinimumScale()) {
	matching_scale = psfs.getScaleSmallerThan(matching_scale);
	history << "# Deconvolution failed (" << flags << "), repeat matching with scale s = " << matching_scale << std::endl;
	match(obj, mo_w);
	history << "# Deweighted moments:\t" << mo << std::endl;
	tmp = mo;
	deconvolve(psfs.getAtScale(matching_scale));
	history << "# Deconvolved moments:\t" << mo << std::endl;
	flags[1] = flagMoments(mo);
      }
      if (flags.any()) {
	history << "# Deconvolution failed (" << flags << "), minimum PSF scale reached. GAME OVER." << std::endl;
	mo = tmp;
      }
      else
	computeCovariances(mo_w);
    }
    else
      history << "# Deweighting failed (" << flags << "), minimum PSF scale reached. GAME OVER." << std::endl;
    
  }

  DEIMOSCircular::DEIMOSCircular(std::string filename) {
    fitsfile* fptr = FITS::openFile(filename);
    tmv::Matrix<data_t> M;
    FITS::readImage(fptr,M);
    N = int(M.nrows()) - 1;
    mo.setOrder(N);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	mo(n1,n2) = M(n1,n2);

    // necessary parameters
    FITS::readKeyword(fptr,"SCALE_M",matching_scale);
    FITS::readKeyword(fptr,"C",C);

    try {
      FITS::readKeyword(fptr,"SCALEFAC",scale_factor);
    } catch (std::invalid_argument) {
      scale_factor = 1;
    }
    try {
      FITS::readKeyword(fptr,"SCALE",scale);
    } catch (std::invalid_argument) {
      scale = matching_scale * scale_factor;
    }
    try {
      FITS::readKeyword(fptr,"SN",SN[matching_scale]);
    } catch (std::invalid_argument) {}
    try {
      FITS::readKeyword(fptr,"R2",R2);
    } catch (std::invalid_argument) {
      R2 = 0;
    }
    try {
      unsigned long flagnum;
      FITS::readKeyword(fptr,"FLAGS",flagnum);
      flags = std::bitset<2>(flagnum);
    } catch (std::invalid_argument) {}

    // read covariance
    try {
      FITS::moveToExtension(fptr, 2);
      FITS::readImage(fptr,S);
    } catch (std::runtime_error) {}
	
    FITS::closeFile(fptr);
  }

  void DEIMOSCircular::save(std::string filename) const {
    fitsfile* fptr = FITS::createFile(filename);
    tmv::Matrix<data_t> M(N+1,N+1);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	M(n1,n2) = mo(n1,n2);
    FITS::writeImage(fptr,M);
    FITS::updateKeyword(fptr,"C",C,"deweighting correction order");
    FITS::updateKeyword(fptr,"SCALE_M",matching_scale,"matching_scale [pixel]");
    FITS::updateKeyword(fptr,"SCALE",scale,"actual weighting function width");
    FITS::updateKeyword(fptr,"SCALEFAC",scale_factor,"avg. WCS units/pixel");
    FITS::updateKeyword(fptr,"CENTROID",std::complex<data_t>(centroid(0), centroid(1)),"weight funciton centroid");
    std::map<data_t, data_t>::const_iterator iter = SN.find(matching_scale);
    if (iter != SN.end())
      FITS::updateKeyword(fptr,"SN",iter->second,"S/N");
    FITS::updateKeyword(fptr,"R2", R2, "resolution of the galaxy");
    FITS::updateKeyword(fptr,"FLAGS",flags.to_ulong(),"matching and deconvolution flags");
    FITS::writeImage(fptr,S,"VARIANCE");
    FITS::closeFile(fptr);
  }

  // determine optimal weighting parameters, centroid and size
  void DEIMOSCircular::match(const Object& obj, Moments& mo_w) {
    // obj.centroid is in pixel coordinates by convention
    // while shape measurement (including scale, centroid, ellipticity)
    // needs to be done in WCS coordinates.
    centroid = obj.centroid;
    if (Config::USE_WCS)
      obj.grid.getWCS().transform(centroid);
    scale = matching_scale * scale_factor;
    computeDeweightingMatrix(mo_w);

    flags.reset();
    int iter = 0, maxiter = 8, iter_initial = 0;
    data_t SN_initial = 0;
    Point<data_t> centroid_shift;
    history << "# Matching weight function: s = " << matching_scale << std::endl;
    history << "# iter\tscale\tcentroid\t\tS/N" << std::endl;
    history << "# " + std::string(70, '-') << std::endl;
    while (iter < maxiter) {
      // measure moments under weight function
      GaussianWeightFunction w(scale, centroid);
      mo_w = Moments (obj, w, N + C, &centroid);  
      deweight(mo_w);
      flags[1] = flagMoments(mo);

      if (FIX_CENTROID)
	break;

      SN[matching_scale] = computeSN(mo_w);
      history << "# " << iter+1 << "\t" << scale/scale_factor << "\t" << centroid;
      if (noise.size() > 0) { // only then SN_ is meaningful
	history << "\t";
	history << SN[matching_scale];
      }
      history << std::endl;

      // centroid correction
      centroid_shift(0) = mo(1,0)/mo(0,0);
      centroid_shift(1) = mo(0,1)/mo(0,0);
      data_t shift = sqrt(centroid_shift(0)*centroid_shift(0) + centroid_shift(1)*centroid_shift(1)) / scale_factor;
      if (shift > 5) {
	flags[0] = 1;
	break;
      }

      // shift centroid
      centroid += centroid_shift;

      // centroiding has converged: stop it
      if (shift < 1e-2)
	break;
      iter++;
    }
  }

  void DEIMOSCircular::computeDeweightingMatrix(const Moments& mo_w) {
    data_t s2 = scale*scale;
    data_t s4 = s2*s2;

    unsigned int i,j;
    for (int n=0; n <= N; n++) {
      for (int m=0; m <= n; m++) {
	i = mo.getIndex(m,n-m);
	j = mo_w.getIndex(m,n-m);
	D(i,j) = 1;
	if (C >= 2) {
	  j = mo_w.getIndex(m+2,n-m);
	  D(i,j) = 1./(2*s2);
	  // since the moments are of same order n
	  // and ordered wrt to (first/last) index
	  // moment (m+1,n-m+1) is just directly following in mo
	  j+=2;
	  D(i,j) = 1./(2*s2);
	}
	if (C >= 4) {
	  j = mo_w.getIndex(m+4,n-m);
	  D(i,j) = 1./(8*s4);
	  j+=2;
	  D(i,j) = 1./(4*s4);
	  j+=2;
	  D(i,j) = 1./(8*s4);
	}
	if (C >= 6) {
	  data_t s6 = s2*s2*s2;
	  j = mo_w.getIndex(m+6,n-m);
	  D(i,j) = 1./(48*s6);
	  j+=2;
	  D(i,j) = 1./(16*s6);
	  j+=2;
	  D(i,j) = 1./(16*s6);
	  j+=2;
	  D(i,j) = 1./(48*s6);
	}
	if (C >= 8) {
	  data_t s8 = s4*s4;
	  j = mo_w.getIndex(m+8,n-m);
	  D(i,j) = 1./(384*s8);
	  j+=2;
	  D(i,j) = 1./(96*s8);
	  j+=2;
	  D(i,j) = 1./(64*s8);
	  j+=2;
	  D(i,j) = 1./(96*s8);
	  j+=2;
	  D(i,j) = 1./(384*s8);
	}
	if (C >= 10) {
	  data_t s10 = s4*s4*s2;
	  j = mo_w.getIndex(m+10,n-m);
	  D(i,j) = 1./(3840*s10);
	  j+=2;
	  D(i,j) = 1./(768*s10);
	  j+=2;
	  D(i,j) = 1./(384*s10);
	  j+=2;
	  D(i,j) = 1./(384*s10);
	  j+=2;
	  D(i,j) = 1./(768*s10);
	  j+=2;
	  D(i,j) = 1./(3840*s10);
	}
      }
    }
  }
  void DEIMOSCircular::deweight(const Moments& mo_w) {
    mo = D*mo_w;
  }

  data_t DEIMOSCircular::computeSN(const Moments& mo_w) {
    // measure S/N (of flux)
    data_t SN_ = 0;
    // only do it if we have a noise image
    if (noise.size() > 0) {
      GaussianWeightFunction w2(scale/M_SQRT2, centroid);
      Moments mo_w2(noise, w2, 0, &centroid);
      SN_ = mo_w(0,0) / sqrt(mo_w2(0,0));
    }
    return SN_;
  }

  void DEIMOSCircular::computeCovariances(const Moments& mo_w) {
    // only do something if we have a noise image
    if (noise.size() > 0) {
      // covariances of a weighted moment mo(i,j) is given by
      // <mo(i,j) mo_(k,l)> = sigma_n^2 \int d^2x W^2(x) x_1^{i+k} x_2^{j+l},
      // i.e. moments up to (i*2,j*2) measured with square of weighting function
      // (square of Gaussian: sigma -> sigma/sqrt(2));
      GaussianWeightFunction w2(scale/M_SQRT2, centroid);
      Moments mo_noise(noise, w2, 2*(N+C), &centroid);
      tmv::Matrix<data_t> S_(mo_w.size(), mo_w.size());
      data_t det = 1;
      for (int n = 0; n < mo_w.size(); n++) {
	std::pair<int, int> p = mo_w.getPowers(n);
	for (int m = 0; m < mo_w.size(); m++) {
	  std::pair<int, int> p_ = mo_w.getPowers(m);
	  if (p.first + p.second + p_.first + p_.second <= 2*(N+C))
	    S_(n,m) = mo_noise(p.first + p_.first, p.second + p_.second);
	}
	det *= S_(n,n);
      }
      S = D*S_*D.transpose();
    }
  }
} // end namespace

