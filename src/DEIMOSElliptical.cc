#include "../include/DEIMOSElliptical.h"
#include "../include/LensHelper.h"
#include "../include/FITS.h"

namespace shapelens {

  bool DEIMOSElliptical::FIX_CENTROID = false;

  DEIMOSElliptical::WeightFunction::WeightFunction(data_t scale, const Point<data_t>& centroid_, const complex<data_t>& eps) :
    GaussianWeightFunction(scale,Point<data_t>(0,0)),
    centroid(centroid_),
    T(0,eps) { 
  }

  data_t DEIMOSElliptical::WeightFunction::operator()(const Point<data_t>& P_) const {
    Point<data_t> P = P_;
    P -= centroid;
    T.inverse_transform(P);
    return GaussianWeightFunction::operator()(P);
  }


  DEIMOSElliptical::DEIMOSElliptical(int N_, int C_) : 
    DEIMOS(N_), C(C_), eps(0,0), scale_factor(1),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2, 0) { 
  }

  DEIMOSElliptical::DEIMOSElliptical(const Object& obj, int N_, int C_, data_t matching_scale_) : 
    DEIMOS(N_), C(C_), eps(0,0), matching_scale(matching_scale_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2, 0) {

    // measure moments within optimized weighting function
    Moments mo_w(N+C);
    scale_factor = obj.grid.getScaleFactor();
    match(obj, mo_w);

    setNoiseImage(obj);
    computeCovariances(mo_w);
    SN[matching_scale] = computeSN(mo_w);
 }

  DEIMOSElliptical::DEIMOSElliptical (const Object& obj, const DEIMOSElliptical::PSFMultiScale& psfs, int N_, int C_) :
    DEIMOS(N_), eps(0,0), C(C_), matching_scale(psfs.getMaximumScale()),
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
      while (matching_scale > psfs.getMinimumScale() && SN[matching_scale] < 80) {
	data_t SN_current = SN[matching_scale];
	data_t matching_scale_current = matching_scale;
	complex<data_t> eps_current = eps;
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
	  eps = eps_current;
	  scale = getEpsScale();
	  centroid = centroid_current;
	  break;
	}
      }
      history << "# Deweighted moments:\t" << mo << std::endl;
      Moments tmp = mo;
      deconvolve(psfs.getAtScale(matching_scale));
      flags[3] = flagMoments(mo);
      while (flags.any() && matching_scale > psfs.getMinimumScale()) {
	matching_scale = psfs.getScaleSmallerThan(matching_scale);
	history << "# Deconvolution failed, repeat matching with scale s = " << matching_scale << std::endl;
	// FIXME: if we store centroid, mo, eps, scale,
	// we may avoid matching in cases where smaller scales
	// have already been tried to improve S/N
	match(obj, mo_w);
	history << "# Deweighted moments:\t" << mo << std::endl;
	tmp = mo;
	deconvolve(psfs.getAtScale(matching_scale));
	flags[3] = flagMoments(mo);
      }
      history << "# Deconvolved moments:\t" << mo << std::endl;
      if (flags.any()) {
	history << "# Deconvolution failed, minimum PSF scale reached. GAME OVER." << std::endl;
	mo = tmp;
      }
      else
	computeCovariances(mo_w);
    }
    else
      history << "# Deweighting failed (" << flags << "), minimum PSF scale reached. GAME OVER." << std::endl;
    
  }

  DEIMOSElliptical::DEIMOSElliptical(std::string filename) {
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
      FITS::readKeyword(fptr,"EPS",eps);
    } catch (std::invalid_argument) {
      eps = complex<data_t>(0,0);
    }
    try {
      FITS::readKeyword(fptr,"SCALEFAC",scale_factor);
    } catch (std::invalid_argument) {
      scale_factor = 1;
    }
    try {
      FITS::readKeyword(fptr,"SCALE",scale);
    } catch (std::invalid_argument) {
      scale = getEpsScale();
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
      flags = std::bitset<4>(flagnum);
    } catch (std::invalid_argument) {}

    // read covariance
    try {
      FITS::moveToExtension(fptr, 2);
      FITS::readImage(fptr,S);
    } catch (std::runtime_error) {}
	
    FITS::closeFile(fptr);
  }

  void DEIMOSElliptical::save(std::string filename) const {
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
    FITS::updateKeyword(fptr,"EPS",eps,"weighting function ellipticity");
    std::map<data_t, data_t>::const_iterator iter = SN.find(matching_scale);
    if (iter != SN.end())
      FITS::updateKeyword(fptr,"SN",iter->second,"S/N");
    FITS::updateKeyword(fptr,"R2", R2, "resolution of the galaxy");
    FITS::updateKeyword(fptr,"FLAGS",flags.to_ulong(),"matching and deconvolution flags");
    FITS::writeImage(fptr,S,"VARIANCE");
    FITS::closeFile(fptr);
  }

  data_t DEIMOSElliptical::getEpsScale() const {
    // this scale preserves the size of semi-minor axis of ellipsoid
    // as the smallest scale (= matching_scale)
    data_t abs_eps = abs(eps);
    return scale_factor*matching_scale/sqrt(1 + abs_eps*abs_eps - 2*abs_eps);
  }

  // determine optimal weighting parameters, centroid, and ellipticity
  void DEIMOSElliptical::match(const Object& obj, Moments& mo_w) {
 
    // (re-)set centroid, ellipticity and scale
    centroid = obj.centroid;
    eps = complex<data_t>(0,0);
    flags.reset();

    // obj.centroid is in pixel coordinates by convention
    // while shape measurement (including scale, centroid, ellipticity)
    // needs to be done in WCS coordinates.
    if (Config::USE_WCS)
      obj.grid.getWCS().transform(centroid);

    // use a smaller scale for centroiding
    // but maintain a sensible width (of 1 pixel) otherwise the
    // cenroid could become dominated by local noise fluctuations 
    scale = std::max(matching_scale*0.66, 1.0) * scale_factor;

    int iter = 0, maxiter = 18, maxiter_centroid = 4, iter_initial = 0;
    bool centroiding = true;
    data_t SN_initial = 0;
    Point<data_t> centroid_shift;
    history << "# Matching weight function: s = " << matching_scale << std::endl;
    history << "# iter\tscale\tcentroid\t\tepsilon\t\t\tS/N" << std::endl;
    history << "# " + std::string(70, '-') << std::endl;
    while (iter < maxiter) {
      // measure moments under weight function
      WeightFunction w(scale, centroid, eps);
      mo_w = Moments (obj, w, N + C, &centroid);
      data_t SN_ = computeSN(mo_w);

      history << "# " << iter+1 << "\t" << scale/scale_factor << "\t" << centroid << "\t" << eps;
      if (noise.size() > 0) { // only then SN_ is meaningful
	history << "\t";
	if (centroiding)
	  history << "\t\t";
	history << SN_;
      }
      history << std::endl;

      // deweight and estimate new centroid and ellipticity
      deweight(mo_w);
      bool trouble = flagMoments(mo);

      // first: centroiding
      if (centroiding && iter < maxiter_centroid - 1 && FIX_CENTROID == false)  {
	centroid_shift(0) = mo(1,0)/mo(0,0);
	centroid_shift(1) = mo(0,1)/mo(0,0);
	data_t shift = sqrt(centroid_shift(0)*centroid_shift(0) + centroid_shift(1)*centroid_shift(1)) / scale_factor;
	if (shift > 5) {
	  flags[0] = 1;
	  mo(0,0) = mo_w(0,0);
	  mo(0,2) = mo_w(0,2);
	  mo(1,1) = mo_w(1,1);
	  mo(2,0) = mo_w(2,0);
	  break;
	}
	if (iter == maxiter_centroid - 1 && trouble) {
	  flags[1] = trouble;
	  mo(0,0) = mo_w(0,0);
	  mo(0,2) = mo_w(0,2);
	  mo(1,1) = mo_w(1,1);
	  mo(2,0) = mo_w(2,0);
	  break;
	}

	// shift centroid
	centroid += centroid_shift;

	// centroiding has converged: stop it
	if (shift < 1e-2)
	  centroiding = false;
      }

      // then: match ellipticity
      else if (iter < maxiter - 1) {
	complex<data_t> eps_ = epsilon();

	// abort for non-sensical ellipticities
	if (trouble) {
	  flags[1] = trouble;
	  mo(0,0) = mo_w(0,0);
	  mo(0,2) = mo_w(0,2);
	  mo(1,1) = mo_w(1,1);
	  mo(2,0) = mo_w(2,0);
	  break;
	}
	// convergence test
	if (iter_initial == 0) {
	  iter_initial = iter;
	  centroiding = false;
	}
	if (iter == iter_initial + 1)
	  SN_initial = SN_;
	if (iter > iter_initial + 1) {
	  if (fabs(SN_/SN[matching_scale] - 1) < 1e-5 && 
	      abs(eps - eps_) < 1e-3) {
	    SN[matching_scale] = SN_;
	    break;
	  }
	}
	eps = eps_;
	scale = getEpsScale();
      }
      SN[matching_scale] = SN_;
      iter++;
    }
  }

  void DEIMOSElliptical::deweight(const Moments& mo_w) {
    data_t e1 = real(eps);
    data_t e2 = imag(eps);
    data_t c1 = (1-e1)*(1-e1) + e2*e2;
    data_t c2 = (1+e1)*(1+e1) + e2*e2;
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
	  D(i,j) = c1/(2*s2);
	  // since the moments are of same order n
	  // and ordered wrt to (first/last) index
	  // moment (m+1,n-m+1) is just directly following in mo
	  j++;
	  D(i,j) = - 2*e2/s2;
	  j++;
	  D(i,j) = c2/(2*s2);
	}
	if (C >= 4) {
	  j = mo_w.getIndex(m+4,n-m);
	  D(i,j) = c1*c1/(8*s4);
	  j++;
	  D(i,j) = -c1*e2/s4;
	  j++;
	  D(i,j) = (c1*c2/4 + 2*e2*e2)/s4;
	  j++;
	  D(i,j) = -c2*e2/s4;
	  j++;
	  D(i,j) = c2*c2/(8*s4);
	}
	if (C >= 6) {
	  data_t s6 = s2*s2*s2;
	  j = mo_w.getIndex(m+6,n-m);
	  D(i,j) = c1*c1*c1/(48*s6);
	  j++;
	  D(i,j) = -c1*c1*e2/(4*s6);
	  j++;
	  D(i,j) = (c1*c1*c2/16 + c1*e2*e2)/s6;
	  j++;
	  D(i,j) = -(c1*c2*e2/2 + 4*e2*e2*e2/3)/s6;
	  j++;
	  D(i,j) = (c1*c2*c2/16 + c2*e2*e2)/s6;
	  j++;
	  D(i,j) = -c2*c2*e2/(4*s6);
	  j++;
	  D(i,j) = c2*c2*c2/(48*s6);
	}
	if (C >= 8) {
	  data_t s8 = s4*s4;
	  j = mo_w.getIndex(m+8,n-m);
	  D(i,j) = c1*c1*c1*c1/(384*s8);
	  j++;
	  D(i,j) = -c1*c1*c1*e2/(24*s8);
	  j++;
	  D(i,j) = (c1*c1*c1*c2/96 + c1*c1*e2*e2/4)/s8;
	  j++;
	  D(i,j) = -(c1*c1*c2*e2/8 + 2*c1*e2*e2*e2/3)/s8;
	  j++;
	  D(i,j) = (c1*c1*c2*c2/64 + c1*c2*e2*e2/2 + 2*e2*e2*e2*e2/3)/s8;
	  j++;
	  D(i,j) = -(c1*c2*c2*e2/8 + 2*c2*e2*e2*e2/3)/s8;
	  j++;
	  D(i,j) = (c1*c2*c2*c2/96 + c2*c2*e2*e2/4)/s8;
	  j++;
	  D(i,j) = -c2*c2*c2*e2/(24*s8);
	  j++;
	  D(i,j) = c2*c2*c2*c2/(384*s8);
	}
	if (C >= 10) {
	  data_t s10 = s4*s4*s2;
	  j = mo_w.getIndex(m+10,n-m);
	  D(i,j) = c1*c1*c1*c1*c1/(3840*s10);
	  j++;
	  D(i,j) = -c1*c1*c1*c1*e2/(192*s10);
	  j++;
	  D(i,j) = (c1*c1*c1*c1*c2/768 + c1*c1*c1*e2*e2/24)/s10;
	  j++;
	  D(i,j) = -(c1*c1*c1*c2*e2/48 + c1*c1*e2*e2*e2/6)/s10;
	  j++;
	  D(i,j) = (c1*c1*c1*c2*c2/384 + c1*c1*c2*e2*e2/8 + c1*e2*e2*e2*e2/3)/s10;
	  j++;
	  D(i,j) = -(c1*c1*c2*c2*e2/32 + c1*c2*e2*e2*e2/3 + 4*e2*e2*e2*e2*e2/15)/s10;
	  j++;
	  D(i,j) = (c1*c1*c2*c2*c2/384 + c1*c2*c2*e2*e2/8 + c2*e2*e2*e2*e2/3)/s10;
	  j++;
	  D(i,j) = -(c1*c2*c2*c2*e2/48 + c2*c2*e2*e2*e2/6)/s10;
	  j++;
	  D(i,j) = (c1*c2*c2*c2*c2/768 + c2*c2*c2*e2*e2/24)/s10;
	  j++;
	  D(i,j) = -c2*c2*c2*c2*e2/(192*s10);
	  j++;
	  D(i,j) = c2*c2*c2*c2*c2/(3840*s10);
	}
      }
    }

    mo = D*mo_w;
  }

  data_t DEIMOSElliptical::computeSN(const Moments& mo_w) {
    // measure S/N (of flux)
    data_t SN_ = 0;
    // only do it if we have a noise image
    if (noise.size() > 0) {
      WeightFunction w2(scale/M_SQRT2, centroid, eps);
      Moments mo_w2(noise, w2, 0, &centroid);
      SN_ = mo_w(0,0) / sqrt(mo_w2(0,0));
    }
    return SN_;
  }

  void DEIMOSElliptical::computeCovariances(const Moments& mo_w) {
    // only do something if we have a noise image
    if (noise.size() > 0) {
      // covariances of a weighted moment mo(i,j) is given by
      // <mo(i,j) mo_(k,l)> = sigma_n^2 \int d^2x W^2(x) x_1^{i+k} x_2^{j+l},
      // i.e. moments up to (i*2,j*2) measured with square of weighting function
      // (square of Gaussian: sigma -> sigma/sqrt(2));
      WeightFunction w2(scale/M_SQRT2, centroid, eps);
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

