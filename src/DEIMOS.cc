#include "../include/DEIMOS.h"
#include "../include/LensHelper.h"
#include "../include/FITS.h"
#include "../include/MathHelper.h"

namespace shapelens {

  DEIMOS::DEIMOS (int N_) : N(N_), mo(N_), R2(0), S(((N_+1)*(N_+2))/2, ((N_+1)*(N_+2))/2) {
  }

  DEIMOS::DEIMOS(std::string filename) {
    fitsfile* fptr = FITS::openFile(filename);
    tmv::Matrix<data_t> M;
    FITS::readImage(fptr,M);
    N = int(M.nrows()) - 1;
    mo.setOrder(N);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	mo(n1,n2) = M(n2,n1);
    try {
      FITS::readKeyword(fptr, "R2", R2);
    } catch (std::runtime_error) {}
    // read covariance
    try {
      FITS::moveToExtension(fptr, 2);
      FITS::readImage(fptr,S);
    } catch (std::runtime_error) {}
	
    FITS::closeFile(fptr);
  }

  void DEIMOS::save(std::string filename) const {
    fitsfile* fptr = FITS::createFile(filename);
    tmv::Matrix<data_t> M(N+1,N+1);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	M(n2,n1) = mo(n1,n2); // transpose to have correctly oriented image
    FITS::writeImage(fptr,M);
    FITS::updateKeyword(fptr, "R2", R2, "resolution");
    FITS::writeImage(fptr,S,"VARIANCE");
    FITS::closeFile(fptr);
  }

  void DEIMOS::setNoiseImage(const Object& obj) {
    // only do something if noise properties are set
    if (obj.noise_rms > 0) {
      noise.resize(obj.size());
      noise.grid = obj.grid;
      noise.centroid = obj.centroid;
      noise.noise_rms = obj.noise_rms;
      for (unsigned int i=0; i < noise.size(); i++)
	noise(i) = obj.noise_rms*obj.noise_rms;
      if (obj.segmentation.size() == obj.size())
	noise.segmentation = obj.segmentation;
      noise.id = obj.id;
    }
  }

  Moments DEIMOS::getMomentErrors() const { 
    Moments mo_noise(N);
    if (S(0,0) != 0) {
      tmv::Matrix<data_t> S_1 = S.inverse();
      for (int i=0; i < mo_noise.size(); i++)
	mo_noise(i) = 1./S_1(i,i);
    }
    return mo_noise;
  }

  void DEIMOS::convolve(const Moments& p) {
    // create matrix representation of convolution eq. 9
    tmv::Matrix<data_t> P (mo.size(), mo.size());
    for (int n = 0; n <= N; n++) {
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
    mo = P*mo;
    R2 = shapelens::R2(mo, p);
  }

  void DEIMOS::deconvolve(const Moments& p) {
    Moments& g = mo;
    int Nmin = std::min(p.getOrder(),g.getOrder());
    R2 = shapelens::R2(mo, p);

    //use explicit relations for up to 2nd moments
    g(0,0) /= p(0,0);
    if (Nmin >= 1) {
      g(0,1) -= g(0,0)*p(0,1);
      g(0,1) /= p(0,0);
      g(1,0) -= g(0,0)*p(1,0);
      g(1,0) /= p(0,0);
      if (Nmin >= 2) {
	g(0,2) -= g(0,0)*p(0,2) + 2*g(0,1)*p(0,1);
	g(0,2) /= p(0,0);
	g(1,1) -= g(0,0)*p(1,1) + g(0,1)*p(1,0) + g(1,0)*p(0,1);
	g(1,1) /= p(0,0);
	g(2,0) -= g(0,0)*p(2,0) + 2*g(1,0)*p(1,0);
	g(2,0) /= p(0,0);
	if (Nmin >= 3) {
	  //use general formula (9)
	    for (int n=3; n <= Nmin; n++) {
	      for (int i=0; i <= n; i++) {
		int j = n-i;
		for (int k=0; k <= i-1; k++)
		  for (int l=0; l <= j-1; l++)
		    g(i,j) -= binomial(i,k)*binomial(j,l)*g(k,l)*p(i-k,j-l);
		for (int k=0; k <= i-1; k++)
		  g(i,j) -= binomial(i,k)*g(k,j)*p(i-k,0);
		for (int l=0; l <= j-1; l++)
		  g(i,j) -= binomial(j,l)*g(i,l)*p(0,j-l);
		g(i,j) /= p(0,0);
	      }
	    }
	}
      }
    }
  }
  
  // check whether moments make sense
  bool DEIMOS::flagMoments(const Moments& M) const {
    // FIXME: for higher order moments could also check M(4,0) ...
    if (M(0,0) < 0 || M(2,0) < 0 || M(0,2) < 0 || M(1,1)*M(1,1) > M(2,0)*M(0,2))
      return true;
    else
      return false;
  }

  data_t DEIMOS::size() const {
    return shapelens::scale(mo);
  }

  complex<data_t> DEIMOS::epsilon() const {
    return shapelens::epsilon(mo);
  }
  
  complex<data_t> DEIMOS::chi() const {
    return shapelens::chi(mo);
  }
 
  complex<data_t> DEIMOS::zeta() const {
    return shapelens::zeta(mo);
  }

  complex<data_t> DEIMOS::delta() const {
    return shapelens::delta(mo);
  }


  // PSF multiscale helper class
  void DEIMOS::PSFMultiScale::insert(data_t scale, const Moments& mo) {
    std::map<data_t, Moments>::insert(std::pair<data_t, Moments>(scale, mo));
  }
  const Moments& DEIMOS::PSFMultiScale::getAtScale(data_t scale) const {
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::find(scale);
    if (iter == std::map<data_t, Moments>::end())
      throw std::invalid_argument("DEIMOS::PSFMultiScale: scale not available!");
    return iter->second;
  }

  data_t DEIMOS::PSFMultiScale::getScaleSmallerThan(data_t scale) const {
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::find(scale);
    if (iter == std::map<data_t, Moments>::end())
      throw std::invalid_argument("DEIMOS::PSFMultiScale: scale not available!");
    if (iter == std::map<data_t, Moments>::begin())
      throw std::runtime_error("DEIMOS::PSFMultiScale: no smaller scale available!");
    iter--;
    return iter->first;

  }

  data_t DEIMOS::PSFMultiScale::getMinimumScale() const {
    return std::map<data_t, Moments>::begin()->first;
  }
  
  data_t DEIMOS::PSFMultiScale::getMaximumScale() const {
    return std::map<data_t, Moments>::rbegin()->first;
  }

  data_t DEIMOS::PSFMultiScale::getScaleClosestTo(data_t scale) const {
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::lower_bound(scale);
    // iter points to element >= scale;
    data_t diff = (iter->first)-scale;
    if (diff < 0) { // bigger than largest element
      return std::map<data_t, Moments>::rbegin()->first;
    }
    else {
      // check whethere there is a smaller scale
      if (iter != std::map<data_t, Moments>::begin()) {
	iter--;
	data_t diff_small = scale - (iter->first);
	if (diff < diff_small)
	  iter++;
      }
      return iter->first;
    }
  }

} // end namespace
