#include "../include/PSFEx.h"
#include "../include/MathHelper.h"

namespace shapelens {
  PSFEx::PSFEx(std::string filename) :
    polzero(2), polscale(2), psfaxis(3) {

    fitsfile* fptr = FITS::openFile(filename);
    FITS::moveToExtension(fptr, 2);
    int polnaxis;
    FITS::readKeyword(fptr, "POLNAXIS", polnaxis);
    if (polnaxis != 2)
      throw std::invalid_argument("PSFEx: POLNAXIS != 2");
    FITS::readKeyword(fptr, "POLZERO1", polzero[0]);
    FITS::readKeyword(fptr, "POLSCAL1", polscale[0]);
    FITS::readKeyword(fptr, "POLZERO2", polzero[1]);
    FITS::readKeyword(fptr, "POLSCAL2", polscale[1]);
    std::vector<std::string> polname(2);
    FITS::readKeyword(fptr, "POLNAME1", polname[0]);
    FITS::readKeyword(fptr, "POLNAME2", polname[1]);
    int psfnaxis;
    FITS::readKeyword(fptr, "PSFNAXIS", psfnaxis);
    if (psfnaxis != 3)
      throw std::invalid_argument("PSFEx: PSFNAXIS != 3");
    FITS::readKeyword(fptr, "PSFAXIS1", psfaxis[0]);
    FITS::readKeyword(fptr, "PSFAXIS2", psfaxis[1]);
    FITS::readKeyword(fptr, "PSFAXIS3", psfaxis[2]);
    FITS::readKeyword(fptr, "POLDEG1", poldeg);
    if (psfaxis[2] != ((poldeg+1)*(poldeg+2))/2)
      throw std::invalid_argument("PSFEx: POLDEG and PSFAXIS3 disagree");
    FITS::readKeyword(fptr, "PSF_SAMP", psf_samp);

    // read basis functin shapes into images
    int status = 0, anynull;
    basis = std::vector<Image<data_t> >(psfaxis[2], Image<data_t>(psfaxis[0], psfaxis[1]));
    for (int k = 0; k < psfaxis[2]; k++)
      fits_read_col(fptr, FITS::getDataType(data_t(0)), 1, 1, 1 + k*psfaxis[0]*psfaxis[1], psfaxis[0]*psfaxis[1], NULL, basis[k].ptr(), &anynull, &status);

    FITS::closeFile(fptr);
  }

  data_t PSFEx::maxsize() {
    return ((psfaxis[0]-1)/2.-INTERPFAC)*psf_samp;
  }

  void PSFEx::fillObject(Object& psf) {
    // psf in sampling of PSFEx model
    Image<data_t> tmp = basis[0]; // set with zero-th mode
    data_t dx = (psf.centroid(0)-polzero[0])/polscale[0];
    data_t dy = (psf.centroid(0)-polzero[1])/polscale[1];

    // iterate through all orders from 1 to n (zero done above already)
    for(int n = 1; n <= poldeg; n++) {
      for(int ny =0 ; ny <= n; ny++) { // orders in y
	int nx = n - ny; // order in x so x^(nx)*y^(ny) is of order n=nx+ny
	int k = nx + ny*(poldeg+1) - (ny*(ny-1))/2;
	tmp += (pow_int(dx,nx) * pow_int(dy,ny)) * basis[k];
      }
    }

    // Lanczos3 interpolation to target pixel resolution
    Point<data_t> P;
    data_t sum = 0; // for PSF flux normalization
    for (int l=0; l < psf.size(); l++) {
      P = psf.grid(l);
      data_t relxp = (P(0) - psf.centroid(0))/psf_samp;
      data_t relyp = (P(1) - psf.centroid(1))/psf_samp;
      
      if (fabs(relxp) > (psfaxis[0]-1)/2.-INTERPFAC || fabs(relyp) > (psfaxis[1]-1)/2.-INTERPFAC)
	psf(l) = 0;
      else {
	// iterate over all neighboring pixels in PSFEx sampled grid
	for(int i=0; i<psfaxis[0]; i++) { // x coordinate on psf_sample scale
	  for(int j=0; j<psfaxis[1]; j++) { // y coordinate on psf_sample scale
	    data_t dx = fabs(i - 0.5*psfaxis[0] - relxp);
	    data_t dy = fabs(j - 0.5*psfaxis[1] - relyp);
	    if (dx<INTERPFAC && dy<INTERPFAC) {
	      data_t interpolant = sinc(dx)*sinc(dx/INTERPFAC)*sinc(dy)*sinc(dy/INTERPFAC);
	      psf(l) += tmp(i,j)*interpolant;
	    }
	  }
	}
	psf(l) /= pow2(psf_samp);
	sum += psf(l);
      }
    }

    // flux normalization
    psf *= 1./sum;
  }

  double PSFEx::sinc(double x) { // normalized sinc function, see PSFEx manual
    if (x<1e-5 && x>-1e-5)
      return 1.;
    return sin(x*M_PI)/(x*M_PI);
  }
}
