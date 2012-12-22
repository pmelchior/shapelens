#ifdef HAS_WCSLIB
#include "../include/WCSTransformation.h"
#include "../include/FITS.h"
#include "../include/MathHelper.h"
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>
#include <wcslib/wcs.h>
#include <iostream>
#include <iomanip>

namespace shapelens {

  WCSTransformation::WCSTransformation(fitsfile* fptr, bool intermediate_) : 
    intermediate(intermediate_)  {
    int status = 0, nkeyrec, nreject, nwcs;
    char *header;
    
    // Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
    if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
      throw std::runtime_error("WCSTransformation: Cannot read header of file " + FITS::getFileName(fptr));
    }

    // Interpret the WCS keywords
    struct wcsprm* wcss; // use pointer to read multiple wcs structs if necessary
    status = wcspih(header, nkeyrec, WCSHDR_all, 0, &nreject, &nwcs, &wcss);
    free(header);
    if (status)
      throw std::runtime_error("WCSTransformation: Cannot read WCS header keywords (" + std::string(wcshdr_errmsg[status]) + ")");

    // check there is one (and one only) WCS with 2 coordinate axes
    if (wcss == NULL) {
      throw std::runtime_error("WCSTransformation: No world coordinate systems found in " + FITS::getFileName(fptr));
    }
    else if (nwcs > 1) {
      wcsvfree(&nwcs, &wcss);
      throw std::runtime_error("WCSTransformation: More than one world coordinate systems found in " + FITS::getFileName(fptr));
    }
    else if (wcss->naxis != 2) {
      wcsvfree(&nwcs, &wcss);
      throw std::runtime_error("WCSTransformation: WCS does not have 2 axes in " + FITS::getFileName(fptr));
    }

    // initialize this wcs structure and copy it from first (and only)
    // entry of wcss
    wcs.flag = -1;           // wcslib implementation detail
    wcsini(1, 2, &wcs);      // 1: allocate memory for 2: axes
    wcscopy(0, wcss, &wcs);
    status = wcsset(&wcs);   // set remaining wcs fields from read information
    wcsvfree(&nwcs, &wcss);  // free the read-in structure

    if (status)
      throw std::runtime_error("WCSTransformation: wcsset error (" + std::string(wcs_errmsg[status]) + ")");

    // initialize coordinate containers
    world  = (double*) realloc(NULL,  2 * sizeof(double));
    imgcrd = (double*) realloc(NULL, 2 * sizeof(double));
    pixcrd = (double*) realloc(NULL, 2 * sizeof(double));
    stat   = (int*) realloc(NULL,   2 * sizeof(int));

    // since wcslib does not deal with SIP distortions, we have to read
    // the distortion coefficients if present
    std::string ctype1;
    FITS::readKeyword(fptr, "CTYPE1", ctype1);
    if (ctype1.find("TAN-SIP") != std::string::npos) {
      has_sip = true;
      std::ostringstream key;
      int A_order, B_order;
      FITS::readKeyword(fptr, "A_ORDER", A_order);
      FITS::readKeyword(fptr, "B_ORDER", B_order);
      A.resize(A_order+1, A_order+1);
      B.resize(B_order+1, B_order+1);
      A.setZero();
      B.setZero();
      for (int i=0; i <= A_order; i++) {
	for (int j=0; j <= A_order; j++) {
	  key.str("");
	  key << "A_" << i << "_" << j;
	  try { // not all coefficients need to be present
	    FITS::readKeyword(fptr, key.str(), A(i,j));
	  } catch (std::invalid_argument) {}
	}
      }
      for (int i=0; i <= B_order; i++) {
	for (int j=0; j <= B_order; j++) {
	  key.str("");
	  key << "B_" << i << "_" << j;
	  try {
	    FITS::readKeyword(fptr, key.str(), B(i,j));
	  } catch (std::invalid_argument) {}
	}
      }
      FITS::readKeyword(fptr, "CRPIX1", crpix1);
      FITS::readKeyword(fptr, "CRPIX2", crpix2);
    }
    else
      has_sip = false;
  }
  
  // explicit definition since we have to allocate containers
  // and perform a deep copy of wcs
  WCSTransformation::WCSTransformation(const WCSTransformation& W) :
    A(W.A), B(W.B), intermediate (W.intermediate), has_sip(W.has_sip), crpix1(W.crpix1), crpix2(W.crpix2) {
    world  = (double*) realloc(NULL,  2 * sizeof(double));
    imgcrd = (double*) realloc(NULL, 2 * sizeof(double));
    pixcrd = (double*) realloc(NULL, 2 * sizeof(double));
    stat   = (int*) realloc(NULL,   2 * sizeof(int));
    wcs.flag = -1;
    wcsini(1, 2, &wcs);
    wcscopy(0, &(W.wcs), &wcs);
    wcsset(&wcs);
  }
  
  // explicit definition to deallocate all structures
  WCSTransformation::~WCSTransformation() {
    int nwcs = 1;
    wcsfree(&wcs);
    free(world);
    free(imgcrd);
    free(pixcrd);
    free(stat);
  }

  data_t WCSTransformation::sip_polynomial(const tmv::Matrix<data_t>& M, data_t& u, data_t& v) const {
    data_t f = 0;
    for (int i=0; i < M.nrows(); i++)
      for (int j=0; j < M.ncols(); j++)
	f += M(i,j) * pow_int(u,i) * pow_int(v,j);
    return f;
  }

  void WCSTransformation::f(Point<data_t>& P) const {
    // apply sip transform if necessary
    if (has_sip) {
      // WCS pixels start at (1/1), so add 1
      double u = P(0) + 0 - crpix1, v = P(1) + 0 - crpix2;
      double f = sip_polynomial(A,u,v), g = sip_polynomial(B, u, v);
      *pixcrd = u + f + crpix1;
      *(pixcrd+1) = v + g + crpix2;
    } 
    else {
      *pixcrd = P(0);    
      *(pixcrd+1) = P(1);
    }
    
    // use intermediate world coordinates
    // rather then celestial
    if (intermediate) {
      linp2x(const_cast<linprm*>(&(wcs.lin)), 1, 2, pixcrd, imgcrd);
      P(0) = *imgcrd;
      P(1) = *(imgcrd+1);
    } else {
      double phi, theta;
      wcsp2s(const_cast<wcsprm*>(&wcs), 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
      P(0) = *world;
      P(1) = *(world+1);
    }
    stack_transform(P);
  }

  void WCSTransformation::f_1(Point<data_t>& P) const {
    // inverse: this trafo comes last
    stack_inverse_transform(P);

    // use intermediate world coordinates (as input)
    // rather than celestial
    if (intermediate) {
      *imgcrd = P(0);    
      *(imgcrd+1) = P(1);
      linx2p(const_cast<linprm*>(&(wcs.lin)), 1, 2, imgcrd, pixcrd);
    }
    else {
      *world = P(0);    
      *(world+1) = P(1);
      double phi, theta;
      wcss2p (const_cast<wcsprm*>(&wcs), 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
    }
    P(0) = (*pixcrd);
    P(1) = *(pixcrd+1);

    // need to invert non-linear SIP distortion if necessary
    if (has_sip) {
      data_t u, u_, u__, v, v_, v__;
      u = u_ = P(0) - crpix1;
      v = v_ = P(1) - crpix2;
      do { // simple solver, assumes distortion to be small
	u__ = u;
	v__ = v;
	u = u_ - sip_polynomial(A, u, v);
	v = v_ - sip_polynomial(B, u, v);
      } while (fabs(u - u__) > 1e-2 || fabs(v - v__) >1e-2);
      P(0) = u + crpix1;
      P(1) = v + crpix2;
    }

    // since WCS pixels start at (1/1), subtract of 1 to conform with shapelens
    P(0) -= 0;
    P(1) -= 0;
  }

  boost::shared_ptr<CoordinateTransformation> WCSTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new WCSTransformation(*this));
  }
} // end namespace

#endif // HAS_WCSLIB
