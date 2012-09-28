#ifndef SHAPELENS_WCSTRANSFORMATION_H
#define SHAPELENS_WCSTRANSFORMATION_H

#include "CoordinateTransformation.h"
#include <fitsio.h>
#ifdef HAS_WCSLIB
  #include <wcslib/wcs.h>
#endif

namespace shapelens {
  /// Class for WCS coordinate transformations.
  /// This class depends on the presence of libwcs (and its headers).
  /// To make use of it, define \p -DHAS_WCSLIB at compile time.
  class WCSTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    WCSTransformation(fitsfile* fptr, bool intermediate = true);
    /// Copy constructor.
    WCSTransformation(const WCSTransformation& wcs);
    ///Destructor.
    ~WCSTransformation(); 
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    bool intermediate;
    virtual void f(Point<data_t>& P) const;
    virtual void f_1(Point<data_t>& P) const;
#ifdef HAS_WCSLIB  
    struct wcsprm wcs;
    double *imgcrd, *pixcrd, *world;
    int *stat;
    tmv::Matrix<data_t> A,B;
    data_t crpix1, crpix2;
    bool has_sip;
    data_t sip_polynomial(const tmv::Matrix<data_t>& M, data_t& u, data_t& v) const;
#endif
  };
} // end namespace shapelens

#endif // SHAPELENS_WCSTRANSFORMATION_H

