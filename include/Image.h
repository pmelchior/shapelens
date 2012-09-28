#ifndef SHAPELENS_IMAGE_H
#define SHAPELENS_IMAGE_H

#include <fitsio.h>
#include <string>
#include "ShapeLens.h"
#include "History.h"
#include "FITS.h"
#include "Grid.h"

namespace shapelens {

  /// Class for image data.
  /// Internally, the class represents image data as a vector, but the 
  /// interface enables two-dimensional access.\n
  /// The class can read in data directly from and save it to FITS files.
  template <class T> 
    class Image : public tmv::Vector<T> {
  public:
    /// Default constructor.
  Image() : tmv::Vector<T>() {
      history.clear();
    }
    /// Constructor with given image size \p NxM.
  Image(unsigned int N, unsigned int M) : tmv::Vector<T>(N*M) {
    }
    /// Argumented constructor for reading from a FITS file.
    /// The extension can be given in the standard cfitsio way as
    /// <tt>filename[extension]</tt>.
  Image(const std::string& filename) : tmv::Vector<T>() {
      read(filename);
    }
    /// Copy constructor from base class.
  Image(const tmv::Vector<T>& v) :  tmv::Vector<T>(v) {
    }
    
    /// Copy operator.
    Image<T>* operator=(const Image<T>& im) {
      tmv::Vector<T>::operator=(im);
      return this;
    } 
    /// Copy operator from base-class.
    Image<T>* operator=(const tmv::Vector<T>& v) {
      tmv::Vector<T>::operator=(v);
      return this;
    }

    /// Access operator using pixel index \p i.
    T& operator()(unsigned long i) {
      return tmv::Vector<T>::operator()(i);
    }
    /// const access operator using pixel index \p i.
    T operator()(unsigned long i) const {
      return tmv::Vector<T>::operator()(i);
    }
    /// Access operator using pixel coordinates.
    /// \b NOTE: Image is stored internally in row-major scheme, 
    /// but \p P emulates column-major as its first entry specifies the 
    /// column in the image. Thus, loops should iterate over columns first,
    /// then over rows.
    T& operator()(const Point<int>& P) {
      return tmv::Vector<T>::operator()(grid.getPixel(P));
    }
    /// const Access operator using pixel coordinates.
    /// \b NOTE: Image is stored internally in row-major scheme, 
    /// but \p P emulates column-major as its first entry specifies the 
    /// column in the image. Thus, loops should iterate over columns first,
    /// then over rows.
    T operator()(const Point<int>& P) const {
      return tmv::Vector<T>::operator()(grid.getPixel(P));
    }
    /// Matrix-style access operator.
    /// <tt>x,y</tt> are given as pixel offsets from the left-lower corner.\n\n
    /// \b NOTE: Image is stored internally in row-major scheme, 
    /// but <tt>x,y</tt> emulates column-major as \p x specifies the 
    /// column in the image. Thus, loops should iterate over \p x first,
    /// then over \p y.
    T& operator()(unsigned int x, unsigned int y) {
      Point<int>P(grid.getStartPosition(0)+int(x), grid.getStartPosition(1) + int(y));
      return tmv::Vector<T>::operator()(grid.getPixel(P));
    }
    /// const Matrix-style access operator.
    /// <tt>x,y</tt> are given as pixel offsets from the left-lower corner.\n\n
    /// \b NOTE: Image is stored internally in row-major scheme, 
    /// but <tt>x,y</tt> emulates column-major as \p x specifies the 
    /// column in the image. Thus, loops should iterate over \p x first,
    /// then over \p y.
    T operator()(unsigned int x, unsigned int y) const {
      Point<int>P(grid.getStartPosition(0)+int(x), grid.getStartPosition(1) + int(y));
      return tmv::Vector<T>::operator()(grid.getPixel(P));
    }
  
    /// Get value at image coordinates \p P.
    /// If \p P is outside the image, returns \p 0.
    T get(const Point<int>& P) const {
      long index = grid.getPixel(P);
      if (index == -1)
	return T(0);
      else
	return Image<T>::operator()(index);
    }
    /// Get value at arbitrary World coordinates \p P.
    /// If \p P is outside the image, returns \p 0.
    T get(const Point<data_t>& P) const {
      Point<int> IC = grid.getCoords(P);
      long index = grid.getPixel(IC);
      if (index == -1)
	return T(0);
      else
	return Image<T>::operator()(index);
    }

    /// Bi-linear interpolation at arbitrary World coordinate \p P.
    /// If \p P is outside the image, returns \p 0.
    T interpolate(const Point<data_t>& P) const {
      Point<int> IC = grid.getCoords(P);
      data_t x,y;
      const CoordinateTransformation& ct = grid.getWCS();
      Point<data_t> P_ = P;
      if (&ct != NULL)
	ct.inverse_transform(P_); // World -> pixel
      x = P_(0);
      y = P_(1);
      int x0 = IC(0), y0 = IC(1);
      int x1 = x0+1, y1 = y0+1; // neighborhood in image coords
      T f00,f01,f10,f11;
      long index = grid.getPixel(IC);
      if (index == -1) 
	return f00 = T(0);
      else 
	f00 = Image<T>::operator()(index);
      // right
      IC(0) = x1;
      index = grid.getPixel(IC);
      if (index == -1)
	f10 = T(0);
      else
	f10 = Image<T>::operator()(index);
      // top
      IC(0) = x0;
      IC(1) = y1;
      index = grid.getPixel(IC);
      if (index == -1)
	f01 = T(0);
      else
	f01 = Image<T>::operator()(index);
      // top-right
      IC(0) = x1;
      index = grid.getPixel(IC);
      if (index == -1)
	f11 = T(0);
      else
	f11 = Image<T>::operator()(index);

      return f00*T(x1-x)*T(y1-y) + f01*T(x1-x)*T(y-y0) + f10*T(x-x0)*T(y1-y) + f11*T(x-x0)*T(y-y0);
    }

    /// Slice a sub-image, specified by edge-points \p P1 and \p P2.
    template <class R>
      void slice(Image<R>& sub, const Point<int>& P1, const Point<int>& P2) const {
      int xmin = P1(0), xmax = P2(0), ymin = P1(1), ymax = P2(1);
      int axis0 = xmax-xmin;
      if (sub.grid.getSize(0) != (xmax-xmin) || sub.grid.getSize(1) != (ymax-ymin))  
	sub.resize((xmax-xmin)*(ymax-ymin));

      // loop over all object pixels
      R* ptr = sub.ptr(); 
      int step = sub.step();
      for (int i =0; i < sub.size(); i++) {
	// old coordinates derived from new pixel index i
	int x = i%axis0 + xmin;
	int y = i/axis0 + ymin;
	if (x>=0 && x < Image<T>::grid.getSize(0) && y >= 0 && y < Image<T>::grid.getSize(1))
	  *(ptr +i*step) = Image<T>::operator()(x,y);
      }
      sub.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
      sub.basefilename = basefilename;
      sub.history.setSilent();
      sub.history << "# Slice from " << basefilename << " in the area (";
      sub.history << xmin << "/" << ymin << ") -> (";
      sub.history << xmax << "/" << ymax << ")" << std::endl;
    }

  
    /// Get axis size of the whole image in given direction.
    unsigned int getSize(bool direction) const {
      return grid.getSize(direction);
    }

    /// Get the filename of the FITS file.
    std::string getFilename() const {
      return basefilename;
    }

    /// Save as FITS file.
    void save(std::string filename) const {
      fitsfile* fptr = FITS::createFile(filename);
      FITS::writeImage(fptr,*this);
      // if history is not empty, append history to FITS header
      if (!history.isEmpty())
	FITS::appendHistory(fptr,history.str());
      FITS::closeFile(fptr);
    }

    /// The Grid this Image is defined on.
    /// The grid is defined to range from \p 0 to getAxisSize(i)-1.
    Grid grid;
    /// The image history.
    History history;
    /// The file this data is derived from.
    /// The string should not exceed 80 characters, otherwise it becomes 
    /// eventually truncated during save().
    std::string basefilename;

  private:
    void read(const std::string& filename) {
      basefilename = filename;
      fitsfile *fptr = FITS::openFile(basefilename);
      FITS::readImage(fptr,*this);
      history << "# Reading FITS image " + basefilename << std::endl;
      FITS::closeFile(fptr);
    }
  };
} // end namespace
#endif
