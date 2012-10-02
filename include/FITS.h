#ifndef SHAPELENS_FITS_H
#define SHAPELENS_FITS_H

#include <fitsio.h>
#include <iostream>
#include <string>
#include <map>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include "ShapeLens.h"
#include "Grid.h"
#include "WCSTransformation.h"

namespace shapelens {
  template <class T>
    class Image;

  /// Class for FITS-related functions.
  /// \b NOTE: The functions will throw a instance of \p std::exception 
  /// in case of a failure.
  class FITS {
    // helper function to get pointer to data section fits functions
    template<class T>
      inline static  T* data(T& val) { 
      return &val;
    }    
    template<class T>
      inline static  T* data(Point<T>& p) {
      return p.c_array();
    }
    inline static char* data(std::string s) {
      return const_cast<char*>(s.c_str());
    }
  
  public:
 
    /// Return FITS image format definition, based on the type of \p entry.
    /// - <tt>char, unsigned char -> BYTE_IMG</tt>
    /// - <tt>int -> SHORT_IMG</tt>
    /// - <tt>unsigned int -> USHORT_IMG</tt>
    /// - <tt>long -> LONG_IMG</tt>
    /// - <tt>unsigned long -> ULONG_IMG</tt>
    /// - <tt>float -> FLOAT_IMG</tt>
    /// - <tt>double -> DOUBLE_IMG</tt>
    template <typename T> inline static
      int getImageFormat(const T& entry) {
      // default type, uses template specialization for other types
      // see below
      return BYTE_IMG;
    }

    /// Return FITS image format definition, based on the type of \p entry.
    /// - <tt>bool, char, unsigned char -> TBYTE</tt>
    /// - <tt>int -> TINT</tt>
    /// - <tt>unsigned int -> TUINT</tt>
    /// - <tt>long -> TLONG</tt>
    /// - <tt>unsigned long -> TULONG</tt>
    /// - <tt>float -> TFLOAT</tt>
    /// - <tt>double -> TDOUBLE</tt>
    /// - <tt>std::complex<float> -> TCOMPLEX</tt>
    /// - <tt>std::complex<double> -> TDBLCOMPLEX</tt>
    /// - <tt>std::string -> TSTRING</tt>
    template <typename T> inline static
      int getDataType(const T& entry) {
      // default type, uses template specialization for other types
      // see below
      return TBYTE;
    }

    /// Open FITS file.
    /// If <tt>write == false</tt>, the file will be opened in read-only mode.
    static fitsfile* openFile(const std::string& filename, bool write=false);
    /// Open FITS file and go to first table extension.
    /// If <tt>write == false</tt>, the file will be opened in read-only mode.
    static fitsfile* openTable(const std::string& filename, bool write=false);
    /// Create new FITS file.
    /// If the file \p filename already exists, it will be overwritten.
    static fitsfile* createFile(const std::string& filename);
    /// Close FITS file pointer.
    static void closeFile(fitsfile* fptr);
    /// Move to extension \p i (starting from 1) in FITS file.
    static void moveToExtension(fitsfile* fptr, unsigned int i);
    /// Move to extension \p name in FITS file.
    static void moveToExtension(fitsfile* fptr, const std::string& name);
    /// Append \p history to FITS header histroy.
    static void appendHistory(fitsfile *fptr, const std::string& history);
    /// Read FITS keyword cards directly.
    static void readKeyCards(fitsfile *fptr, const std::string& key, std::string& value);
    /// Get name of FITS file from its pointer.
    static std::string getFileName(fitsfile *fptr);
    /// Get description error message from status code.
    static std::string getErrorMessage(int status);
    /// Read in keyword from FITS header.
    template <class T>
      static void readKeyword(fitsfile *fptr, const std::string& key, T& val) {
      int status = 0;
      char* comment = NULL;
      fits_read_key (fptr,getDataType(val), const_cast<char*>(key.c_str()),data(val),comment, &status);

      if (status != 0) {
	std::ostringstream note;
	note << "FITS: Cannot read FITS keyword " << key << " from " << getFileName(fptr) << ": " << getErrorMessage(status);
	throw std::invalid_argument(note.str());
      }
    }

    /// Set/update keyword in FITS file header.
    template <class T>
      static void updateKeyword(fitsfile *fptr, const std::string& keyword, const T& value_, std::string comment = "") {
      T& value = const_cast<T&>(value_);
      int status = 0;
      fits_write_key (fptr, getDataType(value), const_cast<char *>(keyword.c_str()), data(value), comment.c_str(), &status);

      if (status != 0) {
	std::ostringstream note;
	note << "FITS: Cannot update FITS keyword " << keyword << " = " << value_ << " in " << getFileName(fptr) << ": " << getErrorMessage(status);
	throw std::runtime_error(note.str());
      }
    }


    /// Write FITS image from an Image<T>.
    /// The datatype will be automatically adjusted, based on the
    /// result of getImageFormat() and getDataType().
    /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
    template <class T>
      static void writeImage(fitsfile *fptr, const Image<T>& image, std::string extname="") {
      int dim0 = image.grid.getSize(0);
      int dim1 = image.grid.getSize(1);
      long naxis = 2;      
      long naxes[2] = { dim0, dim1 };
      long npixels = dim0*dim1;

      // define image format and dataformat according to cfitsio definitions
      int imageformat = getImageFormat(image(0));
      int datatype = getDataType(image(0));
      // create HDU
      int status = 0;
      fits_create_img(fptr, imageformat, naxis, naxes, &status);
      // write pixel data
      long firstpix[2] = {1,1};
      fits_write_pix(fptr,datatype,firstpix,npixels,const_cast<T *>(image.cptr()), &status);
      // insert creator and extname keywords
      if (extname != "")
	updateKeyword(fptr, "EXTNAME", extname);
      updateKeyword(fptr, "CREATOR", std::string("shapelens"));

      if (status != 0)
	throw std::runtime_error("FITS: Cannot write FITS image " + extname + " to " + getFileName(fptr) + ": " + getErrorMessage(status));
    }
  
    /// Write FITS image from a tmv::Matrix<T>.
    /// The datatype will be automatically adjusted, based on the
    /// result of getImageFormat() and getDataType().
    /// If \p extname is non-empty, the keyword \p EXTNAME is set to \p extname.
    template <class T>
      static void writeImage(fitsfile *fptr, const tmv::Matrix<T>& M, std::string extname="") {
      int dim0 = M.ncols();
      int dim1 = M.nrows();
      long naxis = 2;      
      long naxes[2] = { dim0, dim1 };
      long npixels = dim0*dim1;

      // define image format and dataformat according to cfitsio definitions
      int imageformat = getImageFormat(M(0,0));
      int datatype = getDataType(M(0,0));
      // create HDU
      int status = 0;
      fits_create_img(fptr, imageformat, naxis, naxes, &status);
      // write pixel data
      long firstpix[2] = {1,1};
      fits_write_pix(fptr,datatype,firstpix,npixels,const_cast<T *>(M.cptr()), &status);
      // insert creator and extname keywords
      if (extname != "")
	updateKeyword(fptr, "EXTNAME", extname);
      updateKeyword(fptr, "CREATOR", std::string("shapelens"));
      if (status != 0)
	throw std::runtime_error("FITS: Cannot write FITS image " + extname + " to " + getFileName(fptr) + ": " + getErrorMessage(status));
    }

    /// Read FITS image into tmv::Matrix<T>.
    /// \p M is adjusted to hold the contents of the image; the image value are 
    /// automatically casted to the type \p T of \p M.
    template <class T>
      static void readImage(fitsfile *fptr, tmv::Matrix<T>& M) {
      int naxis, status = 0;
      fits_get_img_dim(fptr, &naxis, &status);
      if (naxis!=2)
	throw std::invalid_argument("FITS: naxis != 2. Pointer of " + getFileName(fptr) + " does not provide image");

      long naxes[2] = {1,1};
      fits_get_img_size(fptr, naxis, naxes, &status);
      M.resize(naxes[1],naxes[0]);
      long firstpix[2] = {1,1};
      T val;
      int imageformat = getImageFormat(val);
      int datatype = getDataType(val);
      fits_read_pix(fptr, datatype, firstpix, naxes[0]*naxes[1], NULL, M.ptr(), NULL, &status);
      if (status != 0)
	throw std::runtime_error("FITS: Cannot read FITS image from " + getFileName(fptr) + ": " + getErrorMessage(status));
    }

    /// Read FITS image into Image<T>.
    /// \p im is adjusted to hold the contents of the image; the image value are 
    /// automatically casted to the type \p T of \p im.\n
    /// \p Image<T>::grid is set to Grid(0,0,N,M), where \p N and \p M 
    /// are the row and column numbers of the FITS image.
    template <class T>
      static void readImage(fitsfile *fptr, Image<T>& im) {
      int naxis, status = 0;
      fits_get_img_dim(fptr, &naxis, &status);
      if (naxis!=2)
	throw std::invalid_argument("FITS: naxis != 2. Pointer of " + getFileName(fptr) + " does not provide image");

      long naxes[2] = {1,1};
      fits_get_img_size(fptr, naxis, naxes, &status);
      im.grid.setSize(0,0,naxes[0],naxes[1]);
      // set wcs from fits file if requested
      if (Config::USE_WCS) {
#ifdef HAS_WCSLIB
	im.grid.setWCS(WCSTransformation(fptr));
#else
	throw std::runtime_error("FITS: WCS usage requested, but HAS_WCSToolsLib not specified");
#endif
      }

      im.resize(im.grid.size());
      long firstpix[2] = {1,1};
      T val;
      int datatype = getDataType(val);
      fits_read_pix(fptr, datatype, firstpix, im.size(), NULL, im.ptr(), NULL, &status);
      if (status != 0)
	throw std::runtime_error("FITS: Cannot read FITS image from "+ getFileName(fptr) + ": " + getErrorMessage(status));
    }

    /// Get number of rows in FITS table.
    static long getTableRows(fitsfile* fptr);
    /// Get column number of a FITS table with given \p name.
    /// \p name can contain \p * wildcards and is treated case-insensitive.
    static int getTableColumnNumber(fitsfile* fptr, const std::string& name);
    /// Get data type of column \p colnr from a FITS table.
    static int getTableColumnType(fitsfile* fptr, int colnr);
    /// Read \p val from FITS table at \p row and \p colnr.
    /// If a NULL value was stored at this position, \p nullvalue will be 
    /// returned insted.\n
    /// \p row numbers start with 0 according to regular C-style iterations. 
    template <class T>
      static void readTableValue(fitsfile* fptr, long row, int colnr, T& val, T nullvalue = 0) {
      int status = 0, anynull;
      fits_read_col(fptr, getDataType(val), colnr, row+1, 1, 1, &nullvalue, &val, &anynull, &status);
      if (status != 0) {
	std::ostringstream note;
	note << "FITS: Cannot read value in (row/col) = (" << row << "/" << colnr << ") from FITS table in " << getFileName(fptr) << ": " << getErrorMessage(status);
	throw std::runtime_error(note.str());
      }
    }
  };

  // template specializations
  template<> inline
    int FITS::getImageFormat<int>(const int& entry) {
    return SHORT_IMG;
  }
  template<> inline
    int FITS::getImageFormat<unsigned int>(const unsigned int& entry) {
    return USHORT_IMG;
  }
  template<> inline
    int FITS::getImageFormat<long>(const long& entry) {
    return LONG_IMG;
  }
  template<> inline
    int FITS::getImageFormat<unsigned long>(const unsigned long& entry) {
    return SHORT_IMG;
  }
  template<> inline
    int FITS::getImageFormat<float>(const float& entry) {
    return FLOAT_IMG;
  }
  template<> inline
    int FITS::getImageFormat<double>(const double& entry) {
    return DOUBLE_IMG;
  }
  

  template<> inline
    int FITS::getDataType<int>(const int& entry) {
    return TINT;
  }
  template<> inline
    int FITS::getDataType<unsigned int>(const unsigned int& entry) {
    return TUINT;
  }
  template<> inline
    int FITS::getDataType<long>(const long& entry) {
    return TLONG;
  }
  template<> inline
    int FITS::getDataType<unsigned long>(const unsigned long& entry) {
    return TULONG;
  }
  template<> inline
    int FITS::getDataType<float>(const float& entry) {
    return TFLOAT;
  }
  template<> inline
    int FITS::getDataType<double>(const double& entry) {
    return TDOUBLE;
  }
  template<> inline
    int FITS::getDataType<std::complex<float> >(const std::complex<float>& entry) {
    return TCOMPLEX;
  }
  template<> inline
    int FITS::getDataType<std::complex<double> >(const std::complex<double>& entry) {
    return TDBLCOMPLEX;
  }
  template<> inline
    int FITS::getDataType<Point<float> >(const Point<float>& entry) {
    return TCOMPLEX;
  }
  template<> inline
    int FITS::getDataType<Point<double> >(const Point<double>& entry) {
    return TDBLCOMPLEX;
  }
  template<> inline
    int FITS::getDataType<std::string>(const std::string& entry) {
    return TSTRING;
  }

  template <> inline
    void FITS::readKeyword<std::string>(fitsfile *fptr, const std::string& key, std::string& val) {
    int status = 0;
    char value[FLEN_CARD];
    fits_read_key (fptr,FITS::getDataType(val), const_cast<char *>(key.c_str()),&value, NULL, &status);
    val = std::string(value);
  }

} // end namespace
#endif
