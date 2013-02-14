#ifndef SHAPELENS_SEXFRAME_H
#define SHAPELENS_SEXFRAME_H

#include <vector>
#include <list>
#include "ShapeLens.h"
#include "Grid.h"
#include "Image.h"
#include "Object.h"
#include "Catalog.h"

namespace shapelens {
  /// Wrapper class for SExtractor.
  /// Provides segmentation of an image into various Object instances
  /// from SExtractor's output files.
  /// It is assumed that for each image there's a SExtractor catalog file
  /// (in \p ASCII_HEAD or FITS format) and a segmentation map 
  /// (which SExtractor generates if \p CHECKIMAGE_TYPE \p SEGMENTATION 
  /// is specified in its configuration file).
  /// See Catalog and CatObject for the required catalog parameters.
  ///
  /// The Object entities will have these features:
  /// - The image cutout is quadratic such that the area within the object 
  ///   was detected by SExtractor (<tt>XMIN..XMAX,YMIN..YMAX</tt>) 
  ///   plus additional border area, specified by Config::ADD_BORDER, 
  ///   is included. Areas outside of the original image are set to \p 0.
  /// - The segmentation map corresponds to the same cutout area. Areas outside
  ///   of the original image are set to \p -1, bad pixels to \p -2.
  /// - If a weight map is provided, its cutout also corresponds to the image
  ///   cutout, areas outside of the original image are set to \p 0.
  ///
  /// Example usage:
  /// \code
  /// SExFrame f(fitsfile, catfile, segmap);
  /// const Catalog& cat = f.catalog;
  /// Catalog::const_iterator iter;
  /// for(iter = cat.begin(); iter != cat.end(); iter++) {
  ///   unsigned long id = iter->first;          // just for clarity
  ///   const CatObject& catobj = iter->second;  // just for clarity 
  ///   Object obj;
  ///   f->fillObject(obj,iter);
  ///   if (catobj.FLAGS == 0) {
  ///     // do something useful with obj ...
  ///   }
  /// }
  ///\endcode
  class SExFrame {
  public:
    /// Argumented constructor for including a weight map image.
    /// \p datafile, \p weightfile and \p segmapfile are the names of the data,
    /// weight map and segmenation map FITS files 
    /// (Extensions or other selections can be passed in the standard cfitsio 
    /// way: \p filename[extension] ).
    /// \p catfile is the name of the SExtractor catalog file.
    /// \p weightfile and \p segmapfile can be empty strings (<tt>""</tt>), then
    /// they will be ignored.
    SExFrame(std::string datafile, std::string catfile, std::string segmapfile = "", std::string weightfile = "");
    /// Destructor.
    ~SExFrame();
    /// Fill all necessary information into an Object.
    /// The object is selected by the iterator \p catiter, which must must
    /// be constructed from SExFrame::catalog.\n
    void fillObject(Object& obj, Catalog::const_iterator& catiter);
    /// Subtract a constant background level from each object in fillObject().
    void subtractBackground(data_t bg);
    /// Return number of objects found by SExtractor.
    /// If it returns 0, the catalog file is either empty or its format
    /// is wrongly specified.
    unsigned long getNumberOfObjects();
    /// Grid.
    Grid grid;
    /// The catalog of detected objects.
    Catalog catalog;

  private:
    fitsfile *fptr, *fptr_w, *fptr_s;
    History history;
    std::string basefilename;
    void addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax);
    long axsize0, axsize1;
    data_t bg;
  };
} // end namespace
#endif
