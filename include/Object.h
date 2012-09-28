#ifndef SHAPELENS_OBJECT_H
#define SHAPELENS_OBJECT_H

#include "ShapeLens.h"
#include "Image.h"
#include "SegmentationMap.h"
#include <bitset>

namespace shapelens {

  /// Central object representing class.
  /// The purpose of this class is to facilitate the exchange of object-related
  /// information between different parts of the code. 
  /// By definition, an object is a significant peak of the 
  /// brightness distribution, identified by some image processing code such 
  /// as SExtractor. See SExFrame for details on how to load Object information.
  ///
  /// A instance of this class describes an object in real (pixel) space.
  /// It consists of the pixel values and the Grid on which it is defined.
  /// In addition, it can store information from the preceding processing steps:
  /// - centroid position
  /// - detection flags
  /// - processing history
  /// - weight map or noise_rms to describe the noise properties
  /// - segmentation map
  
  class Object : public Image<data_t> {
  public:
    /// Constructor.
    Object();
    /// Argumented constructor for loading an object from a Fits file.
    /// The Fits file shold have been created by Object::save().
    Object (std::string fitsfile);
    /// Copy constructor from base class.
    Object (const Image<data_t>& base);
    /// Save the object information in a Fits file.
    /// Data and SegmentationMap will go to pHDU and 1st extHDU, respectively. 
    /// All other information goes to the pHDU header.
    /// If a weight map is provided, these will be stored in the 
    /// extension \p WEIGTH.
    void save(std::string fitsfile);
    
    /// The \p id of the Object.
    unsigned long id;
    /// The weight (inverse variance) map in the region of this object.
    /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
    Image<data_t> weight;
    /// The position of the object's centroid.
    Point<data_t> centroid;
    /// The detection flags.
    std::bitset<8> flags;
    /// The noise mean \f$\mu_n\f$.
    data_t noise_mean;
    /// The noise RMS \f$\sigma_n\f$.
    data_t noise_rms;
    /// The segmentation map.
    SegmentationMap segmentation;
  };
} // end namespace
#endif
