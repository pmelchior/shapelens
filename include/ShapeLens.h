#ifndef SHAPELENS_H
#define SHAPELENS_H

#include <complex>
#include <iostream>
#include <TMV.h>

/// Namespace for shapelens
namespace shapelens {
  
  /// Floating point format for all data used in shapelens.
  typedef double data_t;

  /// Complex datatype
  using std::complex;

  /// Class for configuration parameters.
  /// This class stores configuration parameters determining the behaviour of 
  /// several aspect of the data flow.\n
  /// In order to be effective, the parameters have to be set/changed before 
  /// the respective method is called, which often means that the parameters
  /// are set at the beginning of the script, either by direct specfication
  /// or from a configuration file:
  /// \code
  /// int main() {
  ///   Config sconf("shapelens.conf");
  ///   ...
  ///   Config::VERBOSITY = true;
  ///   ...
  /// }
  /// \endcode
  /// In this example, the configuration is initialized form a file, but then
  /// altered later, here to switch on the verbose output.

  class Config {
  public:
    /// Default constructor.
    Config();
    /// Argumented constructor.
    /// <tt>filename</tt> is expected to contain configuration parameters
    /// in ASCII format, one keyword/value pair per line, separated by a single or multiple
    /// tab character(s):
    /// \code
    /// CHECK_OBJECT 1
    /// ADD_BORDER   0.25
    /// VERBOSITY    1
    /// \endcode
    Config(std::string filename);
    /// Whether the History should be printed to stdout during data processing,
    /// default = false.
    static bool VERBOSITY;
    /// Whether the segmented objects should be checked for bad pixels and
    /// boundary truncation, default = true.
    static bool CHECK_OBJECT;
    /// The amount by which the objects area is enlarged on each side relative to the
    /// object's area in the segmentation map, default = 0.5
    /// (50% increase in size an all sides).
    static data_t ADD_BORDER;
    /// Whether coordinates are treated as sky coordinates, default = false;
    /// WCS information is read from FITS file header.
    static bool USE_WCS;
  };
} // end namespace
#endif
