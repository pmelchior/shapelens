#ifndef SHAPELENS_PROPERTY_H
#define SHAPELENS_PROPERTY_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <boost/variant.hpp>
#include "ShapeLens.h"

namespace shapelens {
/// Allowed types for usage in Property.
typedef boost::variant<std::string,data_t,int,std::vector<int>,std::vector<data_t>,std::vector<std::string> > variant_t;

/// Flexible container class.
/// The class provides means of manipulating data content from which neighter names, types or
/// values have to be known <i>a priori</i>.
/// It is very convenient to annotate individual objects which carry a
/// Property structure, or to provide configuration information since serialization via 
/// read() and write() is guaranteed.\n\n
/// The class is derived from \p std::map<std::string, variant_t>,
/// where \p variant_t can be either
/// - \p int
/// - \p data_t
/// - \p std::string
/// - \p std::vector<int>
/// - \p std::vector<data_t>
/// - \p std::vector<std::string>
/// 
/// This allows flexible storage containers, e.g.
/// \code
///  Property p;
///  p["e1"] = 0.1;
///  p["index"] = 1;
///  p["info"] = "Some information...";
///  p["vi"] = std::vector<int>(5,42);
/// \endcode
/// The implementation relies on the concept of \p boost::variant, details can be
/// found at http://www.boost.org/doc/html/variant.html.\n\n
/// In comparisons with or assignment from a \p variant_t, the type has to be known
/// \code
/// std::cout << boost::get<int>(p1["index"]) << std::endl;
/// std::cout << boost::get<data_t>(p1["e1"]) << std::endl;
/// \endcode
/// or checked
/// \code
/// for (Property::iterator iter = p1.begin(); iter != p1.end(); iter++) {
///    int i;
///    if (iter->second.type() == typeid(i))
///      std::cout << iter->first << "\t" << boost::get<int>(iter->second) << std::endl;
///    // ... similar checks for other data types ...
///  }
/// \endcode
/// Alternatively, use the concept of \p boost::static_visitor (described as part
/// of the \p boost::variant documentation). 

class Property : public std::map<std::string, variant_t> {
 public:
  /// Constructor.
  Property();
  /// Access operator.
  /// Throws \p std::invalid_argument if \p key is not element of the map.
  variant_t& operator[](const std::string& key);
  /// Const access operator.
  /// Throws \p std::invalid_argument if \p key is not element of the map.
  const variant_t& operator[](const std::string& key) const ;
  /// Convenient getter function.
  /// Throws <tt>std::invalid_argument</tt> if \p key is not present in Property
  /// and <tt>boost::bad_get</tt> if \p key has different type than \p T.
  template <class T>
    void get(const std::string& key, T& val) const {
    val = boost::get<T>(Property::operator[](key));
  }
  /// Conventient setter function.
  template <class T>
    void set(const std::string& key, const T& val) {
    Property::operator[](key) = val;
  }
  /// Write contents of property to \p out.
  /// If the entries should be written with a different separator than \p std::endl,
  /// specify \p linesep. This should only be used in cases where the line separator
  /// needs to be masked in order to be written to a certain format.\n
  /// \b CAUTION: 
  /// - \p linesep must not be either \p "\t" or \p "," as they are
  /// needed for the serialization.
  /// - The elements of \p std::vector<std::string> must not contain \p "\t" or \p ",".
  /// - The accuracy of \p out may truncate number of type \p data_t. The standard accuracy
  /// is 6 digits and may be overwritten by
  /// \code
  /// #include <iomanip>
  /// out << std::setprecision(12);
  /// \endcode
  /// prior to a call to write().
  void write(std::ostream& out, std::string linesep = std::string()) const;
  /// Read contents of property from \p out.
  void read(std::istream& in);
};
} // end namespace
#endif
