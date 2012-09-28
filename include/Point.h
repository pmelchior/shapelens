#ifndef SHAPELENS_POINT2D_H
#define SHAPELENS_POINT2D_H

#include <TMV_Small.h>

namespace shapelens {

/// Templated 2D Point class.

  template <class T>
    class Point : public tmv::SmallVector<T, 2> {
  public:
    /// Constructor.
    Point () : tmv::SmallVector<T, 2>(0) {}
    /// Constructor with given coordinates.
    template <class R> Point (R x0, R x1) : tmv::SmallVector<T, 2>() {
      tmv::SmallVector<T, 2>::operator()(0) = x0;
      tmv::SmallVector<T, 2>::operator()(1) = x1;
    }
    /// Copy constructor from base class.
    Point (const tmv::SmallVector<T, 2>& r) : tmv::SmallVector<T, 2>(r) { }

    /// Copy operator.
    Point<T>* operator=(const Point<T>& P) {
      tmv::SmallVector<T, 2>::operator=(P);
      return this;
    }
    
    /// Copy operator from base class
    Point<T>* operator=(const tmv::SmallVector<T, 2>& P) {
      tmv::SmallVector<T, 2>::operator=(P);
      return this;
    }
    
    /// Comparison operator.
    /// This is important for sorted containers of the STL.
    /// To ensure efficient lookups for image-type data, the points are
    /// ordered according to their 2nd dimension
    bool operator<(const Point<T>& b) const {
      if (tmv::SmallVector<T, 2>::operator()(1) < b(1) || (tmv::SmallVector<T, 2>::operator()(1) == b(1) && tmv::SmallVector<T, 2>::operator()(0) < b(0)))
	return true;
      else
	return false;
    }

    /// Equality operator.
    bool operator==(const Point<T>& b) const {
      if (tmv::SmallVector<T, 2>::operator()(0) == b(0) && tmv::SmallVector<T, 2>::operator()(1) == b(1))
	return true;
      else
	return false;
    }
    /// Return pointer to data storage.
    T* c_array() {
      return tmv::SmallVector<T, 2>::ptr();
    }
    /// Return pointer to data storage.
    const T* c_array() const {
      return tmv::SmallVector<T, 2>::cprt();
    }
  };
} // end namespace
#endif
