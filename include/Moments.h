#ifndef SHAPELENS_MOMENTS_H
#define SHAPELENS_MOMENTS_H

#include "ShapeLens.h"
#include "Object.h"
#include "WeightFunction.h"

namespace shapelens {
  /// Container class for arbitrary moments of a two-dimensional quantity.
  /// Access is provided by an index or by the moment powers in 
  /// \f$x\f$ and \f$y\f$.
  /// All moments are computed in one sweep over the data.
  class Moments : public tmv::Vector<data_t> {
  public:
    /// Default constructor
    Moments();
    /// Constructor for moments up to order \p N.
    Moments(int N);
    /// Copy operator.
    Moments& operator=(const tmv::GenVector<data_t>& m);
    /// Constructor for moments up to order \p N.
    /// The moments are populated from \p obj. If \p centroid is set, it will
    /// be used to center the moments, otherwise Object::centroid is used.
    Moments(const Object& obj, const WeightFunction& w, int N, const Point<data_t>* centroid = NULL);
    /// Access operator for vector index.
    data_t& operator()(unsigned int i);
    /// Access operator for vector index.
    data_t operator()(unsigned int i) const;
    /// Access operator for \f$\langle x^{p_x}\, y^{p_y}\rangle\f$.
    data_t& operator()(unsigned int px, unsigned int py);
    /// Access operator for \f$\langle x^{p_x}\, y^{p_y}\rangle\f$.
    data_t operator()(unsigned int px, unsigned int py) const;
    /// Get maximum moment order.
    int getOrder() const;
    /// Set the maximum order to \p N.
    void setOrder(int N);
    /// Get vector index of moment \f$\langle x^{p_x}\, y^{p_y}\rangle\f$ from 
    /// the powers.
    int getIndex(unsigned int px, unsigned int py) const;
    /// Get moment powers from the vector index \p i.
    std::pair<int, int> getPowers(int i) const;
  private:
    unsigned int pyramid_num(int n) const;
    int N;
  };



} // end namespace
#endif
