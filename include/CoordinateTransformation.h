#ifndef SHAPELENS_COORDINATETRANSFORMATION_H
#define SHAPELENS_COORDINATETRANSFORMATION_H

#include "ShapeLens.h"
#include "Point.h"
#include <boost/shared_ptr.hpp>
#include <list>

namespace shapelens {
  /// Base class for coordinate transformations in 2D.
  class CoordinateTransformation {
  public:
    /// Destructor.
    ~CoordinateTransformation();
    /// Chain a succeding transformation.
    /// The new transformation will be made of the product
    /// <tt>(*this) -> (*this) * C</tt>
    void operator*=(const CoordinateTransformation& C);
    /// Invert the transformation.
    void invert();
    /// Apply transformation to \p P.
    void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    void inverse_transform(Point<data_t>& P) const;
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const = 0;
  protected:
    /// Apply all transformations from the stack.
    void stack_transform(Point<data_t>& P) const;
    /// Apply all inverse transformation from the stack in reverse order.
    void stack_inverse_transform(Point<data_t>& P) const;
    /// Stack of transformation to be applied after another
    std::list<boost::shared_ptr<CoordinateTransformation> > stack;
    CoordinateTransformation();
  private:
    bool inverted;
    // actual virtual functions for forward and inverse transformation
    virtual void f(Point<data_t>& P) const = 0;
    virtual void f_1(Point<data_t>& P) const = 0;
  };

  /// Empty/Null transformations in 2D.
  /// This transformations will leave all points unchanged.
  class NullTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    NullTransformation();
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    /// Apply transformation to \p P.
    virtual void f(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void f_1(Point<data_t>& P) const;
  };

  /// Class for rescaling transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=s\cdot P\f$ with some scalar value \f$s\f$.
  class ScalarTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    ScalarTransformation(data_t scale);
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    data_t s;
    virtual void f(Point<data_t>& P) const;
    virtual void f_1(Point<data_t>& P) const;
  };

  /// Class for translation transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=P+ \Delta P\f$ with some shift direction \f$\Delta P\f$.
  class ShiftTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    ShiftTransformation(const Point<data_t>& dP_);
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    Point<data_t> dP;
    virtual void f(Point<data_t>& P) const;
    virtual void f_1(Point<data_t>& P) const;
  };

  /// Class for linear transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f[P^\prime= M\cdot\ P,\f]
  /// where \f$M\f$ is a 2x2 matrix 
  class LinearTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    LinearTransformation(const tmv::Matrix<data_t>& M_);
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    tmv::Matrix<data_t> M, M_1;
    virtual void f(Point<data_t>& P) const;
    virtual void f_1(Point<data_t>& P) const;
  };

  /// Class for lensing transformations in 2D.
  class LensingTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    LensingTransformation(data_t kappa_, std::complex<data_t> gamma);
    /// Constructor for 2nd order lensing transformation.
    LensingTransformation(data_t kappa_, std::complex<data_t> gamma, std::complex<data_t> F, std::complex<data_t> G);
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    bool flex;
    data_t kappa, gamma1, gamma2, D111, D112, D121, D122, D211, D212, D221, D222;
    virtual void f(Point<data_t>& P) const;
    virtual void f_1(Point<data_t>& P) const;    
  };
} // end namespace

#endif
