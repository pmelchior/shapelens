#ifndef SHAPELENS_WEIGHTFUNCTION_H
#define SHAPELENS_WEIGHTFUNCTION_H

#include "ShapeLens.h"
#include "Point.h"

namespace shapelens {
  /// Abstract base class for weight functions.
  class WeightFunction {
  public:
    /// Get value of weight function at position \p P.
    virtual data_t operator()(const Point<data_t>& P) const = 0;
  };

  /// Default implementation of flat weight function.
  class FlatWeightFunction : public WeightFunction {
  public:
    /// Constructor.
    FlatWeightFunction();
    /// Get value of weight function at position \p P.
    /// Returns 1.
    virtual data_t operator()(const Point<data_t>& P) const;
  };

  /// Localized weight function.
  class LocalWeightFunction : public WeightFunction {
  public:
    /// Constructor.
    LocalWeightFunction(const Point<data_t>& P);
    /// Set the centroid.
    virtual void setCentroid(const Point<data_t>& centroid);
    /// Get centroid position.
    virtual const Point<data_t>& getCentroid() const;
    /// Get value of weight function at position \p P.
    virtual data_t operator()(const Point<data_t>& P) const = 0;
    protected:
    /// Centroid/reference position.
    Point<data_t> C;  
  };

  /// Gaussian weight function.
  /// Returns \f$W(x) = \exp\Bigl(\frac{-x^{\prime 2}}{2s^s}\Bigr)\f$, i.e.
  /// the weight are drawn from a circular Gaussian of scale \f$s\f$,
  /// centered at some predefined centroid position.\n
  /// Derivatives (w.r.t. \f$s\f$ or \f$s^2\f$) can be accessed by
  /// choosing the repective derivative order with setDerivative().
  class GaussianWeightFunction : public LocalWeightFunction {
  public:
    /// Constructor.
    GaussianWeightFunction(data_t scale, const Point<data_t>& centroid);
    /// Get value of weight function.
    virtual data_t operator()(const Point<data_t>& P) const;
    /// Get scale.
    data_t getScale() const;
    /// Set scale.
    void setScale(data_t scale);
    /// Set derivative of weight function.
    /// Positive \p n denotes derivatives w.r.t. \f$s\f$, 
    /// negative w.r.t. \f$s^2\f$. Orders up to 3 are implemented, order 0
    /// recovers the original weight function.
    void setDerivative(int n);
    /// Get the current derivative.
    int getDerivative() const;
  protected:
    int n;
    data_t scale, sigma2;
    data_t (GaussianWeightFunction::*fptr) (const Point<data_t>&) const;
    data_t Gauss(data_t r) const;
    data_t Gauss(const Point<data_t>& P) const;
    data_t Gauss_(const Point<data_t>& P) const;
    data_t Gauss__(const Point<data_t>& P) const;
    data_t Gauss___(const Point<data_t>& P) const;
    data_t Gauss_2(const Point<data_t>& P) const;
    data_t Gauss__2(const Point<data_t>& P) const;
    data_t Gauss___2(const Point<data_t>& P) const;
  };

  inline data_t GaussianWeightFunction::operator()(const Point<data_t>& P) const {
    return (*this.*fptr)(P);
  }

  /// Weight function for power-law distance weighting.
  class PowerLawWeightFunction : public LocalWeightFunction {
  public:
    /// Constructor
    PowerLawWeightFunction(const Point<data_t>& centroid, data_t index);
    /// Return \f$r^n\f$.
    /// \f$ r\f$ is the Euclidean distance between \p P and the reference point
    /// and \f$ n\f$ the power-law index.
    virtual data_t operator() (const Point<data_t>& P) const;
  private:
    data_t n;
  };

} // end namespace
#endif
