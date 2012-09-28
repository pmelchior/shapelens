#include "../include/WeightFunction.h"
#include "../include/MathHelper.h"
#include <stdexcept>

namespace shapelens {

  FlatWeightFunction::FlatWeightFunction() {
  }
  
  data_t FlatWeightFunction::operator()(const Point<data_t>& P) const {
    return 1;
  }

  LocalWeightFunction::LocalWeightFunction(const Point<data_t>& centroid) :
    C(centroid) {
  }

  void LocalWeightFunction::setCentroid(const Point<data_t>& centroid) {
    C = centroid;
  }
  const Point<data_t>& LocalWeightFunction::getCentroid() const {
    return C;
  }  

  GaussianWeightFunction::GaussianWeightFunction(data_t scale_, const Point<data_t>& centroid):
    LocalWeightFunction(centroid),
    scale(scale_), n(0), sigma2(scale_*scale_) {
    fptr = &GaussianWeightFunction::Gauss;
  }

  data_t  GaussianWeightFunction::Gauss(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return Gauss(r);
  }
  data_t  GaussianWeightFunction::Gauss_(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return Gauss(r)*(-r/sigma2);
  }
  data_t  GaussianWeightFunction::Gauss__(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return Gauss(r)*(pow2(r/sigma2) - 1./sigma2);
  }
  data_t  GaussianWeightFunction::Gauss___(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return Gauss(r)*(pow3(-r/sigma2) + 3*r/(sigma2*sigma2));
  }
  data_t  GaussianWeightFunction::Gauss_2(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return -Gauss(r)/(2*sigma2);
  }
  data_t  GaussianWeightFunction::Gauss__2(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return Gauss(r)/pow2(2*sigma2);
  }
  data_t  GaussianWeightFunction::Gauss___2(const Point<data_t>& P) const {
    data_t r = sqrt(pow2(P(0)-C(0)) + pow2(P(1)-C(1)));
    return -Gauss(r)/pow3(2*sigma2);
  }
  data_t  GaussianWeightFunction::Gauss(data_t r) const {
    return exp(-r*r/(2*sigma2));
  }

  void GaussianWeightFunction::setDerivative(int n_) {
    n = n_;
    switch (n) {
      case 0: fptr = &GaussianWeightFunction::Gauss; break;
      case 1: fptr = &GaussianWeightFunction::Gauss_; break;
      case -1: fptr = &GaussianWeightFunction::Gauss_2; break;
      case 2: fptr = &GaussianWeightFunction::Gauss__; break;
      case -2: fptr = &GaussianWeightFunction::Gauss__2; break;
      case 3: fptr = &GaussianWeightFunction::Gauss___; break;
      case -3: fptr = &GaussianWeightFunction::Gauss___2; break;
      default: throw std::invalid_argument("WeightFunction: derivative of Gaussian invalid");
    }
  }
  int GaussianWeightFunction::getDerivative() const {
    return n;
  }

  data_t GaussianWeightFunction::getScale() const {
    return scale;
  }
  void GaussianWeightFunction::setScale(data_t scale_) {
    scale = scale_;
    sigma2 = scale*scale;
  }

  PowerLawWeightFunction::PowerLawWeightFunction(const Point<data_t>& centroid, data_t index) :
    LocalWeightFunction(centroid), n(index) {
  }
  
  data_t PowerLawWeightFunction::operator() (const Point<data_t>& P) const {
    data_t r = sqrt((P(0) - C(0))*(P(0) - C(0)) + 
		    (P(1) - C(1))*(P(1) - C(1)));
    return pow(r,n);
  }
}
