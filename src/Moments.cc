#include "../include/Moments.h"
#include "../include/MathHelper.h"

namespace shapelens {

  Moments::Moments() : tmv::Vector<data_t>(), N(0) { }

  Moments::Moments(int N_) : 
    tmv::Vector<data_t>(pyramid_num(N_+1), 0),
    N(N_) { }

  Moments& Moments::operator=(const tmv::GenVector<data_t>& v) {
    // FIXME: set N!
    tmv::Vector<data_t>::operator=(v);
    return *this;
  }

  Moments::Moments(const Object& obj, const WeightFunction& w, int N_, const Point<data_t>* centroid) :
    tmv::Vector<data_t>(pyramid_num(N_+1), 0), N(N_) {
    data_t w_, diff_x, diff_y, val;
    if (centroid == NULL)
      centroid = &obj.centroid;

    tmv::Vector<data_t> pow_x(N+1), pow_y(N+1);

    for (long i=0; i< obj.grid.size(); i++) {
      val = obj(i);
      w_ = w(obj.grid(i));
      diff_x = obj.grid(i,0) - (*centroid)(0);
      diff_y = obj.grid(i,1) - (*centroid)(1);

      for (int j=0; j <= N; j++) {
	pow_x(j) = pow_int(diff_x,j);
	pow_y(j) = pow_int(diff_y,j);
      }

      // mask bad pixels
      if (Config::CHECK_OBJECT && obj.segmentation.size() > 0)
	if (obj.segmentation(i) != 0 && obj.segmentation(i) != obj.id)
	  w_ = 0;

      for(int n=0; n <= N; n++)
	for(int m=0; m <= n; m++)
	  operator()(m,n-m) += pow_x(m) * pow_y(n-m) * val * w_;
    }
  }

  data_t& Moments::operator()(unsigned int px, unsigned int py) {
    return tmv::Vector<data_t>::operator()(pyramid_num(px+py)+py);
  }

  data_t Moments::operator()(unsigned int px, unsigned int py) const {
    return tmv::Vector<data_t>::operator()(pyramid_num(px+py)+py);
  }

  data_t& Moments::operator()(unsigned int i) {
    return tmv::Vector<data_t>::operator()(i);
  }

  data_t Moments::operator()(unsigned int i) const {
    return tmv::Vector<data_t>::operator()(i);
  }

  int Moments::getOrder() const {
    return N;
  }

  void Moments::setOrder(int N_) {
    N = N_;
    tmv::Vector<data_t>::resize(pyramid_num(N+1));
  }
  
  int Moments::getIndex(unsigned int px, unsigned int py) const {
    return pyramid_num(px+py)+py;
  }

  std::pair<int, int> Moments::getPowers(int i) const {
    // turns out that inverting getIndex() boils down to finding
    // the largest number 0 <= n <= i, for which t = sqrt(8*n + 1) is an integer
    // then, py = i - n and px = (-py -1 + t)/2
    // btw, n = pyramid_num(px + py)
    int n = i;
    data_t t;
    while (n >= 0) {
      t = sqrt(8*n + 1);
      if (t == int(t))
	break;
      n--;
    }
    int py = i - n;
    int px = (-2*py - 1 + int(t))/2;
    return std::pair<int, int> (px, py);
  }
  unsigned int Moments::pyramid_num(int n) const {
    return (n*(n+1))/2;
  }

} // end namespace

