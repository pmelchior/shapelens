#include "../include/Grid.h"

namespace shapelens {

  Grid::Grid() :
    start0(0), start1(0), N0(0), N1(0) {}

  Grid::Grid(int start0_, int start1_, unsigned int N0_, unsigned int N1_) :
    start0(start0_), start1(start1_), N0(N0_), N1(N1_) {}
  
  void Grid::setSize(int start0_, int start1_, unsigned int N0_, unsigned int N1_) {
    start0 = start0_;
    start1 = start1_;
    N0 = (int) N0_;
    N1 = (int) N1_; 
  }

  // this is the most often used access mode, so make it fast
  data_t Grid::operator() (long i, bool direction) const {
    if (ct.use_count() != 0) {
      Point<data_t> P(start0 + i%N0, start1 + i/N0);
      ct->transform(P);
      return P(direction);
    } else {
      int offset;
      if (direction) {
	offset = i/N0;
	return data_t(start1 + offset);
      }
      else {
	offset = i%N0;
	return data_t(start0 + offset);
      }
    }
  }

  Point<data_t> Grid::operator() (long i) const {
    if (ct.use_count() != 0) {
      Point<data_t> P(start0 + i%N0, start1 + i/N0);
      ct->transform(P);
      return P;
    }
    else
      return Point<data_t>(operator()(i,0),operator()(i,1));
  }

  void Grid::setWCS(const CoordinateTransformation& C) {
    if (typeid(C)!=typeid(nt)) // does not set a NullTransformation
      ct = C.clone();
  }

  void Grid::resetWCS() {
    ct.reset();
  }

  const CoordinateTransformation& Grid::getWCS() const {
    if (ct.use_count() != 0)
      return *ct;
    else
      return nt; // return valid WCS in all circumstances
  }

  data_t Grid::getScaleFactor() const {
    if (ct.use_count() != 0) {
      data_t scale_factor = getSupport().getArea();
      scale_factor /= getBoundingBox().getArea();
      return sqrt(scale_factor);
    } else
      return 1;
  }

  Polygon<data_t> Grid::getSupport() const {
    Rectangle<int> bb = getBoundingBox();
    std::list<Point<data_t> > pl;
    // lower-left in image coords
    Point<data_t> P(bb.ll(0), bb.ll(1));
    if (ct.use_count() != 0)
      ct->transform(P);
    pl.push_back(P);
    // lower-right
    P(0) = bb.tr(0);
    P(1) = bb.ll(1);
    if (ct.use_count() != 0)
      ct->transform(P);
    pl.push_back(P);
    // top-left
    P(0) = bb.ll(0);
    P(1) = bb.tr(1);
    if (ct.use_count() != 0)
      ct->transform(P);
    pl.push_back(P);
    // top-right
    P = bb.tr;
    if (ct.use_count() != 0)
      ct->transform(P);
    pl.push_back(P);
    Polygon<data_t> p(pl);
    // as support may have flipped edges, we need to remap it
    if (ct.use_count() != 0)
      p.remap();
      return p;
  }

  Rectangle<int> Grid::getBoundingBox() const {
    Rectangle<int> bbox;
    bbox.ll(0) = start0;
    bbox.ll(1) = start1;
    bbox.tr(0) = start0 + N0;
    bbox.tr(1) = start1 + N1;
    return bbox;
  }

  int Grid::getStartPosition(bool direction) const {
    if (direction)
      return start1;
    else
      return start0;
  }

  int Grid::getStopPosition(bool direction) const {
    if (direction)
      return start1 + N1;
    else
      return start0 + N0;
  }

  unsigned int Grid::getSize(bool direction) const {
    if (direction)
      return N1;
    else
      return N0;
  }

  unsigned long Grid::size() const {
    return N0*N1;
  }

  Point<int> Grid::getCoords(long pixel) const {
    return Point<int>(start0 + pixel%N0,start1 + pixel/N0);
  }

  Point<int> Grid::getCoords(const Point<data_t>& P) const {
    if (ct.use_count() != 0) {
      Point<data_t> P_ = P;
      ct->inverse_transform(P_);
      return Point<int>((int)floor(P_(0)),(int)floor(P_(1)));
    } else
      return Point<int>((int)floor(P(0)),(int)floor(P(1)));
  }

  long Grid::getPixel(const Point<int>& P) const {
    if (P(0) >= start0 && P(0) < start0 + int(N0) && P(1) >= start1 && P(1) < start1 + int(N1))
      return (P(0)-start0) + long(P(1)-start1)*int(N0);
    else
      return -1;
  }

  long Grid::getNeighborPixel(const Point<int>& P, unsigned char direction) const {
    long index;
    int x = P(0), y = P(1);
    switch(direction) {
    case 0: 
      // the pixel itself
      index = y*N0 + x;
      break;
    case 1: 
      if (x>0) index = y*N0 + x - 1; // left
      else index = -1;
      break;
    case 2:
      if (x<N0-1) index = y*N0 + x + 1;  // right neighbour
      else index = -1;
      break;
    case 3: 
      if (y>0 && x>0) index = (y-1)*N0 + x - 1;  // bottom left
      else index = -1;
      break;   
    case 4: 
      if (y>0) index = (y-1)*N0 + x;  // bottom
      else index = -1;
      break;
    case 5: 
      if (y>0 && x<N0-1) index = (y-1)*N0 + x + 1;  // bottom right
      else index = -1;
      break;  
    case 6: 
      if (y<N1-1 && x>0) index = (y+1)*N0 + x - 1;  // top left
      else index = -1;
      break;  
    case 7: 
      if (y<N1-1) index = (y+1)*N0 + x ;  // top
      else index = -1;
      break;
    case 8:
      if (y<N1-1 && x<N0-1) index = (y+1)*N0 + x + 1;  // top right
      else index = -1;
      break;
    }
    return index;
  }

  long Grid::getNeighborPixel(long pixel, unsigned char direction) const {
    if (direction==0)
      return pixel;
    else
      return getNeighborPixel(getCoords(pixel),direction);
  }

  std::ostream& operator<<(std::ostream& os, const Grid& grid) {
    os << "(" << grid.start0 << "/" << grid.start1 << ") .. ";
    os << "(" << grid.start0 + grid.N0 << "/" << grid.start1 + grid.N1 << ")";
    if (grid.ct.use_count() != 0)
      os << " with WCS";
    return os;
  }
} // end namespace
