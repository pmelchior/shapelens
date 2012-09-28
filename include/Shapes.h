#ifndef SHAPELENS_SHAPES_H
#define SHAPELENS_SHAPES_H

#include <list>
#include <stdexcept>
#include "ShapeLens.h"
#include "Point.h"
#include "CoordinateTransformation.h"

namespace shapelens {
  /// Rectangular patch.
  template <class T>
    class Rectangle {
  public:
    /// Lower-left boundary point.
    Point<T> ll;
    /// Top-right boundary point.
    Point<T> tr;
    /// Apply coordinate transformation.
    /// The result is again a Rectangle, defined by the bounding box of
    /// the transformed four points that initially defined the Rectangle.
    void apply(const CoordinateTransformation& C) {
      Point<T> lr(tr(0),ll(1)), tl(ll(0),tr(1));
      C.transform(ll);
      C.transform(tr);
      C.transform(lr);
      C.transform(tl);

      // find bounding box of all four corner points
      Point<T> min= ll;
      if (tr(0) < min(0))
	min(0) = tr(0);
      if (lr(0) < min(0))
	min(0) = lr(0);
      if (tl(0) < min(0))
	min(0) = tl(0);

      if (tr(1) < min(1))
	min(1) = tr(1);
      if (lr(1) < min(1))
	min(1) = lr(1);
      if (tl(1) < min(1))
	min(1) = tl(1);

      Point<T> max = tr;
      if (ll(0) > max(0))
	max(0) = ll(0);
      if (lr(0) > max(0))
	max(0) = lr(0);
      if (tl(0) > max(0))
	max(0) = tl(0);

      if (ll(1) > max(1))
	max(1) = ll(1);
      if (lr(1) > max(1))
	max(1) = lr(1);
      if (tl(1) > max(1))
	max(1) = tl(1);
      ll = min;
      tr = max;
    }

    /// Compute area of rectangle.
    inline T getArea() const {
      return (tr(0)-ll(0))*(tr(1)-ll(1));
      
    }
    /// Check whether Point \p P is inside Rectangle.
    /// Points on any edge are included.  
    inline bool contains(const Point<T>& P) const {
      if (P(0) >= ll(0) && P(0) <= tr(0) && P(1) >= ll(1) && P(1) <= tr(1))
	return true;
      else
	return false;
    }
    /// Add \p P to \p ll and \p tr.
    inline Rectangle<T>& operator+= (const Point<T>& P) {
      ll += P;
      tr += P;
      return *this;
    }
    /// Subtract \p P from \p ll and \p tr.
    inline Rectangle<T>& operator-= (const Point<T>& P) {
      ll -= P;
      tr -= P;
      return *this;
    }
    /// Ostream from Rectangle.
    friend std::ostream& operator<<(std::ostream& os, const Rectangle<T>& p) {
      os << p.ll << ".." << p.tr;
      return os;
    }
  };
  
  /// Edge of Polygon.
  template <class T>
    class Edge {
  public:
    /// First point of edge.
    Point<T> p1;
    /// Second point of edge.
    Point<T> p2;
    /// Move edge to next point \p p.
    /// Replace p1 by p2 and sets p2 to \p p.
    void moveTo(const Point<T>& p) {
      p1 = p2;
      p2 = p;
    }
    /// Apply coordinate transformation.
    void apply(const CoordinateTransformation& C) {
      C.transform(p1);
      C.transform(p2);
    }
  };

  /// Class for simple polygons.
  /// The class does not check whether the polygon is in fact simple or complex,
  /// (a complex polygon is one with crossing edges);
  /// this has to be done by the user when constructing it and can be checked
  /// with checkEdges().
  ///
  /// \b CAUTION: For integer-typed polygons, edge-artefacts can occur
  /// when testing contains().
  template <class T>
  class Polygon {
  public:
    /// Constructor.
    Polygon() {
    };
    /// Constructor from a list of Point.
    /// It is assumed that points form a closed polygon, that means
    /// the end of the edge chain is defined by the frist point in the list.
    Polygon(const std::list<Point<T> >& points) {
      if (points.size() < 3)
	throw std::invalid_argument("Polygon: List of edge points must contain at least 3 points");
      Edge<T> e;
      typename std::list<Point<T> >::const_iterator iter = points.begin();
      e.p1 = *iter;
      iter++;
      e.p2 = *iter;
      iter++;
      edges.push_back(e);
      while (iter != points.end()) {
	e.moveTo(*iter);
	edges.push_back(e);
	iter++;
      }
      e.moveTo(points.front());
      edges.push_back(e);
    }
    /// Add edge to polygon.
    /// The starting point of \p e must be identical to the end-point of the
    /// prior edge.
    void addEdge(const Edge<T>& e) {
      if (edges.size() == 0)
	edges.push_back(e);
      else if (e.p1 == edges.back().p2)
	edges.push_back(e);
      else
	throw std::invalid_argument("Polygon: Edges do not match");
    }
    /// Compute area of polygon.
    /// Returns \f$\frac{1}{2}\sum_{i=1}^{N}\bigl(x_i y_{i+1} - x_{i+1} y_i\bigr)\f$,
    /// where the \f$N\f$ edge points are given by \f$(x_i,y_i)\f$.
    data_t getArea() const {
      data_t a = 0;
      for (typename std::list<Edge<T> >::const_iterator iter = edges.begin(); iter != edges.end(); iter++) {
	a += (iter->p1(0)*iter->p2(1)) - (iter->p2(0)*iter->p1(1));
      }
      return fabs(0.5*a);
    }
    /// Checks whether \p p is inside the polygon.
    /// Uses crossing test: If a ray originating from \p into positive y-direction
    /// crosses an odd number of edges, its inside the polygon.
    bool contains(const Point<T>& p) const {
      unsigned int crossings = 0;
      for (typename std::list<Edge<T> >::const_iterator iter = edges.begin(); iter != edges.end(); iter++) {
	// x-cordinate of edge-points above and below p
	if ((iter->p1(0) - p(0))*(iter->p2(0) - p(0)) <= 0) {
	  // both y-coordinates above: must cross
	  if ((iter->p1(1) > p(1)) && (iter->p2(1) > p(1)))
	    crossings++;
	  // y-coords differ: crosses only if p is below edge line
	  else if ((iter->p1(1) - p(1))*(iter->p2(1) - p(1)) <= 0) {
	    data_t m = (iter->p2(1) - iter->p1(1))/(iter->p2(0) - iter->p1(0));
	    data_t dy = m*(p(0) - iter->p1(0));
	    if (iter->p1(1) + dy > p(1))
	      crossings++;
	  }
	}
      }
      return bool(crossings%2);
    }
    /// Return Rectangle which bounds the polygon.
    Rectangle<T> getBoundingBox() const {
      shapelens::Point<double> min,max;
      T mx,my,MX,MY;
      for (typename std::list<Edge<T> >::const_iterator iter = edges.begin(); iter != edges.end(); iter++) {
	mx = std::min(iter->p1(0),iter->p2(0));
	MX = std::max(iter->p1(0),iter->p2(0));
	my = std::min(iter->p1(1),iter->p2(1));
	MY = std::max(iter->p1(1),iter->p2(1));
	if (iter == edges.begin()) {
	  min(0) = mx;
	  min(1) = my;
	  max(0) = MX;
	  max(1) = MY;
	} else {
	  if (mx < min(0))
	    min(0) = mx;
	  if (my < min(1))
	    min(1) = my;
	  if (MX > max(0))
	    max(0) = MX;
	  if (MY > max(1))
	    max(1) = MY;
	}
      }
      shapelens::Rectangle<T> r;
      r.ll = min;
      r.tr = max;
      return r;
    }
    /// Checks whether polygonial chain is closed.
    /// Should be called after completion of the chain to ensure that first
    /// and last edge-points are identical.
    bool checkEdges() {
      if (edges.front().p1 == edges.back().p2)
	return true;
      else
	return false;
    }
    /// Apply coordinate transformation.
    void apply(const CoordinateTransformation& C) {
      for (typename std::list<Edge<T> >::iterator iter = edges.begin(); iter != edges.end(); iter++)
	iter->apply(C);
      //remap();
    }

    /// Remap the support after application of a CoordinateTransformation.
    /// As a CoordinateTransformation can create crossing edges (and thus
    /// complex polygons), the list of points must be reconnected with
    /// a new set of edges.\n
    /// \b CAUTION: remap() assumes the Polygon to be convex, i.e. the 
    /// chain of edges has minimal length.
    void remap() {
      // list of edge points
      std::list<Point<T> > points;
      for (typename std::list<Edge<T> >::iterator iter = edges.begin(); iter != edges.end(); iter++)
	points.push_back(iter->p1);
      edges.clear();
      // find shortest path thru list of edge points
      Edge<T> e;
      e.p1 = points.front();
      typename std::list<Point<T> >::iterator iter = points.begin(), nearest_iter;
      points.erase(iter);
      while (points.size() > 0) {
	data_t dist = distance(e.p1,points.front());
	iter = nearest_iter = points.begin();
	iter++;
	for(iter; iter != points.end(); iter++) {
	  data_t d = distance(e.p1,*iter);
	  if (d < dist) {
	    dist = d;
	    nearest_iter = iter;
	  }
	}
	e.p2 = *nearest_iter;
	addEdge(e);
	e.p1 = e.p2;
	points.erase(nearest_iter);
      }
      e.p2 = edges.front().p1;
      addEdge(e);
    }

    /// Add \p P to all Polygon edge-points
    Polygon<T>& operator+=(const Point<T>& P) {
      for (typename std::list<Edge<T> >::iterator iter = edges.begin(); iter != edges.end(); iter++) {
	iter->p1+=(P);
	iter->p2+=(P);
      }
      return *this;
    }
    /// Subtract \p P from all Polygon edge-points
    Polygon<T>& operator-=(const Point<T>& P) {
      for (typename std::list<Edge<T> >::iterator iter = edges.begin(); iter != edges.end(); iter++) {
	iter->p1-=(P);
	iter->p2-=(P);
      }
      return *this;
    }

    /// Ostream from Polygon
    friend std::ostream& operator<<(std::ostream& os, const Polygon<T>& p) {
      os << "//" << std::endl;
      for (typename std::list<Edge<T> >::const_iterator iter = p.edges.begin(); iter != p.edges.end(); iter++)
	os << iter->p1(0) << " " << iter->p1(1) << std::endl;
      return os;
    }
    /// Istream to Polygon
    friend std::istream& operator>>(std::istream& is, Polygon<T>& p) {
      if (is.bad())
	std::ios_base::failure("Polygon: Cannot read from istream");
      
      std::list<Point<T> > points;
      std::string delimiter;
      T x,y;
      is >> delimiter;
      while (is >> x >> y)
	points.push_back(Point<T>(x,y));
      if(points.size() > 0) {
        is.clear();
        p = Polygon<T>(points);
      }
      return is;
    }
    /// List of edges.
    std::list<Edge<T> > edges;
  private:
    data_t distance(const Point<T>& p1, const Point<T>& p2) {
      return (p1(0) - p2(0))*(p1(0) - p2(0)) + 
	(p1(1) - p2(1))*(p1(1) - p2(1));
    }
  };

  
} // end namespace

#endif
