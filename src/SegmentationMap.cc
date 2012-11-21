#include "../include/SegmentationMap.h"

using namespace shapelens;
using namespace std;

SegmentationMap::SegmentationMap() : 
Image<int>() {
}

SegmentationMap::SegmentationMap(string segMapFile) : 
Image<int>(segMapFile) {
}

int SegmentationMap::getNumberOfObjects() {
  set<int> objects;
  for (unsigned long i=0; i<SegmentationMap::size(); i++)
    // if pixel belongs to object
    if (SegmentationMap::operator()(i) != 0)
      // inserts object number in set if it is not present yet
      objects.insert(SegmentationMap::operator()(i));
  return objects.size();
}


// find list of pixels due to given object from segmentation map 
void SegmentationMap::findObjectPixels(std::set<unsigned long>& pixelset, unsigned long objectnr, long xmin, long xmax, long ymin, long ymax) {
  long axsize0 = SegmentationMap::grid.getSize(0), axsize1 = SegmentationMap::grid.getSize(1);
  pixelset.clear();
  Point<int> P;
  for (P(1) = std::max<int>(ymin,0); P(1) <= std::min<int>(ymax,axsize1-1); P(1)++) {
    for (P(0) = std::max<int>(xmin,0); P(0) <= std::min<int>(xmax,axsize0-1); P(0)++) {
      int j = SegmentationMap::grid.getPixel(P);
      if (SegmentationMap::operator()(j) == objectnr)
	pixelset.insert(j);
    }
  }
}

// draws a rectangular frame with the given limits in the segmenation map
void SegmentationMap::setSegmentBorder(int tag, long xmin, long xmax, long ymin, long ymax) {
  long axsize0 = SegmentationMap::grid.getSize(0), axsize1 = SegmentationMap::grid.getSize(1);
  // check if corners are within frame
  if (xmin<0) xmin=0;
  if (xmax>=axsize0) xmax=axsize0-1;
  if (ymin<0) ymin=0;
  if (ymax>=axsize1) xmax=axsize1-1;
  // low border
  for (long i=(xmin+ymin*axsize0); i<=(xmax+ymin*axsize0); i++)
    if (SegmentationMap::operator()(i) == 0)
      SegmentationMap::operator()(i) = tag;
  // high border
  for (long i=(xmin+ymax*axsize0); i<=(xmax+ymax*axsize0); i++)
    if (SegmentationMap::operator()(i) == 0)
      SegmentationMap::operator()(i) = tag;
  // left border
  for (long i=(xmin+(ymin+1)*axsize0);i<=(xmin+(ymax-1)*axsize0); i+=axsize0)
    if (SegmentationMap::operator()(i) == 0)
      SegmentationMap::operator()(i) = tag;
  // right border
  for (long i=(xmax+(ymin+1)*axsize0);i<=(xmax+(ymax-1)*axsize0); i+=axsize0)
    if (SegmentationMap::operator()(i) == 0)
      SegmentationMap::operator()(i) = tag;
}
