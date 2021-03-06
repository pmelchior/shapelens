#include "../include/SExFrame.h"
#include "../include/MathHelper.h"
#include <set>
#include <vector>
#include <bitset>

using namespace std;
using namespace shapelens;
typedef unsigned int uint;

SExFrame::SExFrame (std::string datafile, std::string catfile, std::string segmapfile, std::string weightfile) : 
  catalog(catfile), basefilename(datafile), bg(0) {

  // open datafile
  fptr = FITS::openFile(datafile);
  long naxis[2];
  int status = 0;
  int dimensions;
  fits_get_img_dim (fptr, &dimensions, &status);
  if (dimensions == 2) {
    fits_get_img_size(fptr, dimensions, naxis, &status);
    axsize0 = naxis[0];
    axsize1 = naxis[1];
    grid.setSize(0,0,axsize0, axsize1);
    if (Config::USE_WCS) {
#ifdef HAS_WCSLIB
      grid.setWCS(WCSTransformation(fptr));
#else
      throw std::runtime_error("FITS: WCS usage requested, but HAS_WCSLIB not specified");
#endif
    }
  } else
    throw std::invalid_argument("SExFrame: FITS file " + datafile + " does not provide valid image!");
  
  // open weight and segmentation maps if specified
  if (weightfile != "") {
    fptr_w = FITS::openFile(weightfile);
    fits_get_img_dim (fptr_w, &dimensions, &status);
    if (dimensions != 2)
      throw std::invalid_argument("SExFrame: FITS file " + weightfile + " does not provide valid image!");
  }
  else
    fptr_w = NULL;
  if (segmapfile != "") {
    fptr_s = FITS::openFile(segmapfile);
    fits_get_img_dim (fptr_s, &dimensions, &status);
    if (dimensions != 2)
      throw std::invalid_argument("SExFrame: FITS file " + segmapfile + " does not provide valid image!");
  }
  else
    fptr_s = NULL;
}

SExFrame::~SExFrame() {
  FITS::closeFile(fptr);
  if (fptr_w != NULL)
    FITS::closeFile(fptr_w);
  if (fptr_s != NULL)
    FITS::closeFile(fptr_s);
}

unsigned long SExFrame::getNumberOfObjects() {
  return catalog.size();
}

void SExFrame::fillObject(Object& O, Catalog::const_iterator& catiter) {
  if (catiter != catalog.end()) {
    O.id = catiter->first;
    int xmin, xmax, ymin, ymax;
    xmin = catiter->second.XMIN;
    xmax = catiter->second.XMAX;
    ymin = catiter->second.YMIN;
    ymax = catiter->second.YMAX;

    O.history << "# Extracting Object " << (*catiter).first;
    O.history << " found in the area (" << xmin << "/" << ymin << ") to (";
    O.history << xmax << "/" << ymax << ")" << std::endl;
  
    // check if outer sizes of the object are identical to the image
    // boundary, since then the objects is cutted 
    bool cutflag = 0;
    if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
      O.history << "# Object cut off at the image boundary!" << endl;
      cutflag = 1;
    }
  
    addFrameBorder(Config::ADD_BORDER, xmin,xmax,ymin,ymax);
    O.history << "# Extending the area around object to (" << xmin << "/" << ymin << ") to (";
    O.history << xmax << "/" << ymax << ")" << std::endl;

    // prepare containers
    O.resize((xmax-xmin)*(ymax-ymin));
    // Grid will be changed but not shifted (all pixels stay at their position)
    O.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
    O.centroid = Point<data_t>(catiter->second.XCENTROID,
			       catiter->second.YCENTROID);
    if (Config::USE_WCS)
      O.grid.setWCS(grid.getWCS());

    // copy pixel data
    int status = 0;
    data_t nullval = 0;
    int anynull = 0;
    long firstpix[2] = {xmin+1,ymin+1}, lastpix[2] = {xmax, ymax}, inc[2] = {1,1};

    // read image, and subtract background if requested
    fits_read_subset(fptr, FITS::getDataType(data_t(0)), firstpix, lastpix, inc, &nullval, O.ptr(), &anynull, &status);
    if (bg != 0)
      O.addToAll(-bg);

    // read weightmap and segmentation map if provided
    if (fptr_w != NULL) {
      O.weight.resize((xmax-xmin)*(ymax-ymin));
      O.weight.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
      if (Config::USE_WCS)
	O.weight.grid.setWCS(grid.getWCS());
      fits_read_subset(fptr_w, FITS::getDataType(data_t(0)), firstpix, lastpix, inc, &nullval, O.weight.ptr(), &anynull, &status);
    }
    // even if there is not segmentation map, 
    // create one to store bad pixels
    if (fptr_s != NULL || Config::CHECK_OBJECT) {
      O.segmentation.resize((xmax-xmin)*(ymax-ymin));
      O.segmentation.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
      if (Config::USE_WCS)
	O.segmentation.grid.setWCS(grid.getWCS());

      if (fptr_s != NULL)
	fits_read_subset(fptr_s, FITS::getDataType(int(0)), firstpix, lastpix, inc, &nullval, O.segmentation.ptr(), &anynull, &status);
      else
	O.segmentation.setZero();
    }

    // check image for bad pixels
    if (Config::CHECK_OBJECT) {
      if (fptr_w != NULL) {
	for (int i = 0; i < O.size(); i++) {
	  if (O.weight(i) <= 0)
	    O.segmentation(i) = -2;
	}
      }
    }

    // Fill other quantities into Object
    O.flags = std::bitset<8>(catiter->second.FLAGS);
    O.basefilename = basefilename;
  }
  else {
    std::ostringstream mess;
    mess << "SExFrame: Object " << O.id << " does not exist!";
    throw std::invalid_argument(mess.str());
  }
}

// now extend to region around the object by a factor > 0, minimum 16 pixels
void SExFrame::addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
  if (factor > 0) {
    long int xrange, yrange, xborder, yborder;
    xrange = xmax - xmin;
    yrange = ymax - ymin;
    if (xrange%2 == 1) {
      xmax++;
      xrange++;
    }
    if (yrange%2 == 1) {
      ymax++;
      yrange++;
    }
    // make the object frame square to avoid anisotropies
    if (xrange < yrange) {
      yborder = std::max((int)floor(yrange*factor), 16);
      xborder = yborder + (yrange - xrange)/2;
    } else {
      xborder = std::max((int)floor(xrange*factor), 16);
      yborder = xborder + (xrange - yrange)/2;
    }
    // limit corner points to within frame
    xmin = std::max((long int) 0, xmin - xborder);
    xmax = std::min(axsize0, xmax + xborder);
    ymin = std::max((long int) 0, ymin - yborder);
    ymax = std::min(axsize1, ymax + yborder);
  }
}

void SExFrame::subtractBackground(data_t bg_) {
  bg = bg_;
}
