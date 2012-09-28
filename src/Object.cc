#include "../include/Object.h"
#include <sstream>

using namespace shapelens;

Object::Object() : Image<data_t>(), segmentation() {
  id = 0;
  flags = 0;
  centroid(0) = centroid(1) = 0;
  noise_mean = noise_rms = 0;
}

Object::Object (const Image<data_t>& base) : Image<data_t>(base), segmentation()  {
  id = 0;
  flags = 0;
  centroid(0) = centroid(1) = 0;
  noise_mean = noise_rms = 0;
}

Object::Object(std::string objfile) : Image<data_t>(), segmentation() {
  int status, nkeys, keypos, hdutype;
  char card[FLEN_CARD];
  char comment[FLEN_CARD];

  history << "# Loading object from Fits file " << objfile << ":" << std::endl;
  fitsfile* fptr = FITS::openFile(objfile);

  // reading objects pixel data
  history << "# Reading object's pixel data";
  FITS::readImage(fptr,*this);
  
  // recover object information from header keywords
  FITS::readKeyword(fptr,"BASEFILE",Image<data_t>::basefilename);
  FITS::readKeyword(fptr,"ID",id);
  int xmin,ymin;
  FITS::readKeyword(fptr,"XMIN",xmin);
  FITS::readKeyword(fptr,"YMIN",ymin);
  Image<data_t>::grid = Grid(xmin,ymin,grid.getSize(0),grid.getSize(1));
  std::complex<data_t> xc;
  FITS::readKeyword(fptr,"CENTROID",xc);
  centroid(0) = real(xc);
  centroid(1) = imag(xc);
  FITS::readKeyword(fptr,"BG_MEAN",noise_mean);
  FITS::readKeyword(fptr,"BG_RMS",noise_rms);
  unsigned long f;
  FITS::readKeyword(fptr,"FLAG",f);
  flags = std::bitset<8>(f);
  
  // read history
  std::string hstr;
  FITS::readKeyCards(fptr,"HISTORY",hstr);

  // check whether grid has same size as object
  if (Object::size() != Image<data_t>::grid.size())
    throw std::invalid_argument("Object: Grid size from header keywords wrong!");
  
  history << ", segmentation map";
  // move to 1st extHDU for the segmentation map
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  FITS::readImage(fptr, segmentation);

  // check if there is 2nd extHDU: the weight map or correlation
  if (!fits_movabs_hdu(fptr, 3, &hdutype, &status)) {
    std::string extname;
    FITS::readKeyword(fptr, "EXTNAME", extname);
    if (extname == "WEIGHT") {
      history << " and weight map";
      FITS::readImage(fptr, weight);
    } else if (extname == "CORRELATFITSN") {
      history << " and correlation function";
      tmv::Matrix<data_t> corr;
      FITS::readImage(fptr, corr);
    }
  }
  history << std::endl;

  // append pHDUs history
  history.setSilent();
  history << hstr;
  history.unsetSilent();

  FITS::closeFile(fptr);
}

void Object::save(std::string filename) {
  // write pixel data
  fitsfile *outfptr = FITS::createFile(filename);
  FITS::writeImage(outfptr,*this);

  // add object information to header
  FITS::updateKeyword(outfptr,"BASEFILE",Image<data_t>::basefilename,"name of source file");
  FITS::updateKeyword(outfptr,"ID",id,"object id");
  FITS::updateKeyword(outfptr,"XMIN",Image<data_t>::grid.getStartPosition(0),"min(X) in image pixels");
  FITS::updateKeyword(outfptr,"YMIN",Image<data_t>::grid.getStartPosition(1),"min(Y) in image pixels");
  std::complex<data_t> xc(centroid(0),centroid(1));
  FITS::updateKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");
  FITS::updateKeyword(outfptr,"BG_MEAN",noise_mean,"mean of background noise");
  FITS::updateKeyword(outfptr,"BG_RMS",noise_rms,"rms of background noise");
  FITS::updateKeyword(outfptr,"FLAG",flags.to_ulong(),"extraction flags");
  FITS::appendHistory(outfptr,Image<data_t>::history.str());

  // save segmentation
  if (segmentation.size() != 0) {
    FITS::writeImage(outfptr,segmentation,"SEGMAP");
    FITS::appendHistory(outfptr,segmentation.history.str());
  }

  //if weight map provided, save it too
  if (weight.size() != 0)
    FITS::writeImage(outfptr,weight,"WEIGHT");

  FITS::closeFile(outfptr);
}
