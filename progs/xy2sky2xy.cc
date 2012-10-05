#include <shapelens/FITS.h>
#include <shapelens/WCSTransformation.h>
#include <tclap/CmdLine.h>
#include <iomanip>

using namespace shapelens;

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Transform coordinates from pixel to sky (and back)", ' ', "0.1");
  TCLAP::SwitchArg reverse("r","reverse","Use reverse direction (sky to pixel)", cmd, false);
  TCLAP::SwitchArg intermediate("i","intermediate","Use intermediate transformation on tangent plane", cmd, false);
  TCLAP::ValueArg<std::string> filename("f","file", "FITS file with WCS header information", true, "", "string", cmd);
  TCLAP::UnlabeledValueArg<data_t> p1("p1","First coordinate of point",true,0,"data_t", cmd);
  TCLAP::UnlabeledValueArg<data_t> p2("p2","Second coordinate of point",true,0,"data_t", cmd);
  cmd.parse(argc,argv);

  Point<data_t> P(p1.getValue(), p2.getValue());
  fitsfile *fptr = FITS::openFile(filename.getValue());
  WCSTransformation wcs(fptr, intermediate.getValue());
  std::cout.flags(std::ios_base::fixed);
  if (reverse.isSet()) {
    wcs.inverse_transform(P);
    // shapelens uses (0/0) for left-lower corner, WCS (1/1)
    P(0) += 1;
    P(1) += 1;
    std::cout << std::setprecision(2) << P(0) << "\t" << P(1) << std::endl;
  } else {
    // shapelens uses (0/0) for left-lower corner, WCS (1/1)
    P(0) -= 1;
    P(1) -= 1;
    wcs.transform(P);
    std::cout << std::setprecision(6) << P(0) << "\t" << P(1) << std::endl;
  }
}
