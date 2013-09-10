#include <shapelens/ShapeLens.h>
#include <shapelens/PSFEx.h>
#include <tclap/CmdLine.h>

using namespace shapelens;

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Create PSF at specified location from PSFEx model", ' ', "0.1");
  TCLAP::ValueArg<std::string> psffile("p","psf_file","PSF moment file",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outfile("o","output","PSF image FITS file",true,"","string",cmd);
  TCLAP::ValueArg<int> patch("P","patch_size", "Size of simulated patch",false,0,"int",cmd);
  TCLAP::ValueArg<data_t> x("x","x_pos", "Position x of object centroid",true,0,"data_t",cmd);
  TCLAP::ValueArg<data_t> y("y","y_pos", "Position y of object centroid",true,0,"data_t",cmd);
  TCLAP::SwitchArg verbose("v","verbose","Verbose output", cmd, false);
  cmd.parse(argc,argv);

  Config::VERBOSITY = verbose.getValue();
  PSFEx p(psffile.getValue());

  Object psf;
  // size of PSF postage stamp:
  // this can be smaller than max(p.getSize(0), p.getSize(1))
  // since the PSF model only has finite power 
  int P = int(2*ceil(p.maxsize()));;
  if (patch.isSet())
    P = patch.getValue();
  psf.resize(P*P);
  psf.setZero();
  // build PSF model at position of galaxy
  psf.centroid(0) = x.getValue();
  psf.centroid(1) = y.getValue();
  psf.grid.setSize(int(floor(psf.centroid(0) - P/2)), int(floor(psf.centroid(1) - P/2)), P, P);
  p.fillObject(psf);
  psf.save(outfile.getValue());
}
