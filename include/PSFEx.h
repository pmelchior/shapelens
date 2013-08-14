#include "ShapeLens.h"
#include "Object.h"
#include <vector>
#include <string>

namespace shapelens {
  /// Wrapper class for PSFEx.
  /// Inspired by code from Daniel Gruen (LSW Munich), this allows to sample a 
  /// given PSFEx model on the pixel grid.
  class PSFEx {
  public:
    /// Argumented constructor.
    /// \p filename needs to point to a valid PSFEx \p psf file.
    PSFEx(std::string filename);
    /// The maximum size (in pixel) where the model has power.
    data_t maxsize();
    /// Fill the Object container with a pixellated PSFEx model.
    /// The extent of the model are taken from Object::grid, the nominal
    /// center of the model is set by Object::centroid.\n
    /// \b CAUTION: PSFEx models do not guarantee that the centroid of the
    /// pixellated model coincides with the given centroid position.
    void fillObject(Object& psf);
  private:
    double sinc(double x);
    static const data_t INTERPFAC = 3.0;
    std::vector<data_t> polzero, polscale;
    std::vector<int> psfaxis;
    int poldeg;
    data_t psf_samp;
    std::vector<Image<data_t> > basis;
  };
}
