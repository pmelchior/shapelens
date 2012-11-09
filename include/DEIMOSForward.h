#ifndef SHAPELENS_DEIMOSFORWARD_H
#define SHAPELENS_DEIMOSFORWARD_H
#include "DEIMOSElliptical.h"
#include <vector>

namespace shapelens {
  typedef std::vector<Object> MultiExposureObject;
  typedef std::vector<Moments> MultiExposureMoments;
  class DEIMOSForward : public DEIMOS {
  public:
    DEIMOSForward(const MultiExposureObject& meo, const std::vector<DEIMOS::PSFMultiScale>& mePSFMultiScale, int N, int C);
    std::vector<DEIMOSElliptical> meD;
    /// Centroid of the object (in WCS coordinates).
    Point<data_t> centroid;
    /// S/N of moment measurement.
    data_t SN;
    data_t eta;
    /// History of modelling.
    History history;
    /// Modelling flags.
    /// - <tt>flags[0]</tt>: non-sensical moments occurred during one of the iterations
    /// - <tt>flags[1]</tt>: non-sensical moments occurred during the last iterations
    std::bitset<2> flags;

    /// Whether to fix the centroid of the object in all exposures.
    static bool FIX_CENTROID;
    static data_t ETA_MAX;
    
  protected:
    void computeMomentsFromGuess();
    data_t getWeightFunctionScale(unsigned int k) const;
    const MultiExposureObject& meo;
    const std::vector<DEIMOS::PSFMultiScale>& mePSFMultiScale;
    MultiExposureMoments mem;
    Moments mo0;
    std::vector<tmv::Matrix<data_t> > meP;
    void initialize(int C);
    void minimize();
    void convolveExposure(unsigned int k);
    int K;
    int t;
  };
} // end namespace
#endif
