#ifndef SHAPELENS_DEIMOSFORWARD_H
#define SHAPELENS_DEIMOSFORWARD_H
#include "DEIMOSElliptical.h"
#include <vector>

namespace shapelens {
  typedef std::vector<Object> MultiExposureObject;
  typedef std::vector<Moments> MultiExposureMoments;
  class DEIMOSForward : public DEIMOS {
  public:
    DEIMOSForward(const MultiExposureObject& meo, const MultiExposureObject& mepsf, int N, int C, const std::set<data_t>& scales);
    //DEIMOSForward(const MultiExposureObject& meo, const std::vector<DEIMOS::PSFMultiScale>& mePSFMultiScale, int N, int C);
    std::vector<DEIMOSElliptical> meD;
    /// Centroid of the object (in WCS coordinates).
    Point<data_t> centroid;
    /// S/N of moment measurement.
    std::map<data_t, data_t> SN;
    std::map<data_t, data_t> eta;
    data_t matching_scale;
    /// History of modelling.
    History history;
    /// Modelling flags.
    /// - <tt>flags[0]</tt>: non-sensical moments occurred during one of the iterations
    /// - <tt>flags[1]</tt>: non-sensical moments occurred during the last iterations
    std::bitset<2> flags;

    /// Get PSF moments for each exposure
    const MultiExposureMoments& getPSFMoments() const;

    /// Whether to fix the centroid of the object in all exposures.
    static bool FIX_CENTROID;
    static data_t ETA_MAX;
    
  protected:
    void computeMomentsFromGuess();
    const MultiExposureObject& meo, mepsf;
    //std::vector<DEIMOS::PSFMultiScale> mePSFMultiScale;
    MultiExposureMoments mem, mem_psf;
    DEIMOS::PSFMultiScale memom;
    std::vector<tmv::Matrix<data_t> > meP;
    void initialize(const std::set<data_t>& scales, int C);
    void minimize(const std::set<data_t>& scales);
    void minimizeAtFixedScale();
    data_t computeEta() const;
    void convolveExposure(unsigned int k);
    int K;
    int t;
    void rollBack(data_t scale, const Moments& best_mo, const Point<data_t>& best_centroid, const std::vector<DEIMOSElliptical>& best_meD, const tmv::Matrix<data_t>& best_S, const MultiExposureMoments& best_mem_psf);
  };
} // end namespace
#endif
