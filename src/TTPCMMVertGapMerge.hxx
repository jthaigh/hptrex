#ifndef TTPCMMVertGapMerge_hxx_seen
#define TTPCMMVertGapMerge_hxx_seen

#include "TTPCSeed2CluMatchMerge.hxx"

namespace ND {
  class TTPCMMVertGapMerge;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class ND::TTPCMMVertGapMerge : private TTPCSeed2CluMatchMerge {
  public:
    /// Default constructor.
    TTPCMMVertGapMerge();
    /// Default destructor.
    ~TTPCMMVertGapMerge() {};

    /// Perform the matching and merging
    void Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns);

  private:
    bool HitSelected(ND::THandle<ND::TTPCHitPad> hitPad, int &Zone);
    void FindLocation(ND::THandle<ND::TTPCPattern> pattern, int &tpc, int &rp, int &updown);

    bool IsInDriftVolume(ND::THandle<ND::TTPCPattern> pattern, int tpc, int rp, int updown);

    bool fRun;
    int fNbWantedColumns;

};


#endif
