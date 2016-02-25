#ifndef TTPCMMHoriGapMerge_hxx_seen
#define TTPCMMHoriGapMerge_hxx_seen

#include "TTPCSeed2CluMatchMerge.hxx"

namespace ND {
  class TTPCMMHoriGapMerge;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class ND::TTPCMMHoriGapMerge : private TTPCSeed2CluMatchMerge {
  public:
    /// Default constructor.
    TTPCMMHoriGapMerge();
    /// Default destructor.
    ~TTPCMMHoriGapMerge() {};

    /// Perform the matching and merging
    void Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns);

  private:
    bool IsClusterUsable(ND::THandle<ND::TTPCHVCluster> Cluster);
    bool HitSelected(ND::THandle<ND::TTPCHitPad> hitPad, int &Zone);
    void FindLocation(ND::THandle<ND::TTPCPattern> pattern, int &tpc, int &rp, std::map<int,int> &mm);

    bool IsInDriftVolume(ND::THandle<ND::TTPCPattern> pattern, int tpc, int rp, int updown);

    bool fRun;
    int fNbWantedRows;

};


#endif
