#ifndef TTPCCathCrosserMerge_hxx_seen
#define TTPCCathCrosserMerge_hxx_seen

#include <TReconBase.hxx>

#include "TTPCPattern.hxx"
#include "TTPCPath.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCT0Finder.hxx"

#define MAXNBCATHMATCHCANDIDATE 100
#define MAXNBCATHPROPAGCANDIDATE 100

namespace ND {
  class TTPCCathCrosserMerge;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class ND::TTPCCathCrosserMerge {
  public:
    /// Default constructor.
    TTPCCathCrosserMerge(ND::TTPCT0Finder *T0Finder);
    /// Default destructor.
    ~TTPCCathCrosserMerge() {};

    /// Perform the matching and merging
    void Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns);

  private:
//    bool HitSelected(ND::THandle<ND::TTPCHitPad> hitPad, int &Zone);
    void FindLocation(ND::THandle<ND::TTPCPattern> pattern, int &tpc, int &rp);

    void FindPropagCandidates(ND::THandle<ND::TTPCPattern> pattern, int tpc, int matchRP);
//    bool IsInDriftVolume(ND::THandle<ND::TTPCPattern> pattern, int tpc, int rp, int updown);
    void CleanUp();
    void CleanUpMatchCand();
    void CleanUpPropagCand();
    void MatchPatternsAWithB(std::vector< ND::THandle<ND::TReconBase> > &rawPat, int zoneA, int zoneB, std::vector< ND::THandle<ND::TReconBase> > &resPat);
    bool FindMatchCandidates(ND::THandle<ND::TTPCPattern> pattern, int zone);
    void FindMatches(int matchId, ND::THandle<ND::TTPCPattern> MatchCandPat);
    ND::THandle<ND::TTPCPattern> MergePatterns(ND::THandle<ND::TTPCPattern> origPattern, int matchId);
    void MigratePatterns(ND::THandle<ND::TTPCPattern> startPattern, unsigned int oldPathId, ND::THandle<ND::TTPCPath> newPath, ND::THandle<ND::TTPCPattern> endPattern);

    bool fRun;
    ND::TTPCT0Finder *fT0Finder;

    /// Defines the clusters or hit pads that will be matched to a path.
    struct MatchCandidate {
      public:
        ND::THandle<ND::TTPCPath> Path;
        ND::THandle<ND::THitSelection> ClustersToMatch;
        int MatchedPathId;
        double ReduChi2;
    };

    /// Defines the path's state that will be propagated and match to clusters
    /// of another path or hit pads in a junction.
    struct PropagCandidate {
      public:
        ND::THandle<ND::TTPCPattern> Pattern;
        ND::THandle<ND::TTPCPath> Path;
        State propagState[2];
    };

    // Have two groups: Upstream and downstream of the MM gap.
    MatchCandidate fMatchCand[MAXNBCATHMATCHCANDIDATE];
    int fNbMatchCand;

    PropagCandidate fPropagCand[MAXNBCATHPROPAGCANDIDATE];
    int fNbPropagCand;

    int fMaxNbMatchCand;
    int fMaxNbPropagCand;

    double fMaxReduChi2Match;
    double fMaxFirstResidual;

    Surface fCathSurf[3][2];
    double fRangeInZ[3][2];
    double fRangeInY[2];
};


#endif
