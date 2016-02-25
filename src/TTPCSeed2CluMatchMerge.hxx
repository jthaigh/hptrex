#ifndef TTPCSeed2CluMatchMerge_hxx_seen
#define TTPCSeed2CluMatchMerge_hxx_seen

#include <TReconBase.hxx>

#include "TTPCPattern.hxx"
#include "TTPCPath.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCSeeding.hxx"

#define MAXNBMATCHCANDIDATE 100
#define MAXNBPROPAGCANDIDATE 100

namespace ND {
  class TTPCSeed2CluMatchMerge;
}

enum TTPCSeed2CluMatchMergeAlgo {kNoSeed2Clu = 0, kMMHoriGapMerge, kMMVertGapMerge};

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class ND::TTPCSeed2CluMatchMerge {
  public:
    /// Default constructor.
    TTPCSeed2CluMatchMerge();
    /// Default destructor.
    virtual ~TTPCSeed2CluMatchMerge() {};

    /// Perform the matching and merging
    virtual void Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns) {};

  protected:

    /// Overload this if you need to select special type of clusters.
    /// For example for the MMHoriGapMerge, we don't want clusters with hits
    /// from two different MM.
    virtual bool IsClusterUsable(ND::THandle<ND::TTPCHVCluster> Cluster) {return true;};

    /// MUST BE OVERLOADED
    virtual bool HitSelected(ND::THandle<ND::TTPCHitPad> hitPad, int &Zone) {return false;};

    /// Should probably remain the same
    void CleanUp();
    void CleanUpMatchCand();
    void CleanUpPropagCand();
    void MatchPatternsAWithB(std::vector< ND::THandle<ND::TReconBase> > &rawPat, int zoneA, int zoneB, std::vector< ND::THandle<ND::TReconBase> > &resPat);
    void GetPathEndsAtZoneEdge(ND::THandle<ND::TTPCPath> path, int &FirstLastAll, int &zone);
    void FindPropagCandidates(ND::THandle<ND::TTPCPattern> pattern, int zone);
    bool FindMatchCandidates(ND::THandle<ND::TTPCPattern> pattern, int zone);
    void FindMatches(int matchId, ND::THandle<ND::TTPCPattern> MatchCandPat);

    ND::THandle<ND::TTPCPattern> MergePatterns(ND::THandle<ND::TTPCPattern> origPattern, int matchId);
    void MigrateJunctions(ND::THandle<ND::TTPCPattern> startPattern, unsigned int oldPathId, ND::THandle<ND::TTPCPath> newPath, ND::THandle<ND::TTPCPattern> endPattern);

    /// Define the minimum distance to have a match between 2 paths
    double fMinDistForMatchP2P;
    /// Define the minimum distance to have a match between a path and a junction
    double fMinDistForMatchP2J;

    ND::TReconBase::StateBits fTPCBits[3];

    /// Redo the seeding after tracks were merged
    ND::TTPCSeeding* fNewSeeding;

    /// Just to know which algorithm is using the methods.
    TTPCSeed2CluMatchMergeAlgo fAlgo;

    /// Flag to create fake junctions instead of merging two paths.
    bool fFakeJunctionForStudies;

    /// Defines the clusters or hit pads that will be matched to a path.
    struct MatchCandidate {
      public:
        ND::THandle<ND::TReconBase> PathOrJunction;
        ND::THandle<ND::TTPCHVCluster> ClusterToMatch[2];
        unsigned int NbClusterToMatch;
        TVector3 JuncPosToMatch;
        int MatchedPathId;
        TVector3 Residuals;
        int MatchedClusterIdx;
    };

    /// Defines the path's state that will be propagated and match to clusters
    /// of another path or hit pads in a junction.
    struct PropagCandidate {
      public:
        ND::THandle<ND::TTPCPattern> Pattern;
        ND::THandle<ND::TTPCPath> Path;
        State propagState;
        int FrontOrBack;
    };

    // Have two groups: Upstream and downstream of the MM gap.
    MatchCandidate fMatchCand[MAXNBMATCHCANDIDATE];
    int fNbMatchCand;

    PropagCandidate fPropagCand[MAXNBPROPAGCANDIDATE];
    int fNbPropagCand;

    int fMaxNbMatchCand;
    int fMaxNbPropagCand;
};


#endif
