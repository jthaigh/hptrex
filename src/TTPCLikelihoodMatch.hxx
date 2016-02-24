#ifndef TTPCLikelihoodMatch_hxx_seen
#define TTPCLikelihoodMatch_hxx_seen

#include <TReconBase.hxx>

#include "TTPCPattern.hxx"
#include "TTPCPath.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCLikFitPath.hxx"

namespace ND {
  class TTPCLikelihoodMatch;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class ND::TTPCLikelihoodMatch {
  public:
    TTPCLikelihoodMatch();
    ~TTPCLikelihoodMatch();

    void Process( ND::TReconObjectContainer *allPatterns);
    void MatchBrokenPaths(std::vector< ND::THandle<ND::TTPCPattern> > inPatterns);
    void MatchAcrossJunctions(ND::THandle<ND::TTPCPattern> Pattern);

  private:
    void MatchPathsAtJunction(ND::THandle<ND::TTPCPath> Path1, ND::THandle<ND::TTPCPath> Path2, int JunctionId);
    bool CheckClusterMatch(State &propState, ND::THandle<ND::TTPCHVCluster> Target, int targetSense, bool ForcePropagSense);
    bool CheckClusterMatch(State &propState, ND::THandle<ND::TTPCHVCluster> Target, int targetSense, bool ForcePropagSense, double &matchDistance);
    TTPCLogLikelihood GetMatchLikelihood(State trkState, ND::THandle<ND::THitSelection> cluToMatch, double trkLength);

    TTPCLikFitPath *fLklhdCalc;
};


#endif
