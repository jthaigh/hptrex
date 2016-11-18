#ifndef TTPCLikelihoodMatch_hxx_seen
#define TTPCLikelihoodMatch_hxx_seen

#include "TTRExPattern.hxx"
#include "TTRExPath.hxx"
#include "TTRExHVCluster.hxx"
#include "TTPCLikFitPath.hxx"

namespace trex {
  class TTPCLikelihoodMatch;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class trex::TTPCLikelihoodMatch {
  public:
    TTPCLikelihoodMatch();
    ~TTPCLikelihoodMatch();

  void Process( std::vector<trex::TTRExPattern>& allPatterns);
    void MatchBrokenPaths(std::vector< trex::TTRExPattern>& inPatterns);
  void MatchAcrossJunctions(trex::TTRExPattern& Pattern);

  private:
  void MatchPathsAtJunction(trex::TTRExPath& Path1, trex::TTRExPath& Path2, int JunctionId);
  bool CheckClusterMatch(std::vector<double> propState, trex::TTRExHVCluster& Target, int targetSense, bool ForcePropagSense);
  bool CheckClusterMatch(std::vector<double> propState, trex::TTRExHVCluster& Target, int targetSense, bool ForcePropagSense, double &matchDistance);
  trex::TTPCLogLikelihood GetMatchLikelihood(std::vector<double> trkState, std::vector<trex::TTRExHVCluster*>& cluToMatch, double trkLength);

    TTPCLikFitPath *fLklhdCalc;
};


#endif
