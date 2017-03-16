#ifndef TTPCLinearMerge_hxx_seen
#define TTPCLinearMerge_hxx_seen

#include "TTRExPattern.hxx"
#include "TTRExPath.hxx"
#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"
#include "TTPCSeeding.hxx"

namespace trex {
  class TTPCLinearMerge;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class trex::TTPCLinearMerge {
  public:
    /// Default constructor.
    TTPCLinearMerge();
    /// Default destructor.
    virtual ~TTPCLinearMerge() {};

    /// Perform the matching and merging
  void Process(std::vector<trex::TTRExPattern>& inputPatterns, std::vector<trex::TTRExPattern>& mergedPatterns);

  private:

    void CleanUp();

  double CalcChi2PerDOF( trex::TTRExPath& PathA, trex::TTRExPath& PathB);

  //    double PathToPatternMatch( trex::TTRExPath& PathA, trex::TTRExPattern& PatternB, trex::TTRExPath* bestPathB);
  //  void MatchBrokenPaths(std::vector< trex::TTRExPattern>& patterns);
        
    void MatchThroughJunctions(trex::TTRExPattern& pattern);

    void MergeAll(trex::TTRExPattern& inputPattern, std::vector<trex::TTRExPattern>& outputVector);

    class MatchingTracker{
      private:
      trex::TTRExPath* fRawPath[2];
      trex::TTRExPath* fMergedPath;
      trex::TTRExJunction* fJunction;
      bool fMergeMe;
      int fChain;

      public: 
        /// Default constructor.
      MatchingTracker(trex::TTRExPath& Path1, trex::TTRExPath& Path2, trex::TTRExJunction& Junction){
	fRawPath[0] = &Path1;
	fRawPath[1] = &Path2;
	fJunction = &Junction;
	fMergeMe = true;
	fChain = -1;
          
      }
      MatchingTracker(trex::TTRExPath& Path1, trex::TTRExPath& Path2){
	fRawPath[0] = &Path1;
	fRawPath[1] = &Path2;
	fJunction = 0;
	fMergeMe = true;
	fChain = -1;
        
      }
      MatchingTracker( trex::TTRExJunction& Junction){
	fRawPath[0] = 0;
	fRawPath[1] = 0;
	fMergedPath = 0;
	fJunction = &Junction;
	fMergeMe = false;
	fChain = -1;
      }
      /// Default destructor.
      virtual ~MatchingTracker() {};
      bool NeedsMerging(){return (fMergeMe && fChain<0) ;};
      bool HasThisPath(trex::TTRExPath& testPath){ return (fRawPath[0] == &testPath || fRawPath[1] == &testPath);};
      trex::TTRExPath* GetRawPath(int i) {return fRawPath[i];};
      trex::TTRExPath* GetMergedPath() {return fMergedPath;};
      void SetMergedPath(trex::TTRExPath& merged) {fMergedPath = &merged;};
      bool IsMergingChainOk(int ChainId) {return (fChain == ChainId);};
      void SetMergingChain(int ChainId) {fChain = ChainId; fMergeMe = false;};
      trex::TTRExJunction* GetJunction() {
	return fJunction;
      };
    };

    std::vector<MatchingTracker> fMTracker;

    bool fRunThroughGoingMerging;
  //bool fRunBrokenTracksMerging;
    double fThruGoDeltaLogLklhdCut;
  //double fBrkTrkDeltaLogLklhdCut;

};


#endif

