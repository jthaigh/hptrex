#ifndef TTPCLikelihoodMerge_hxx_seen
#define TTPCLikelihoodMerge_hxx_seen

#include <TReconBase.hxx>

#include "TTRExPattern.hxx"
#include "TTRExPath.hxx"
#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"
#include "TTPCSeeding.hxx"

#define MAXPATTERNCHAIN 10

namespace trex {
  class TTPCLikelihoodMerge;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class trex::TTPCLikelihoodMerge {
  public:
    /// Default constructor.
    TTPCLikelihoodMerge();
    /// Default destructor.
    virtual ~TTPCLikelihoodMerge() {};

    /// Perform the matching and merging
  void Process(std::vector<trex::TTRExPattern>& inputPatterns, std::vector<trex::TTRExPatterns>& mergedPatterns);

  private:

    void CleanUp();
    double PathToPatternMatch( trex::TTRExPath& PathA, trex::TTRExPattern& PatternB, trex::TTRExPath* bestPathB);
    void MatchBrokenPaths(std::vector< trex::TTRExPattern>& patterns);
        
    void MatchThroughJunctions(trex::TTRExPattern& pattern);
    trex::TTRExPattern* MergeAll();


    class MatchingTracker{
      private:
      trex::TTRExPath* fRawPath[2];
      trex::TTRExPath* fMergedPath;
      std::vector<trex::TTPCHitPad>* fJunction;
      bool fMergeMe;
      int fChain;

      public: 
        /// Default constructor.
      MatchingTracker(trex::TTRExPath& Path1, trex::TTRExPath& Path2, std::vector<trex::TTPCHitPad>& Junction){
	fRawPath[0] = &Path1;
	fRawPath[1] = &Path2;
	fJunction = &Junction;
	fMergeMe = true;
	fChain = -1;
          
      }
      MatchingTracker(trex::TTRExPath& Path1, trex::TTRExPath& Path2){
	fRawPath[0] = *Path1;
	fRawPath[1] = *Path2;
	fJunction = 0;
	fMergeMe = true;
	fChain = -1;
        
      }
      MatchingTracker( std::vector<trex::TTPCHitPad>& Junction){
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
      std::vector<trex::TTPCHitPad>* GetJunction() {
	return fJunction;
      };
    };

    std::vector<MatchingTracker> fMTracker;
    trex::TTRExPattern* fPatternChain[MAXPATTERNCHAIN];
    int fNbPatternChained;


    bool fRunThroughGoingMerging;
    bool fRunBrokenTracksMerging;
    double fThruGoDeltaLogLklhdCut;
    double fBrkTrkDeltaLogLklhdCut;

};


#endif

