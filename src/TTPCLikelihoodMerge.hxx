#ifndef TTPCLikelihoodMerge_hxx_seen
#define TTPCLikelihoodMerge_hxx_seen

#include <TReconBase.hxx>

#include "TTPCPattern.hxx"
#include "TTPCPath.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCSeeding.hxx"

#define MAXPATTERNCHAIN 10

namespace ND {
  class TTPCLikelihoodMerge;
}

/// Algorithm to match and merge the patterns on each side of the vert MM gap
class ND::TTPCLikelihoodMerge {
  public:
    /// Default constructor.
    TTPCLikelihoodMerge();
    /// Default destructor.
    virtual ~TTPCLikelihoodMerge() {};

    /// Perform the matching and merging
    void Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns);

  private:

    void CleanUp();
    double PathToPatternMatch( ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPattern> PatternB, ND::THandle<ND::TTPCPath> &bestPathB);
    void MatchBrokenPaths(std::vector< ND::THandle<ND::TTPCPattern> > patterns);
        
    void MatchThroughJunctions(ND::THandle<ND::TTPCPattern> pattern);
    ND::THandle<ND::TTPCPattern> MergeAll();


    class MatchingTracker{
      private:
        ND::THandle<ND::TTPCPath> fRawPath[2];
        ND::THandle<ND::TTPCPath> fMergedPath;
        ND::THandle<ND::TTPCJunction> fJunction;
        bool fMergeMe;
        int fChain;

      public: 
        /// Default constructor.
        MatchingTracker(ND::THandle<ND::TTPCPath> Path1, ND::THandle<ND::TTPCPath> Path2, ND::THandle<ND::TTPCJunction> Junction){
          fRawPath[0] = Path1;
          fRawPath[1] = Path2;
          fJunction = Junction;
          fMergeMe = true;
          fChain = -1;
          
        }
        MatchingTracker(ND::THandle<ND::TTPCPath> Path1, ND::THandle<ND::TTPCPath> Path2){
          fRawPath[0] = Path1;
          fRawPath[1] = Path2;
          fJunction = ND::THandle<ND::TTPCJunction> ();
          fMergeMe = true;
          fChain = -1;
          
        }
        MatchingTracker( ND::THandle<ND::TTPCJunction> Junction){
          fRawPath[0] = ND::THandle<ND::TTPCPath> ();
          fRawPath[1] = ND::THandle<ND::TTPCPath> ();
          fMergedPath = ND::THandle<ND::TTPCPath> ();
          fJunction = Junction;
          fMergeMe = false;
          fChain = -1;
        }
        /// Default destructor.
        virtual ~MatchingTracker() {};
        bool NeedsMerging(){return (fMergeMe && fChain<0) ;};
        bool HasThisPath(ND::THandle<ND::TTPCPath> testPath){ return (fRawPath[0] == testPath || fRawPath[1] == testPath);};
        ND::THandle<ND::TTPCPath> GetRawPath(int i) {return fRawPath[i];};
        ND::THandle<ND::TTPCPath> GetMergedPath() {return fMergedPath;};
        void SetMergedPath(ND::THandle<ND::TTPCPath> merged) {fMergedPath = merged;};
        bool IsMergingChainOk(int ChainId) {return (fChain == ChainId);};
        void SetMergingChain(int ChainId) {fChain = ChainId; fMergeMe = false;};
        ND::THandle<ND::TTPCJunction> GetJunction() {
          if (!fJunction) return ND::THandle<ND::TTPCJunction> ();
          return fJunction;
        };
    };

    std::vector<MatchingTracker> fMTracker;
    ND::THandle<ND::TTPCPattern> fPatternChain[MAXPATTERNCHAIN];
    int fNbPatternChained;


    bool fRunThroughGoingMerging;
    bool fRunBrokenTracksMerging;
    double fThruGoDeltaLogLklhdCut;
    double fBrkTrkDeltaLogLklhdCut;

};


#endif

