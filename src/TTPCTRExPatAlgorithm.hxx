#ifndef TTPCTREXPATALGORITHM_HXX
#define TTPCTREXPATALGORITHM_HXX

// c++
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

// ROOT
#include <TFile.h>
#include <TStyle.h>

// TREx
#include "TTrueHit.hxx"
#include "TTRExPattern.hxx"

// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"
#include "TTPCLayout.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroupMan.hxx"
#include "TTPCOrderedVolGroup.hxx"
#include "TTPCVolGroup.hxx"
#include "TTRExHVCluster.hxx"
#include "TTRExPath.hxx"
#include "TTRExJunction.hxx"

namespace trex{
  /// Pattern recognition for TPC tracks based on path finding.
  class TTPCTRExPatAlgorithm {
    public:
      /// Default constructor
    TTPCTRExPatAlgorithm();
      /// Default destructor
      ~TTPCTRExPatAlgorithm();

      /// Clean up values from previous processing
      void CleanUp();

    //Output - reimplement
      /// Run the pattern recognition and return TTPCPattern containing result
      //void GetPatterns(TReconObjectContainer *foundPatterns);

      /// Current processing pattern recognition
    void Process(std::vector<trex::TTPCHitPad*>& hits, std::vector<trex::TTPCHitPad*>& used, std::vector<trex::TTPCHitPad*>& unused, std::vector<TTrueHit*>& trueHits, trex::TTRExEvent* event, trex::TTPCLayout& layout);

      /// Getters
      /// Get iterator to start of set of sub-algorithms for sub-events in event
      std::vector<trex::TTPCTRExPatSubAlgorithm>::iterator GetSubAlgorithmsBegin(){ return fSubAlgorithms.begin(); }
      /// Get iterator to end of set of sub-algorithms for sub-events in event
      std::vector<trex::TTPCTRExPatSubAlgorithm>::iterator GetSubAlgorithmsEnd(){ return fSubAlgorithms.end(); }
      /// Get reference to layout contained by this object
      trex::TTPCLayout* GetLayout(){ return fMasterLayout; }
      /// Get reference to feature finder contained by this object
      std::map<long, trex::TTPCUnitVolume*>& GetHitMap(){ return fMasterHitMap; }

    private:

    /// Add a selection of hits for the first time and work out preliminary t0 and cathode crossing, and set up hits
    void PrepareHits(std::vector<trex::TTPCHitPad*>& hits);

      /// Master layout to use for this event
      trex::TTPCLayout* fMasterLayout;
      /// Manager for handling all vol groups
      trex::TTPCVolGroupMan* fMasterVolGroupMan;

      /// Master map of all hits in this event
      //MDH
      //This class owns the hits in the map and will delete them 
      //on destruction.
      std::map<long, trex::TTPCUnitVolume*> fMasterHitMap;
      /// Object containing all hit pads
      std::vector<trex::TTPCHitPad*> fHits;
      /// Object containing all sub-events

      /// Whether hits have been added
      bool fHasHits;

      /// Set of sub-algorithms for sub-events in event
      std::vector<trex::TTPCTRExPatSubAlgorithm> fSubAlgorithms;

    TFile* fPlotFile;
    
  };
}

#endif
