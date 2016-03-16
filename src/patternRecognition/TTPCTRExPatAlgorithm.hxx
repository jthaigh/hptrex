#ifndef TTPCTREXPATALGORITHM_HXX
#define TTPCTREXPATALGORITHM_HXX

// c++
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

// ROOT
#include <TVector3.h>

// TREx
#include <TTPCPattern.hxx>
#include <TTPCPath.hxx>
#include <TTPCJunction.hxx>
#include <TTPCDebug.hxx>

// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"
#include "TTPCLayout.hxx"
#include "TTPCFeatureFinder.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroupMan.hxx"
#include "TTPCOrderedVolGroup.hxx"
#include "TTPCVolGroup.hxx"

namespace ND{
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

      /// Get drift speed for internal main algorithm
      double GetDriftSpeed(){ return fMasterLayout->GetDriftSpeed(); }
      /// Get x cathode crossing for internal main algorithm
      bool GetXCathodeCross(){ return fMasterLayout->GetXCathodeCross(); }

      /// Current processing pattern recognition
    void Process(const std::vector<ND::TTPCHitPad>& hits, std::vector<ND::TTPCHitPad*>& used, std::vector<ND::TTPCHitPad*>& unused);

      /// Getters
      /// Get iterator to start of set of sub-algorithms for sub-events in event
      std::vector<ND::TTPCTRExPatSubAlgorithm>::iterator GetSubAlgorithmsBegin(){ return fSubAlgorithms.begin(); }
      /// Get iterator to end of set of sub-algorithms for sub-events in event
      std::vector<ND::TTPCTRExPatSubAlgorithm>::iterator GetSubAlgorithmsEnd(){ return fSubAlgorithms.end(); }
      /// Get reference to layout contained by this object
      ND::TTPCLayout* GetLayout(){ return fMasterLayout; }
      /// Get reference to feature finder contained by this object
      std::map<long, ND::TTPCUnitVolume>& GetHitMap(){ return fMasterHitMap; }
      /// Get object containing all sub-events
      std::vector< ND::TTPCVolGroup >& GetSubEvents(){ return fSubEvents; }
      /// Get object containing all attached groups of delta hits
      std::vector< ND::TTPCVolGroup >& GetDeltaHits(){ return fDeltaHits; }

    private:

    /// Add a selection of hits for the first time and work out preliminary t0 and cathode crossing, and set up hits
    void PrepareHits(std::vector<ND::TTPCHitPad> hits);
    /// Populate list of delta ray hits
    void PopulateDeltaHits();

      /// Master layout to use for this event
      ND::TTPCLayout* fMasterLayout;
      /// Manager for handling all vol groups
      ND::TTPCVolGroupMan* fMasterVolGroupMan;

      /// Master map of all hits in this event
      //MDH
      //This class owns the hits in the map and will delete them 
      //on destruction.
      std::map<long, ND::TTPCUnitVolume*> fMasterHitMap;
      /// Object containing all hit pads
      std::vector<const ND::TTPCHitPad*> fHits;
      /// Object containing all sub-events
      std::vector< ND::TTPCVolGroup > fSubEvents;
      /// Object containing all attached groups of delta hits
      std::vector< ND::TTPCVolGroup > fDeltaHits;

      /// Whether hits have been added
      bool fHasHits;

      /// Drift velocity
      double fDriftVelocity;

      /// Set of sub-algorithms for sub-events in event
      std::vector<ND::TTPCTRExPatSubAlgorithm> fSubAlgorithms;

    
  };
}

#endif
