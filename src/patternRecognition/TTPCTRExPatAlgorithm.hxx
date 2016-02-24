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

// nd280
#include <TND280Log.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <THit.hxx>
#include <TComboHit.hxx>
#include <TReconBase.hxx>
#include <TReconTrack.hxx>
#include <TReconHit.hxx>
#include "TTPCCalibration.hxx"
#include <createTReconTrack.hxx>
#include <TChannelId.hxx>
#include <TTPCChannelId.hxx>

// TREx
#include <TTPCPattern.hxx>
#include <TTPCPath.hxx>
#include <TTPCJunction.hxx>
#include <TTPCDebug.hxx>

// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"
#include "TTPCLayout.hxx"
#include "TTPCFeatureFinder.hxx"
#include "TTPCAStar.hxx"
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

      /// Run the pattern recognition and return TTPCPattern containing result
      void GetPatterns(TReconObjectContainer *foundPatterns);

      /// Get drift speed for internal main algorithm
      double GetDriftSpeed(){ return fMasterLayout->GetDriftSpeed(); }
      /// Get x cathode crossing for internal main algorithm
      bool GetXCathodeCross(){ return fMasterLayout->GetXCathodeCross(); }

      /// Current processing pattern recognition
      void Process(const ND::TAlgorithmResult &in, ND::THitSelection* used, ND::THitSelection* unused);
      /// Current processing pattern recognition
      void Process(ND::THandle<ND::THitSelection> hits, ND::THitSelection* used, ND::THitSelection* unused);

      /// Getters
      /// Get iterator to start of set of sub-algorithms for sub-events in event
      std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator GetSubAlgorithmsBegin(){ return fSubAlgorithms.begin(); }
      /// Get iterator to end of set of sub-algorithms for sub-events in event
      std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator GetSubAlgorithmsEnd(){ return fSubAlgorithms.end(); }
      /// Get reference to layout contained by this object
      ND::TTPCLayout* GetLayout(){ return fMasterLayout; }
      /// Get reference to feature finder contained by this object
      ND::TTPCFeatureFinder* GetFeatureFinder(){ return fFeatureFinder; }
      /// Get map of cells contained by this object
      std::map<long, ND::TTPCUnitVolume*> GetHitMap(){ return fMasterHitMap; }
      /// Get object containing all sub-events
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetSubEvents(){ return fSubEvents; }
      /// Get object containing all attached groups of delta hits
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetDeltaHits(){ return fDeltaHits; }

    private:
      /// Master layout to use for this event
      ND::TTPCLayout* fMasterLayout;
      /// Manager for handling all vol groups
      ND::TTPCVolGroupMan* fMasterVolGroupMan;
      /// Manager for handling vol groups passing a high charge cut
      ND::TTPCVolGroupMan* fMasterVolGroupManHigh;
      /// Master A* algorithm for making paths
      ND::TTPCAStar* fMasterAStar;
      /// Feature finder for identifying broad features
      ND::TTPCFeatureFinder* fFeatureFinder;

      /// Master map of all hits in this event
      std::map<long, ND::TTPCUnitVolume*> fMasterHitMap;
      /// Master map of all hits in this event passing a charge cut
      std::map<long, ND::TTPCUnitVolume*> fMasterHitMapHigh;
      /// Object containing all hit pads
      ND::THandle<ND::THitSelection> fHits;
      /// Object containing all sub-events
      std::vector< ND::THandle<ND::TTPCVolGroup> > fSubEvents;
      /// Object containing all attached groups of delta hits
      std::vector< ND::THandle<ND::TTPCVolGroup> > fDeltaHits;

      /// Whether hits have been added
      bool fHasHits;

      /// Drift velocity
      double fDriftVelocity;

      /// Set of sub-algorithms for sub-events in event
      std::vector<ND::TTPCTRExPatSubAlgorithm*> fSubAlgorithms;

      /// Add a selection of hits for the first time and work out preliminary t0 and cathode crossing, and set up hits
      void PrepareHits(ND::THandle<ND::THitSelection> hits);
      /// Populate list of delta ray hits
      void PopulateDeltaHits();
      /// Fill hits with information about whether they were found to be low charge
      void FillInfoLowCharge();
      /// Fill map with information about whether they were associated with an ASIC with too many hits
      void FillInfoFullASIC();
      /// Fill map with information about whether they were found associated with an ASIC with saturated hits
      void FillInfoSatASIC();
      /// Fill map with information about whether they were found associated with early or late negative peaks
      void FillInfoNegativePeak();
      /// Fill list of high quality map with those found not to be pathological
      void FillHighQualityHits();

      /// Fill unused hits given input of used hits
      void FillUsedUnusedHits(ND::THitSelection* usedTREx, ND::THitSelection* used, ND::THitSelection* unused);

      /// Print debug messages for hits based on desired level of verbosity, at pattern recognition level
      void VerifyHitsPattern();
      /// Print debug messages on the status of unused and unused hits, at pattern recognition level
      void VerifyUsedUnusedPattern(ND::THitSelection* used, ND::THitSelection* unused);
      /// Print debug messages for hits based on desired level of verbosity, at final result level
      void VerifyHitsFull(ND::TReconObjectContainer* patterns);
      /// Print debug messages on the status of unused and unused hits, at final result level
      void VerifyUsedUnusedFull(ND::TReconObjectContainer* patterns);
      /// Print debug messages for paths and junctions in a pattern
      void VerifyPatternConstituents(ND::THandle<ND::TTPCPattern> pattern);
  };
}

#endif
