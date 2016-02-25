#ifndef TTPCTREXPATSUBALGORITHM_HXX
#define TTPCTREXPATSUBALGORITHM_HXX

// c++
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>

// ROOT
#include <TVector3.h>

// nd280
#include <TND280Log.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <THit.hxx>
#include <TReconPID.hxx>
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

// TREx
#include <TTPCHitPad.hxx>
#include <TTPCPattern.hxx>
#include <TTPCPath.hxx>
#include <TTPCJunction.hxx>
#include <TTPCDebug.hxx>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCAStar.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroupMan.hxx"
#include "TTPCOrderedVolGroup.hxx"
#include "TTPCVolGroup.hxx"
#include "TTPCCalibration.hxx"

namespace ND{
  /// Main algorithm for pattern recognition for path finding.  Processes a sub-event of connected hits and returns a set of paths.
  class TTPCTRExPatSubAlgorithm {
    public:
      /// Default constructor
      TTPCTRExPatSubAlgorithm(ND::TTPCLayout* layout);
      /// Default destructor
      ~TTPCTRExPatSubAlgorithm();

      /// Set up pattern recognition algorithm cells from an iterable map of previously defined ones
      void SetUpHits(std::map<long, ND::TTPCUnitVolume*> map, ND::TTPCAStar* aStarCopy=0);
      /// Set up high charge cut version of pattern recognition algorithm cells from an iterable map of previously defined ones
      void SetUpHitsHigh(std::map<long, ND::TTPCUnitVolume*> map, ND::TTPCAStar* aStarCopy=0);
      /// Add a new cell to current map, either by creating a new one or incrementing the charge of an existing one depending on wheter position is already occupied
      void AbsorbCell(long id, ND::TTPCUnitVolume* cell);
      /// Add a list of new cells to current map, through a series of AbsorbCell calls on each element of the group
      void AppendHits(ND::THandle<ND::TTPCVolGroup> hits);
      /// Fill containers in this sub-object
      void ProduceContainers();
      /// Endure that none of the containers in this sub-object are empty
      void CleanContainers();
      /// Produce pattern so it can be returned
      void ProducePattern(ND::THitSelection* used);

      /// Return and this object's pattern
      ND::THandle<ND::TTPCPattern> GetPattern();

      /// Get hits from internal group
      ND::THandle<ND::THitSelection> GetHits();
      /// Get hits from algorithm's hit map corresponding to provided path
      ND::THandle<ND::THitSelection> GetHits(ND::THandle<ND::TTPCOrderedVolGroup> path);

      /// Get groups of connected cells for defining sub-events
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetRegions();

      /// Fill used hits from this sub algorithm
      void FillUsedHits(ND::THitSelection* used);

      /// Get iterator for first element in vector of all paths
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator GetPathsBegin(){ return fTracks.begin(); }
      /// Get iterator for last element in  vector of all paths
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator GetPathsEnd(){ return fTracks.end(); }

      /// Get iterator for first element of hit map
      std::map<long, ND::TTPCUnitVolume*>::iterator GetHitMapBegin(){ return fHitMap.begin(); }
      /// Get iterator for last element of hit map
      std::map<long, ND::TTPCUnitVolume*>::iterator GetHitMapEnd(){ return fHitMap.end(); }
      /// Get size of hit map
      int GetHitMapSize(){ return fHitMap.size(); }

      /// Get drift speed for internal layout
      double GetDriftSpeed(){ return fLayout->GetDriftSpeed(); }
      /// Get x cathode crossing for internal layout
      bool GetXCathodeCross(){ return fLayout->GetXCathodeCross(); }
      /// Get whether or not this sub group has valid paths
      bool GetHasValidPaths(){ return fHasValidPaths; }
      /// Get reference to layout contained by this object
      ND::TTPCLayout* GetLayout(){ return fLayout; }
      /// Get map of cells contained by this object
      std::map<long, ND::TTPCUnitVolume*> GetHitMap(){ return fHitMap; }
      /// Get points of interest contained by this object
      std::vector<TVector3> GetPOIs(){ return fPOIs; }
      /// Get groups of interest contained by this object
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetGOIs(){ return fGOIs; }
      /// Get found paths contained by this object
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >GetTracks(){ return fTracks; }
      /// Get found paths contained by this object as unordered groups
      std::vector< ND::THandle<ND::TTPCVolGroup> >GetTrackExtendedHits();

      /// Get whether this is the primary sub group
      bool GetPrimary(){ return fPrimary; }

      /// Set whether this is the primary sub group
      void SetPrimary(bool primary=true){ fPrimary = primary; }

    private:
      /// Map of all cells and their unique IDs
      std::map<long, ND::TTPCUnitVolume*> fHitMap;
      /// Map of all cells passing a charge cut and their unique IDs
      std::map<long, ND::TTPCUnitVolume*> fHitMapHigh;
      /// Manager for handling all vol groups
      ND::TTPCVolGroupMan* fVolGroupMan;
      /// Manager for handling all vol groups that pass a charge threshold cut
      ND::TTPCVolGroupMan* fVolGroupManHigh;
      /// Path finder object for connecting points of interest
      ND::TTPCAStar* fAStar;
      /// Path finder object for connecting points of interest through hits passing a charge threshold cut
      ND::TTPCAStar* fAStarHigh;

      /// Overall pattern produced
      ND::THandle<ND::TTPCPattern> fPattern;

      /// Vector of all paths
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > fTracks;
      /// Vector of all points of interest
      std::vector<TVector3> fPOIs;
      /// Vector of all groups of interest
      std::vector< ND::THandle<ND::TTPCVolGroup> > fGOIs;

      /// Set maximum, minimum and range of cells in x, y and z
      void SetRanges(int rangeX,int minX,int maxX, int rangeY,int minY,int maxY, int rangeZ,int minZ,int maxZ);

      /// Whether the event contains any hits or not
      bool fHasHits;
      /// Whether or not this event contains valid paths
      bool fHasValidPaths;
      /// TPC number for processing algorithm (0 for all TPCs, 1, 2 or 3 to restrict to an individual one)
      int fTPC;
      /// Layout information for the current event
      ND::TTPCLayout* fLayout;

      /// Whether this is the primary sub group
      bool fPrimary;
  };
}

#endif
