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

// TREx
#include <TTPCHitPad.hxx>
//#include <TTPCPattern.hxx>
//#include <TTPCPath.hxx>
//#include <TTPCJunction.hxx>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCAStar.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroupMan.hxx"
#include "TTPCOrderedVolGroup.hxx"
#include "TTPCVolGroup.hxx"
#include "TTRExHVCluster.hxx"
#include "TTRExPath.hxx"
#include "TTRExJunction.hxx"
#include "TTRExPattern.hxx"

namespace trex{
  /// Main algorithm for pattern recognition for path finding.  Processes a sub-event of connected hits and returns a set of paths.
  class TTPCTRExPatSubAlgorithm {
    public:
      /// Default constructor
      TTPCTRExPatSubAlgorithm(trex::TTPCLayout* layout);
      /// Default destructor
      ~TTPCTRExPatSubAlgorithm();

      TTPCTRExPatSubAlgorithm(const TTPCTRExPatSubAlgorithm& in) = delete;

    TTPCTRExPatSubAlgorithm(TTPCTRExPatSubAlgorithm&& in);
      

      /// Set up pattern recognition algorithm cells from an iterable map of previously defined ones
      void SetUpHits(std::map<long, trex::TTPCUnitVolume*>& map, trex::TTPCAStar* aStarCopy=0);

    //MDH
    //never used.
    /*
      /// Add a new cell to current map, either by creating a new one or incrementing the charge of an existing one depending on wheter position is already occupied
      void AbsorbCell(long id, trex::TTPCUnitVolume* cell);
      /// Add a list of new cells to current map, through a series of AbsorbCell calls on each element of the group
      void AppendHits(trex::THandle<trex::TTPCVolGroup> hits);
    */


      /// Fill containers in this sub-object
      void ProduceContainers();
      /// Ensure that none of the containers in this sub-object are empty
      void CleanContainers();
      
    /// Produce pattern so it can be returned
    void ProducePattern(trex::TTRExPattern& output);//std::vector<trex::TTPCHitPad*>& used);




    //NEED TO CHANGE RETURN TYPE HERE

    void ConnectJunctionAndPath(trex::TTRExJunction& junction, trex::TTRExPath& path);

      /// Get hits from internal group
    std::vector<trex::TTPCHitPad*> GetHits();

      /// Get hits from algorithm's hit map corresponding to provided path
    std::vector<trex::TTRExHVCluster*> GetHits(trex::TTPCOrderedVolGroup& path);

      /// Get groups of connected cells for defining sub-events
    void GetRegions(std::vector< trex::TTPCVolGroup >& regions);

    //MDH
    //Not currently used
      /// Fill used hits from this sub algorithm
    //      void FillUsedHits(trex::THitSelection* used);

      /// Get iterator for first element in vector of all paths
      std::vector< trex::TTPCOrderedVolGroup >::iterator GetPathsBegin(){ return fTracks.begin(); }
      /// Get iterator for last element in  vector of all paths
      std::vector< trex::TTPCOrderedVolGroup >::iterator GetPathsEnd(){ return fTracks.end(); }

      /// Get iterator for first element of hit map
      std::map<long, trex::TTPCUnitVolume*>::iterator GetHitMapBegin(){ return fHitMap.begin(); }
      /// Get iterator for last element of hit map
      std::map<long, trex::TTPCUnitVolume*>::iterator GetHitMapEnd(){ return fHitMap.end(); }
      /// Get size of hit map
      int GetHitMapSize(){ return fHitMap.size(); }

      /// Get whether or not this sub group has valid paths
      bool GetHasValidPaths(){ return fHasValidPaths; }
      /// Get reference to layout contained by this object
      trex::TTPCLayout* GetLayout(){ return fLayout; }
      /// Get map of cells contained by this object
      std::map<long, trex::TTPCUnitVolume*> GetHitMap(){ return fHitMap; }
      /// Get found paths contained by this object
      std::vector< trex::TTPCOrderedVolGroup >& GetTracks(){ return fTracks; }
      /// Get found paths contained by this object as unordered groups
      void GetTrackExtendedHits(std::vector<trex::TTPCVolGroup>& extHits);

      /// Get whether this is the primary sub group
      bool GetPrimary(){ return fPrimary; }

      /// Set whether this is the primary sub group
      void SetPrimary(bool primary=true){ fPrimary = primary; }

    private:
      /// Map of all cells and their unique IDs
      std::map<long, trex::TTPCUnitVolume*> fHitMap;
      /// Manager for handling all vol groups
      trex::TTPCVolGroupMan* fVolGroupMan;
      /// Path finder object for connecting points of interest
      trex::TTPCAStar* fAStar;

      /// Overall pattern produced
    //    trex::TTPCPattern* fPattern;

      /// Vector of all paths
      std::vector< trex::TTPCOrderedVolGroup > fTracks;

    //Vector of vertices (never actually read from, just used for persistency
    std::vector< trex::TTPCVolGroup > fVertices;

      /// Set maximum, minimum and range of cells in x, y and z
      void SetRanges(int rangeX,int minX,int maxX, int rangeY,int minY,int maxY, int rangeZ,int minZ,int maxZ);

      /// Whether the event contains any hits or not
      bool fHasHits;
      /// Whether or not this event contains valid paths
      bool fHasValidPaths;
      /// TPC number for processing algorithm (0 for all TPCs, 1, 2 or 3 to restrict to an individual one)
      int fTPC;
      /// Layout information for the current event
      trex::TTPCLayout* fLayout;

      /// Whether this is the primary sub group
      bool fPrimary;

  };
}

#endif
