#ifndef TTPCASTAR_HXX
#define TTPCASTAR_HXX

// c++
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <stdexcept>

// nd280
#include <TND280Log.hxx>
#include <THandle.hxx>

// ROOT
#include <TVector3.h>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCPathVolume.hxx"
#include "TTPCVolGroupMan.hxx"
#include "TTPCOrderedVolGroup.hxx"
#include "TTPCVolGroup.hxx"

namespace ND{
  /// Struct for holding elements for A* algorithm
  struct TTPCAStarPoint{
    ND::TTPCUnitVolume* vol;

    int x;
    int y;
    int z;

    float aStarCost;
    float aStarHeuristic;

    bool aStarHasFriends;
    bool aStarClosed;
    bool aStarOpen;

    ND::TTPCAStarPoint* aStarParent;
    std::map<ND::TTPCAStarPoint*, float> aStarFriends;
  };
  /// Performs an AStar algorithm to connect two points in a selection via a chain of valid cells
  class TTPCAStar{
    public:
      /// Default constructor
      TTPCAStar(ND::TTPCLayout* layout);
      /// Default destructor
      virtual ~TTPCAStar();

      /// Add the hits to use for the path finding
      void AddHits(ND::TTPCVolGroupMan* volGroupMan, std::map<long, ND::TTPCUnitVolume*> hitMap, bool extendMode=false);
      /// Add hits, saving time by copying from a previous A* container
      void AddHits(ND::TTPCAStar* prevAStar, std::map<long, ND::TTPCUnitVolume*> hitMap);

      /// Connect a group of vertices to a group of track ends
      std::vector< ND::THandle<ND::TTPCVolGroup> > ConnectVertexGroups(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices, std::vector< ND::THandle<ND::TTPCVolGroup> > edges, int maxNo=999);
      /// Connect a group of vertices to a group of track ends, holding on to the order of the hits in the path)
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >  ConnectVertexGroupsOrdered(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices, std::vector< ND::THandle<ND::TTPCVolGroup> > edges, int maxNo=999);

      /// Connect pairs of groups to each other
      std::vector< ND::THandle<ND::TTPCVolGroup> > ConnectGroups(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> >, bool allConnections=false, int maxNo=999);
      /// Connect pairs of groups to each other, holding on to the order of the hits in the path
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ConnectGroupsOrdered(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, bool vertices=false, bool allConnections=false, int maxNo=999);

      /// Clear any groups that lie too close to paths between other groups
      std::vector< ND::THandle<ND::TTPCVolGroup> > ClearRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, int maxNo=999);
      /// Experimental function to clear any groups that lie too close to paths between other groups
      std::vector< ND::THandle<ND::TTPCVolGroup> > ExperimentalClearRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, bool forceReconstructable=false, int maxNo=999);
      /// Experimental function to clear any groups that lie too close to paths between other groups
      std::vector< ND::THandle<ND::TTPCVolGroup> > ExperimentalClearRedundancies2(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > groups, int maxNo=999);
      /// Experimental function to clear any vertices that lie too close together
      std::vector< ND::THandle<ND::TTPCVolGroup> > ExperimentalClearVertexRedundancies(std::vector< ND::THandle<ND::TTPCVolGroup> > groups, int maxNo=999);
      /// Clear any groups that pass too close to an intermediate vertex
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ClearVertexConnectionRedundancies(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices);

      /// Get cost to connect one group to another group
      float FindConnectionCost(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, bool fullASICPenalty=false, bool extendMode=false, bool reduced=false, float maxCost=-1.);
      /// Get cost to connect one group to a volume
      float FindConnectionCost(ND::THandle<ND::TTPCVolGroup> group, ND::TTPCUnitVolume* vol, bool fullASICPenalty=false, bool extendMode=false, bool reduced=false, float maxCost=-1.);
      /// Get cost to connect one volume to anothger volume
      float FindConnectionCost(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, bool fullASICPenalty=false, bool extendMode=false, bool reduced=false, float maxCost=-1.);

      /// Connect one group to another group
      ND::THandle<ND::TTPCOrderedVolGroup> ConnectGroupPair(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, bool vertexGroup1=false, bool vertexGroup2=false, bool fullASICPenalty=false, bool extendMode=false);


      /// Associate each hit with the best matched path through pathfinding (heavy-ish function - don't call more than necessary)
      void AssociateBestHits(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > inPaths, bool fullASICPenalty=false, bool extendMode=false, float maxCost=9999);
      /// Merge in hits primary group to their closest groups in the input list
      void MergeBestHits(ND::TTPCVolGroupMan* volGroupMan, ND::THandle<ND::TTPCVolGroup> inGroup, bool fullASICPenalty=false, bool extendMode=false, float maxCost=9999.);
      /// Merge in hits primary group to their closest groups in the input list
      void MergeBestHits(ND::TTPCVolGroupMan* volGroupMan, std::vector< ND::THandle<ND::TTPCVolGroup> > inGroups, bool fullASICPenalty=false, bool extendMode=false, float maxCost=9999.);


      /// Get iterator to beginning of list of points to use for path finding
      std::vector<ND::TTPCAStarPoint*>::iterator begin(){ return fAStarPoints.begin(); }
      /// Get iterator to end of list of points to use for path finding
      std::vector<ND::TTPCAStarPoint*>::iterator end(){ return fAStarPoints.end(); }

    private:
      /// Find connected cells and the cost of those connections for cell vol
      void GetNearHitConnections(ND::TTPCVolGroupMan* volGroupMan, ND::TTPCAStarPoint* point, bool extendMode=false);

      /// Get cost of a connection between two cells point1 and point2
      float GetConnectionCost(ND::TTPCAStarPoint* point1, ND::TTPCAStarPoint* point2, bool penalty=false);
      /// Get heuristic cost between cell point1 and target
      float GetHeuristicCost(ND::TTPCAStarPoint* point1, ND::TTPCAStarPoint* target);

      /// Set up all internal variables resulting from connecting two vols
      int DoConnection(ND::TTPCAStarPoint* pathVolStart, ND::TTPCAStarPoint* pathVolEnd, bool fullASICPenalty=false, bool extendMode=false, float maxCost=-1.);
      /// Merge nearby groups from their A* connection distance
      std::vector< ND::THandle<ND::TTPCVolGroup> > MergeGroupsAStar(std::vector< ND::THandle<ND::TTPCVolGroup> > groups, float mergeDist=0.);

      /// Get modified cost based on a cell's charge and ASIC status
      float GetModifiedCost(float cost, ND::TTPCAStarPoint* point1, ND::TTPCAStarPoint* point2, bool fullASICPenalty=false, bool extendMode=false);

      /// Reset cells for this iteration without losing information which may need recalculating
      void RebootHits();

      /// Default scale for one cell apart connection in x direction
      float fXScale;
      /// Default scale for one cell apart connection in y direction
      float fYScale;
      /// Default scale for one cell apart connection in z direction
      float fZScale;

      /// Factor to weight heuristic by, to alter the performance and effectiveness of the algorithm
      float fHeuristicFactor;

      /// Layout associated with the cells used by this object
      ND::TTPCLayout* fLayout;
      /// Map of the hits to use for path finding and their associated unique ids
      std::map<long, ND::TTPCUnitVolume*> fHitMap;
      /// List of points to use for path finding
      std::vector<ND::TTPCAStarPoint*> fAStarPoints;
  };
}

#endif
