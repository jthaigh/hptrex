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

// ROOT
#include <TVector3.h>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCPathVolume.hxx"
#include "TTPCVolGroupMan.hxx"
#include "TTPCOrderedVolGroup.hxx"
#include "TTPCVolGroup.hxx"

namespace trex{
  /// Struct for holding elements for A* algorithm
  struct TTPCAStarPoint{
    trex::TTPCUnitVolume* vol;

    int x;
    int y;
    int z;

    float aStarCost;
    float aStarHeuristic;

    bool aStarHasFriends;
    bool aStarClosed;
    bool aStarOpen;

    trex::TTPCAStarPoint* aStarParent;
    std::map<trex::TTPCAStarPoint*, float> aStarFriends;
  };
  /// Performs an AStar algorithm to connect two points in a selection via a chain of valid cells
  class TTPCAStar{
    public:
      /// Default constructor
      TTPCAStar(trex::TTPCLayout* layout);
      /// Default destructor
      virtual ~TTPCAStar();

      /// Add the hits to use for the path finding
      void AddHits(trex::TTPCVolGroupMan* volGroupMan, std::map<long, trex::TTPCUnitVolume*> hitMap, bool extendMode=false);
      /// Add hits, saving time by copying from a previous A* container
      void AddHits(trex::TTPCAStar* prevAStar, std::map<long, trex::TTPCUnitVolume*> hitMap);



      /// Connect a group of vertices to a group of track ends, holding on to the order of the hits in the path)
    void ConnectVertexGroupsOrdered(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& vertices, std::vector< trex::TTPCVolGroup >& edges, std::vector< trex::TTPCOrderedVolGroup >& connections, int maxNo=999);

      /// Connect pairs of groups to each other, holding on to the order of the hits in the path
    void ConnectGroupsOrdered(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& groups, std::vector< trex::TTPCOrderedVolGroup >& connections, bool vertices=false, bool allConnections=false, int maxNo=999);
      
      /// Clear any groups that lie too close to paths between other groups
    void ClearRedundancies(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& groups, int maxNo=999);

    void ClearVertexConnectionRedundancies(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCOrderedVolGroup >& paths, std::vector< trex::TTPCVolGroup >& vertices);

      /// Get cost to connect one group to another group
      float FindConnectionCost(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, bool extendMode=false, bool reduced=false, float maxCost=-1.);
      /// Get cost to connect one group to a volume
      float FindConnectionCost(trex::TTPCVolGroup& group, trex::TTPCUnitVolume* vol, bool extendMode=false, bool reduced=false, float maxCost=-1.);
      /// Get cost to connect one volume to anothger volume
      float FindConnectionCost(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, bool extendMode=false, bool reduced=false, float maxCost=-1.);

      /// Connect one group to another group
    void ConnectGroupPair(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, trex::TTPCOrderedVolGroup& connection, bool vertexGroup1=false, bool vertexGroup2=false, bool extendMode=false);


      /// Associate each hit with the best matched path through pathfinding (heavy-ish function - don't call more than necessary)
      void AssociateBestHits(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCOrderedVolGroup >& inPaths, bool extendMode=false, float maxCost=9999);

      /// Merge in hits primary group to their closest groups in the input list
      void MergeBestHits(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::TTPCVolGroup >& inGroups, bool extendMode=false, float maxCost=9999.);


      /// Get iterator to beginning of list of points to use for path finding
      std::vector<trex::TTPCAStarPoint*>::iterator begin(){ return fAStarPoints.begin(); }
      /// Get iterator to end of list of points to use for path finding
      std::vector<trex::TTPCAStarPoint*>::iterator end(){ return fAStarPoints.end(); }

    private:
      /// Find connected cells and the cost of those connections for cell vol
      void GetNearHitConnections(trex::TTPCVolGroupMan* volGroupMan, trex::TTPCAStarPoint* point, bool extendMode=false);

      /// Get cost of a connection between two cells point1 and point2
      float GetConnectionCost(trex::TTPCAStarPoint* point1, trex::TTPCAStarPoint* point2);
      /// Get heuristic cost between cell point1 and target
      float GetHeuristicCost(trex::TTPCAStarPoint* point1, trex::TTPCAStarPoint* target);

      /// Set up all internal variables resulting from connecting two vols
      int DoConnection(trex::TTPCAStarPoint* pathVolStart, trex::TTPCAStarPoint* pathVolEnd, bool extendMode=false, float maxCost=-1.);
      
      /// Get modified cost based on a cell's charge and ASIC status
      float GetModifiedCost(float cost, trex::TTPCAStarPoint* point1, trex::TTPCAStarPoint* point2, bool extendMode=false);

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
      trex::TTPCLayout* fLayout;
      /// Map of the hits to use for path finding and their associated unique ids
      std::map<long, trex::TTPCUnitVolume*> fHitMap;
      /// List of points to use for path finding
      std::vector<trex::TTPCAStarPoint*> fAStarPoints;

  public:
      /// Connect a group of vertices to a group of track ends
      //MDH
      //Not used
      //std::vector< trex::THandle<trex::TTPCVolGroup> > ConnectVertexGroups(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::THandle<trex::TTPCVolGroup> > vertices, std::vector< trex::THandle<trex::TTPCVolGroup> > edges, int maxNo=999);
      /// Connect pairs of groups to each other
      //MDH
      //Not used
      //std::vector< trex::THandle<trex::TTPCVolGroup> > ConnectGroups(trex::TTPCVolGroupMan* volGroupMan, std::vector< trex::THandle<trex::TTPCVolGroup> >, bool allConnections=false, int maxNo=999);
//MDH
    //Not used
    /// Merge nearby groups from their A* connection distance
    //std::vector< trex::THandle<trex::TTPCVolGroup> > MergeGroupsAStar(std::vector< trex::THandle<trex::TTPCVolGroup> > groups, float mergeDist=0.);

    

  };
}

#endif
