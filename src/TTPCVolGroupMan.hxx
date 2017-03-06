#ifndef TTPCVOLGROUPMAN_HXX
#define TTPCVOLGROUPMAN_HXX

// c++
#include <iostream>
#include <utility>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <cmath>

// ROOT
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCPathVolume.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroup.hxx"
#include "TTPCOrderedVolGroup.hxx"

namespace trex{
  /// Enum for filters when grouping hits
  namespace TTPCHitGroupings{
    enum Type {
      all,
      delta,
      nonDelta,
      diff
    };
  };
  /// Class for holding and retrieving groups of TPC pattern recognition volume elements, and cleaning up once deleted
  class TTPCVolGroupMan{
    public:
      /// Constructor
      TTPCVolGroupMan(trex::TTPCLayout* layout);
      /// Destructor
      virtual ~TTPCVolGroupMan(){};

      /// Add the hits to the main group of hits for this manager
      void AddPrimaryHits(std::map<long, trex::TTPCUnitVolume*>& hitMap);

    //MDH
    //Not used
      /// Group delta hits in a two dimensional projection, then extend to 3D
      //std::vector< trex::THandle<trex::TTPCVolGroup> > GroupDeltaHits(trex::THandle<trex::TTPCVolGroup> suspectHits);
      /// Grab one group of delta hits, removing them from an input group
      //trex::THandle<trex::TTPCVolGroup> GrabADeltaHitGroup(trex::THandle<trex::TTPCVolGroup> inputGroup, trex::THandle<trex::TTPCVolGroup> allHits);

      /// Get vector of groups at edge of volumes split into delta and non-delt
      std::vector< trex::TTPCVolGroup > GetAllEdgeGroups();
      /// Get vector of groups on the edge of the MM volume
      std::vector< trex::TTPCVolGroup > GetEdgeGroups();
      /// Get vector of groups on the edge of the MM volume
      std::vector< trex::TTPCVolGroup > GetEdgeGroups(trex::TTPCVolGroup& inGroup, bool fiddlyLeans=false);

      /// Find points along set of paths at which they start to diverge
      void GetFoci(std::vector< trex::TTPCOrderedVolGroup >& paths, std::vector< trex::TTPCVolGroup >& groupsOut, float diffThreshold=0.5);

      /// Find first point along a path that's more than a given distance from a second path
      trex::TTPCPathVolume* GetDetatchmentPoint(trex::TTPCOrderedVolGroup& path1, trex::TTPCOrderedVolGroup& path2, trex::TTPCConnection::Type type=trex::TTPCConnection::vertexFind);
      /// Find point backwards along a path that's a given distance from the provided point on the same path
      trex::TTPCPathVolume* GetProjectionPoint(trex::TTPCOrderedVolGroup& path, trex::TTPCPathVolume* pathVol, int dist=0);
      /// Find closest point between two lines spanning between start and end points
      trex::TTPCUnitVolume* GetClosestPoint(trex::TTPCPathVolume* begin1, trex::TTPCPathVolume* end1, trex::TTPCPathVolume* begin2, trex::TTPCPathVolume* end2);

      /// Make sure a sensible number of paths and vertices is returned
      void CleanUpVertices(std::vector< trex::TTPCVolGroup >& edgeGroups, std::vector< trex::TTPCVolGroup >& vertices);
      

      /// Work back along two paths to find point where they started to diverge 
      trex::TTPCUnitVolume* GetZero(trex::TTPCOrderedVolGroup& path1, trex::TTPCOrderedVolGroup& path2, trex::TTPCPathVolume* vol1=0, trex::TTPCPathVolume* vol2=0);
      /// Work back along first path to find point where it started to diverge from second
      trex::TTPCPathVolume* GetPathZero(trex::TTPCOrderedVolGroup& path1, trex::TTPCOrderedVolGroup& path2, trex::TTPCPathVolume* vol=0);
      /// Get area of a triangle between three points
      double GetTriangleArea(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, trex::TTPCUnitVolume* vol3);
      /// Get distance from one point in a path to the closest point in another
      double GetMinDistance(trex::TTPCOrderedVolGroup& path, trex::TTPCPathVolume* point);
      /// Get distance from one unit volume to the closest point in another
      double GetMinDistance(trex::TTPCOrderedVolGroup& path, trex::TTPCUnitVolume* vol);

      /// Get preliminary group of hits which could be a track end
      void GetFarHitsPreGroup(trex::TTPCOrderedVolGroup& path, trex::TTPCVolGroup& hits);
      /// Get group of hits which could be a track end
      void GetFarHitsGroup(trex::TTPCOrderedVolGroup& path, trex::TTPCVolGroup& hits);

      /// Determine whether a vol overlaps a path 
      bool GetPathVolOverlap(trex::TTPCOrderedVolGroup& path, trex::TTPCUnitVolume* vol, trex::TTPCConnection::Type type = trex::TTPCConnection::edgeMerge);
      /// Determine whether the average position of one group comes near another
      bool GetGroupGroupOverlap(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, trex::TTPCConnection::Type type = trex::TTPCConnection::edgeMerge, bool simple=false, bool checkLean=false);

      /// Merge two groups and return the result
      void MergeGroups(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2,trex::TTPCVolGroup& merged);

      /// Add anything inside this group to it
      void BulkGroup(trex::TTPCVolGroup& group);
      /// Add anything inside each group to it
      void BulkGroups(std::vector< trex::TTPCVolGroup >& groups);

      /// Break input tracks around any kinks found inside them
      void BreakPathsAboutKinks(std::vector< trex::TTPCOrderedVolGroup >& paths);

      /// Get a pointer to the main group of hits associated with this manager
      trex::TTPCVolGroup& GetPrimaryHits(){ return fPrimaryHits; }

      /// Get unordered group corresponding to hits in an ordered group
    void GetUnorderedGroup(trex::TTPCOrderedVolGroup& in, trex::TTPCVolGroup& out);


      /// Process hits in unordered group, associating an unordered group of nearby hits with it
      void BuildGroupFriends(trex::TTPCOrderedVolGroup& in, trex::TTPCConnection::Type type = trex::TTPCConnection::pathHits);
      /// Produce clusters in unordered group
      void ClusterGroupFriends(trex::TTPCOrderedVolGroup& in, bool doClustering=false, bool checkX=false, bool partial=false);
      
    //MDH
    //Not used
/// Return copy of list of paths with cleared empty paths
    //std::vector< trex::THandle<trex::TTPCOrderedVolGroup> > ClearEmptyPaths(std::vector< trex::THandle<trex::TTPCOrderedVolGroup> >& paths);

      /// Uniquely build friends for all paths
      void BuildAllFriends(std::vector< trex::TTPCOrderedVolGroup>& paths);
      /// Separate hits assigned to incomplete paths so that each contains a unique set
      void SeparateHits(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Merge x-paths into junctions
      void SeparateXPathHits(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Separate out any empty hits
      void SeparateEmptyClusters(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Separate hits mis-assigned to one cluster from a nearby one
      void SeparateClusterHits(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Separate hits next to a junction which don't seem to fit in to the rest of the path
      void SeparateAnomHits(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Separate hits assigned to junctions from those assigned to complete paths
      std::vector<trex::TTPCUnitVolume*> SeparateAnomHitsPath(trex::TTPCOrderedVolGroup& path, int checkDir=0);
      /// Separate hits assigned to junctions from those assigned to complete paths
      void SeparateJunctionHits(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Enforce ordering on paths
      void EnforceOrdering(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Get whether two vol groups overlap
      bool GetOverlaps(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, bool checkHits=true);

      /// Set vertex statuses exclusively based on how many paths a junction is connected to
      void ResetVertexStatuses(std::vector< trex::TTPCOrderedVolGroup >& paths, bool partial=false);

      /// Get extended group corresponding to hits in given group, and hits from hitMap which are close enough to be incorporated
      void GetExtendedGroup(trex::TTPCVolGroup& in, trex::TTPCVolGroup& out);

      /// Get hits in fPrimaryHits connected to this one within specified distance
    void GetConnectedHits(std::vector<trex::TTPCVolGroup>& out,trex::TTPCConnection::Type type=trex::TTPCConnection::path, trex::TTPCHitGroupings::Type typeFilter=trex::TTPCHitGroupings::all, bool usabilityCheck=false);

      /// Get hits in given group connected to this one within specified distance
    void GetConnectedHits(trex::TTPCVolGroup& in, std::vector<trex::TTPCVolGroup>& out, trex::TTPCConnection::Type type=trex::TTPCConnection::path, trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all, bool usabilityCheck=false);

      /// Check that a provided group covers the minimum number of pads
      bool CheckUsability(trex::TTPCVolGroup& inGroup);

      /// Get cells in an input group within specified distance of a specified id (including itself if inclusive is set to true)
      void GetNearHits(trex::TTPCVolGroup& in, trex::TTPCVolGroup& cellsOut, long id, trex::TTPCConnection::Type type=trex::TTPCConnection::path, bool inclusive=false, bool singular=false, float distFilter=-2., trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all, bool square=false);
      /// Get cells in an input group within specified distance of a specified cell (including itself if inclusive is set to true)
      void GetNearHits(trex::TTPCVolGroup& in, trex::TTPCVolGroup& cellsOut, trex::TTPCUnitVolume* vol, trex::TTPCConnection::Type type=trex::TTPCConnection::path, bool inclusive=false, bool singular=false, float distFilter=-1., trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all, bool square=false);

      /// Get list of tagged vertices from a list of paths
    std::vector< trex::TTPCVolGroup* > GetJunctionsFromPaths(std::vector< trex::TTPCOrderedVolGroup >& paths);
      /// Get unused hits given an input of paths and junctions
    void GetUnusedHits(std::vector< trex::TTPCOrderedVolGroup >& paths, TTPCVolGroup& unusedHits);

      /// Associate unused hits with path junctions
    void AssociateUnusedWithJunctions(trex::TTPCVolGroup& unused, 
				      std::vector< trex::TTPCOrderedVolGroup >& paths,
				      std::vector< trex::TTPCVolGroup >& junctionsToAddTo);

      /// Filter out any remaining tracks and hits that don't make sense
      void SanityFilter(std::vector< trex::TTPCOrderedVolGroup >& input);
   
      /// Get representation of position using co-ordinates
      TVector3 GetAvgPosRep(trex::TTPCPathVolume* vol);
      /// Get representation of position using co-ordinates
      TVector3 GetAvgPosRep(trex::TTPCUnitVolume* vol, int sign=0);

      /// Get whether a specified cell is within a specified range of another
      bool IsInRange(trex::TTPCPathVolume* point1, trex::TTPCPathVolume* point2, int sizeX=1, int sizeY=1, int sizeZ=1);
      /// Get whether a specified path meets the criteria for being tagged as an x-path
      bool IsXPathCandidate(trex::TTPCOrderedVolGroup& inPath);

      /// Get hits associated with primary path
      std::vector<trex::TTPCHitPad*> GetHits();

      
    private:
      /// TPC layout associated with this group
      trex::TTPCLayout* fLayout;

      /// Primary map of all hits associated with this manager 
      std::map<long, trex::TTPCUnitVolume*> fHitMap; 
      /// Primary group of hits assoicated with this manager
      trex::TTPCVolGroup fPrimaryHits;

      /// Recursively build groups and add them to a given container
      void FillWithSplitGroups(std::vector< trex::TTPCVolGroup >& container, trex::TTPCVolGroup& inputHits, trex::TTPCConnection::Type type=trex::TTPCConnection::path, int maxFilterX=0, int maxFilterY=0, int maxFilterZ=0);
      /// Retrieve split groups from one input
    std::vector< trex::TTPCVolGroup > GetSplitGroups(trex::TTPCVolGroup& inputHits, trex::TTPCConnection::Type type=trex::TTPCConnection::path);

      /// Build a list of connected cells recursively, with type specifying mode
      void RecursiveFriendBuild(long startID, trex::TTPCVolGroup& target, trex::TTPCVolGroup& source, trex::TTPCConnection::Type type=trex::TTPCConnection::path, trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all);
            
    /// Get list of connected cells within a certain distance of a specified list recursively
      void RecursiveFriendSeek(trex::TTPCOrderedVolGroup& inList, trex::TTPCVolGroup& target, float dist=9999., trex::TTPCConnection::Type type=trex::TTPCConnection::path, trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all);
      /// Get list of connected cells within a certain distance of a specified list recursively
      void RecursiveFriendSeek(trex::TTPCVolGroup& inList, trex::TTPCVolGroup& target, float dist=9999., trex::TTPCConnection::Type type=trex::TTPCConnection::path, trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all);
      /// Fill distances for cells in input list to start within specified distance 
      void RecursiveFriendListSeek(long startID, trex::TTPCVolGroup& source, float curDist=9999., float dist=9999., trex::TTPCConnection::Type type=trex::TTPCConnection::path, trex::TTPCHitGroupings::Type typeFilter = trex::TTPCHitGroupings::all);

      /// Check if two unit volumes are within a given distance of each other
      bool HitTest(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, trex::TTPCConnection::Type type);
      /// Check if two unit volumes are within a given distance of each other
      bool HitTest(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, int xDist, int yDist, int zDist);
  };
}

#endif
