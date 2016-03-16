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

// nd280
#include <TND280Log.hxx>
#include <THandle.hxx>

// eddy
#include "TTPCLayout.hxx"
#include "TTPCPathVolume.hxx"
#include "TTPCUnitVolume.hxx"
#include "TTPCVolGroup.hxx"
#include "TTPCOrderedVolGroup.hxx"

namespace ND{
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
      TTPCVolGroupMan(ND::TTPCLayout* layout);
      /// Destructor
      virtual ~TTPCVolGroupMan(){};

      /// Add the hits to the main group of hits for this manager
      void AddPrimaryHits(std::map<long, ND::TTPCUnitVolume*> hitMap);

      /// Group delta hits in a two dimensional projection, then extend to 3D
      std::vector< ND::THandle<ND::TTPCVolGroup> > GroupDeltaHits(ND::THandle<ND::TTPCVolGroup> suspectHits);
      /// Grab one group of delta hits, removing them from an input group
      ND::THandle<ND::TTPCVolGroup> GrabADeltaHitGroup(ND::THandle<ND::TTPCVolGroup> inputGroup, ND::THandle<ND::TTPCVolGroup> allHits);

      /// Get vector of groups at edge of volumes split into delta and non-delt
      void GetAllEdgeGroups(std::vector< ND::TTPCVolGroup >& output);
      /// Get vector of groups on the edge of the MM volume
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetEdgeGroups();
      /// Get vector of groups on the edge of the MM volume
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetEdgeGroups(ND::THandle<ND::TTPCVolGroup> inGroup, bool fiddlyLeans=false);


      /// Find points along set of paths at which they start to diverge
      void GetFoci(std::vector< ND::TTPCOrderedVolGroup > paths, std::vector< ND::TTPCVolGroup > groupsOut, float diffThreshold=0.5);

      /// Find first point along a path that's more than a given distance from a second path
      ND::TTPCPathVolume* GetDetatchmentPoint(ND::THandle<ND::TTPCOrderedVolGroup> path1, ND::THandle<ND::TTPCOrderedVolGroup> path2, ND::TTPCConnection::Type type=ND::TTPCConnection::vertexFind);
      /// Find point backwards along a path that's a given distance from the provided point on the same path
      ND::TTPCPathVolume* GetProjectionPoint(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCPathVolume* pathVol, int dist=0);
      /// Find closest point between two lines spanning between start and end points
      ND::TTPCUnitVolume* GetClosestPoint(ND::TTPCPathVolume* begin1, ND::TTPCPathVolume* end1, ND::TTPCPathVolume* begin2, ND::TTPCPathVolume* end2);

      /// Make sure a sensible number of paths and vertices is returned
      void CleanUpVertices(std::vector< ND::TTPCVolGroup > edgeGroups, std::vector< ND::TTPCVolGroup > vertices);
      

      /// Work back along two paths to find point where they started to diverge 
      ND::TTPCUnitVolume* GetZero(ND::THandle<ND::TTPCOrderedVolGroup> path1, ND::THandle<ND::TTPCOrderedVolGroup> path2, ND::TTPCPathVolume* vol1=0, ND::TTPCPathVolume* vol2=0);
      /// Work back along first path to find point where it started to diverge from second
      ND::TTPCPathVolume* GetPathZero(ND::THandle<ND::TTPCOrderedVolGroup> path1, ND::THandle<ND::TTPCOrderedVolGroup> path2, ND::TTPCPathVolume* vol=0);
      /// Get area of a triangle between three points
      double GetTriangleArea(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, ND::TTPCUnitVolume* vol3);
      /// Get distance from one point in a path to the closest point in another
      double GetMinDistance(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCPathVolume* point);
      /// Get distance from one unit volume to the closest point in another
      double GetMinDistance(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCUnitVolume* vol);

      /// Get preliminary group of hits which could be a track end
      ND::THandle<ND::TTPCVolGroup> GetFarHitsPreGroup(ND::THandle<ND::TTPCOrderedVolGroup> path, bool tryChargeCut=false);
      /// Get group of hits which could be a track end
      ND::THandle<ND::TTPCVolGroup> GetFarHitsGroup(ND::THandle<ND::TTPCOrderedVolGroup> path, bool tryChargeCut=false);
      /// Find one point along set of paths at which discontinuities occur (the most significant)
      ND::THandle<ND::TTPCVolGroup> GetDiscontinuity(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Find points along set of paths at which discontinuities occur
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetDiscontinuities(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, float diffThreshold=.5, float threshold=2.);
      /// Get in-path id of a sudden sharp change 
      int GetDiscontinuityID(ND::THandle<ND::TTPCOrderedVolGroup> path, int dimension=1, float threshold=2.);
      /// Find suddent sharp change in angle along a path 
      void FindDiscontinuity(float& pos, float& step, float size, TH1F* hist, float threshold=2.);

      /// Determine whether a vol overlaps a path 
      bool GetPathVolOverlap(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCUnitVolume* vol, ND::TTPCConnection::Type type = ND::TTPCConnection::edgeMerge);
      /// Determine whether the average position of one group comes near another
      bool GetGroupGroupOverlap(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, ND::TTPCConnection::Type type = ND::TTPCConnection::edgeMerge, bool simple=false, bool checkLean=false);

      /// Merge two groups and return the result
      ND::THandle<ND::TTPCVolGroup> MergeGroups(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2);

      /// Add anything inside this group to it
      void BulkGroup(ND::TTPCVolGroup& group);
      /// Add anything inside each group to it
      void BulkGroups(std::vector< ND::TTPCVolGroup >& groups);

      /// Break input tracks around any kinks found inside them
      void ND::TTPCVolGroupMan::BreakPathsAboutKinks(std::vector< ND::TTPCOrderedVolGroup >& paths, bool tryChargeCut=false);

      /// Get a pointer to the main group of hits associated with this manager
      ND::THandle<ND::TTPCVolGroup> GetPrimaryHits(){ return fPrimaryHits; }

      /// Get unordered group corresponding to hits in an ordered group
      ND::THandle<ND::TTPCVolGroup> GetUnorderedGroup(ND::THandle<ND::TTPCOrderedVolGroup> in);


      /// Process hits in unordered group, associating an unordered group of nearby hits with it
      void BuildGroupFriends(ND::TTPCOrderedVolGroup& in, ND::TTPCConnection::Type type);
      /// Produce clusters in unordered group
      void ClusterGroupFriends(ND::TTPCOrderedVolGroup& in, bool doClustering=false, bool checkX=false, bool partial=false);
      /// Return copy of list of paths with cleared empty paths
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ClearEmptyPaths(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& paths);
      /// Uniquely build friends for all paths
      void BuildAllFriends(std::vector< ND::TTPCOrderedVolGroup>& paths);
      /// Separate hits assigned to incomplete paths so that each contains a unique set
      void SeparateHits(std::vector< ND::TTPCOrderedVolGroup >& paths);
      /// Merge x-paths into junctions
      void SeparateXPathHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& paths);
      /// Separate out any empty hits
      void SeparateEmptyClusters(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Separate hits mis-assigned to one cluster from a nearby one
      void SeparateClusterHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Separate hits next to a junction which don't seem to fit in to the rest of the path
      void SeparateAnomHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Separate hits assigned to junctions from those assigned to complete paths
      std::vector<ND::TTPCUnitVolume*> SeparateAnomHitsPath(ND::THandle<ND::TTPCOrderedVolGroup> path, int checkDir=0);
      /// Separate hits assigned to junctions from those assigned to complete paths
      void SeparateJunctionHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Enforce ordering on paths
      void EnforceOrdering(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Get whether two vol groups overlap
      bool GetOverlaps(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, bool checkHits=true);

      /// Set vertex statuses exclusively based on how many paths a junction is connected to
      void ResetVertexStatuses(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, bool partial=false);

      /// Get extended group corresponding to hits in given group, and hits from hitMap which are close enough to be incorporated
      ND::THandle<ND::TTPCVolGroup> GetExtendedGroup(ND::THandle<ND::TTPCVolGroup> in);


      /// Get hits in fPrimaryHits connected to this one within specified distance
    void GetConnectedHits(std::vector<ND::TTPCVolGroup>& out, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter, bool usabilityCheck);

      /// Get hits in given group connected to this one within specified distance
    void GetConnectedHits(ND::TTPCVolGroup& in, std::vector<ND::TTPCVolGroup>& out, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter, bool usabilityCheck);

      /// Check that a provided group covers the minimum number of pads
      bool CheckUsability(ND::THandle<ND::TTPCVolGroup> inGroup);

      /// Get cells in an input group within specified distance of a specified id (including itself if inclusive is set to true)
      ND::THandle<ND::TTPCVolGroup> GetNearHits(ND::THandle<ND::TTPCVolGroup> in, long id, ND::TTPCConnection::Type type=ND::TTPCConnection::path, bool inclusive=false, bool singular=false, float distFilter=-2., ND::TTPCHitGroupings::Type typeFilter = ND::TTPCHitGroupings::all, bool square=false);
      /// Get cells in an input group within specified distance of a specified cell (including itself if inclusive is set to true)
      ND::THandle<ND::TTPCVolGroup> GetNearHits(ND::THandle<ND::TTPCVolGroup> in, ND::TTPCUnitVolume* vol, ND::TTPCConnection::Type type=ND::TTPCConnection::path, bool inclusive=false, bool singular=false, float distFilter=-1., ND::TTPCHitGroupings::Type typeFilter = ND::TTPCHitGroupings::all, bool square=false);

      /// Get list of tagged vertices from a list of paths
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetJunctionsFromPaths(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);
      /// Get unused hits given an input of paths and junctions
    void GetUnusedHits(std::vector< ND::TTPCOrderedVolGroup >& paths, TTPCVolGroup& unusedHits);

      /// Associate unused hits with path junctions
    void ND::TTPCVolGroupMan::AssociateUnusedWithJunctions(ND::TTPCVolGroup& unused, std::vector< ND::TTPCOrderedVolGroup >& paths);

      /// Filter out any remaining tracks and hits that don't make sense
      void SanityFilter(std::vector< ND::TTPCOrderedVolGroup >& input);
   
      /// Get representation of position using time and co-ordinates
      TVector3 GetAvgPosRep(ND::TTPCPathVolume* vol);
      /// Get representation of position using time and co-ordinates
      TVector3 GetAvgPosRep(ND::TTPCUnitVolume* vol, int sign=0);

      /// Get whether a specified cell is within a specified range of another
      bool IsInRange(ND::TTPCPathVolume* point1, ND::TTPCPathVolume* point2, int sizeX=1, int sizeY=1, int sizeZ=1);
      /// Get whether a specified path meets the criteria for being tagged as an x-path
      bool IsXPathCandidate(ND::THandle<ND::TTPCOrderedVolGroup> inPath);

      /// Get hits associated with primary path
      ND::THandle<ND::THitSelection> GetHits();

      
    private:
      /// TPC layout associated with this group
      ND::TTPCLayout* fLayout;

      /// Primary map of all hits associated with this manager 
      std::map<long, ND::TTPCUnitVolume*> fHitMap; 
      /// Primary group of hits assoicated with this manager
      ND::TTPCVolGroup fPrimaryHits;

      /// Recursively build groups and add them to a given container
      void FillWithSplitGroups(std::vector< ND::THandle<ND::TTPCVolGroup> >& container, ND::THandle<ND::TTPCVolGroup> inputHits, ND::TTPCConnection::Type type=ND::TTPCConnection::path, int maxFilterX=0, int maxFilterY=0, int maxFilterZ=0);
      /// Retrieve split groups from one input
      std::vector< ND::THandle<ND::TTPCVolGroup> > GetSplitGroups(ND::THandle<ND::TTPCVolGroup> inputHits, ND::TTPCConnection::Type type=ND::TTPCConnection::path);

      /// Build a list of connected cells recursively, with type specifying mode
      void RecursiveFriendBuild(long startID, ND::THandle<ND::TTPCVolGroup> target, ND::THandle<ND::TTPCVolGroup> source, ND::TTPCConnection::Type type=ND::TTPCConnection::path, ND::TTPCHitGroupings::Type typeFilter = ND::TTPCHitGroupings::all);
            
    /// Get list of connected cells within a certain distance of a specified list recursively
      void RecursiveFriendSeek(ND::THandle<ND::TTPCOrderedVolGroup> inList, ND::THandle<ND::TTPCVolGroup> target, float dist=9999., ND::TTPCConnection::Type type=ND::TTPCConnection::path, ND::TTPCHitGroupings::Type typeFilter = ND::TTPCHitGroupings::all);
      /// Get list of connected cells within a certain distance of a specified list recursively
      void RecursiveFriendSeek(ND::THandle<ND::TTPCVolGroup> inList, ND::THandle<ND::TTPCVolGroup> target, float dist=9999., ND::TTPCConnection::Type type=ND::TTPCConnection::path, ND::TTPCHitGroupings::Type typeFilter = ND::TTPCHitGroupings::all);
      /// Fill distances for cells in input list to start within specified distance 
      void RecursiveFriendListSeek(long startID, ND::THandle<ND::TTPCVolGroup> source, float curDist=9999., float dist=9999., ND::TTPCConnection::Type type=ND::TTPCConnection::path, ND::TTPCHitGroupings::Type typeFilter = ND::TTPCHitGroupings::all);

      /// Check if two unit volumes are within a given distance of each other
      bool HitTest(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, ND::TTPCConnection::Type type);
      /// Check if two unit volumes are within a given distance of each other
      bool HitTest(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, int xDist, int yDist, int zDist);

  public:
    //MDH
    //Not used
    /// Separate an input into its delta and non-delta hits
      //void GetDeltaNonDelta(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& nonDelta, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& delta, ND::THandle<ND::TTPCOrderedVolGroup> input);
      /// Break long junctions into x-paths
      //void BreakLongJunctions(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& paths);

      /// Count any duplicated hits in paths and junctions
      //std::pair<std::pair<int, int>, std::pair<int, int> > CountDuplicates(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);

      /// Get iterable vector to ids in 3D spheroid of hits in map, given input connection scheme
      //std::vector<long> SpheroidIterable(ND::TTPCCell3D cell, ND::TTPCConnection::Type type);
      /// Get iterable vector to ids in 3D spheroid of hits in map, given input dimensions in x, y and z
      //std::vector<long> SpheroidIterable(ND::TTPCCell3D cell, int diffX, int diffY, int diffZ);
    /// Debugging function to check overlap between path extended hits
    //      void CheckExtendedHitOverlap(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);

    //MDH
    //Not used
      /// Try to fill in gaps between MM volumes with pseudo-hits to tape over difficulties in pattern recognition
      //ND::THandle<ND::TTPCVolGroup> GetGapFillHits(ND::THandle<ND::TTPCVolGroup> in);
      /// Extrapolate this group of hits in a particular direction based on larger set of reference hits
      //ND::THandle<ND::TTPCVolGroup> GetGapProjections(ND::THandle<ND::TTPCVolGroup> in, int dirX, int dirY, int dirZ, std::map<long, ND::TTPCUnitVolume*> refHits);

      /// Find the nearest group of hits to the provided one along a given axis
      //ND::THandle<ND::TTPCVolGroup> GetNearestGroup(ND::THandle<ND::TTPCVolGroup> in, int x, int y, int z, int axis, int maxDist=100);
      /// Find group containing just nearest hit to a given x, y and z from input
      //ND::THandle<ND::TTPCVolGroup> GetNearestHitGroup(ND::THandle<ND::TTPCVolGroup> in, int x, int y, int z, int maxDist=100);

      /// Find group representing portion of a pair of groups which overlaps (or doesn't if swap=true)
      //ND::THandle<ND::TTPCVolGroup> GetOverlap(ND::THandle<ND::TTPCVolGroup> grp1, ND::THandle<ND::TTPCVolGroup> grp2, bool swap=false);
    //MDH
    //Not used
      /// Find one point along set of paths at which they start to diverge (the most significant)
      //ND::THandle<ND::TTPCVolGroup> GetFocus(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths);

    //MDH
    //Not used
      /// Find points along set of paths which could be track ends
      //std::vector< ND::THandle<ND::TTPCVolGroup> > GetFarHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, bool tryChargeCut=false);

    //MDH
    //Not used
      /// Merge input tracks which severely overlap 
      //std::vector< ND::THandle<ND::TTPCVolGroup> > GetFilteredPaths(std::vector< ND::THandle<ND::TTPCVolGroup> > inPaths);

    //MDH
    //Not used
      /// Get extended unordered group corresponding to hits in given group
      //ND::THandle<ND::TTPCVolGroup> GetExtendedUnorderedGroup(ND::THandle<ND::TTPCOrderedVolGroup> in);

  private:

    //MDH
    //Not used
    /// Get list of connected cells within a certain distance of a specified list recursively
    //void FriendConnect(ND::THandle<ND::TTPCOrderedVolGroup> inList, ND::THandle<ND::TTPCVolGroup> inField);

  };
}

#endif
