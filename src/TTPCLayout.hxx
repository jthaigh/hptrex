#ifndef TTPCLAYOUT_HXX
#define TTPCLAYOUT_HXX

// c++
#include <iostream>
#include <map>
#include <algorithm>
#include <sstream>

// ROOT
#include <TVector3.h>

// TRex
#include <TTPCHitPad.hxx>

namespace trex{

  /// Enumeration for specifying types of connection (effects characteristic distances)
  // should work similar way to c++11 enum class
  namespace TTPCConnection{
    enum Type {
      path,
      pathHits,
      extraHits,
      edgeMerge,
      vertexMerge,
      vertexFind,
      vertexHit,
      vertexPath,
      clusterMerge
    };
  };
  /// Structure for 3D cell with x, y and z id
  struct TTPCCellInfo3D{
    int x;
    int y;
    int z;
  };

  /// Class containing relevant information and methods for expressing 3D layout of hits in a way that makes sense to feature finding and path finding algorithms
  class TTPCLayout{
    public:
      /// Default constructor
      TTPCLayout();
      /// Default destructor
      ~TTPCLayout();

      /// Set cell id ranges in x, y and z
      void SetRanges(int minX,int maxX, int minY,int maxY, int minZ,int maxZ);

      /// Get cell id in x, y and z based on a 3D position
      trex::TTPCCellInfo3D GetPadPosID(TVector3 pos, int tpcMask=-1);

      /// Convert cell id in x, y and z to a single unique long integer
      long Mash(int x, int y, int z);
      /// Convert cell id in x, y and z to a single unique long integer, returning -1 if bad ids are supplied
      long SafeMash(int x, int y, int z);

      bool GetHaveX(){ return fHaveX; }

      /// Find appropriate distances for a given type of connection
      void GetTypeDistances(int& distX, int& distY, int& distZ, trex::TTPCConnection::Type type);

      /// Get distance between MM pad centres in y direction
      double GetPadPitchY(){ return fPadPitchY; }
      /// Get distance between MM pad centres in z direction
      double GetPadPitchZ(){ return fPadPitchZ; }

      /// Get size of an individual x cel
      double GetXCellSize(){ return fXCellSize; }

      /// Get minimum number of pads constituting a useful pattern
      int GetMinPatternPads(){ return fMinPatternPads; }
      /// Get minimum number of clusters constituting a useful path
      int GetMinPathClusters(){ return fMinPathClusters; }

      /// Get number of layers at edge of sub event to search for track ends
      int GetEdgeLayers(){ return fEdgeLayers; }

      /// Default A* scale for one cell apart connection in x direction
      float GetAStarXScale(){ return fAStarXScale; }
      /// Default A* scale for one cell apart connection in y direction
      float GetAStarYScale(){ return fAStarYScale; }
      /// Default A* scale for one cell apart connection in z direction
      float GetAStarZScale(){ return fAStarZScale; }
      /// Factor to weight A* heuristic by, to alter the performance and effectiveness of the algorithm
      float GetAStarHeuristicFactor(){ return fAStarHeuristicFactor; }

      /// Get how close cells can be to connect for pattern recognition and path finding, in x
      int GetConnectDistX(){ return fConnectDistX; }
      /// Get how close cells can be to connect for pattern recognition and path finding, in y
      int GetConnectDistY(){ return fConnectDistY; }
      /// Get how close cells can be to connect for pattern recognition and path finding, in z
      int GetConnectDistZ(){ return fConnectDistZ; }
      /// Get scale for specifying rough structure for pattern recognition and path finding, in x
      int GetStructDistX(){ return fStructDistX; }
      /// Get scale for specifying rough structure for pattern recognition and path finding, in y
      int GetStructDistY(){ return fStructDistY; }
      /// Get scale for specifying rough structure for pattern recognition and path finding, in z
      int GetStructDistZ(){ return fStructDistZ; }

      /// Get distance for connecting hits through A*
      int GetAStarConnectDist(){ return fAStarConnectDist; }
      /// Get distance for associating hits with a path
      int GetPathHitConnectDist(){ return fPathHitConnectDist; }
      /// Get distance for identifying extra hits by distance from an existing path
      int GetExtraHitConnectDist(){ return fExtraHitConnectDist; }
      /// Get maximum for merging overlapping path edge groups
      int GetEdgeMergeStructDist(){ return fEdgeMergeStructDist; }
      /// Get maximum for merging overlapping vertex groups
      int GetVertexMergeStructDist(){ return fVertexMergeStructDist; }
      /// Get distance for finding vertices as points of divergence
      int GetVertexFindStructDist(){ return fVertexFindStructDist; }
      /// Get distance for adding hits to vertices
      int GetVertexHitStructDist(){ return fVertexHitStructDist; }
      /// Get distance for tagging a path as redundant when it comes close to a vertex
      int GetVertexPathStructDist(){ return fVertexPathStructDist; }
      /// Get distance for merging two nearby clusters into a vetex
      int GetClusterMergeStructDist(){ return fClusterMergeStructDist; }

      /// Get maximum extent in MM pads for a path counted as an x-path
      int GetXPathMaxPads(){ return fXPathMaxPads; }
      /// Get minimum ratio between x end-start position and y-z end-start position for a path counted as an x-path
      double GetXPathMinEndRatio(){ return fXPathMinEndRatio; }

      /// Get minumum distance to run precheck in for vertex detection
      double GetEdgePreDist(){ return fEdgePreDist; }
      /// Get threshold angle for precheck for vertex detection
      double GetEdgePreAng(){ return fEdgePreAng; }
      /// Get minumum distance to collect hits in for vertex detection
      double GetEdgeMinDist(){ return fEdgeMinDist; }
      /// Get minumum distance to compute offset with for vertex detection
      double GetEdgeMinDistLow(){ return fEdgeMinDistLow; }
      /// Get fractional minimum hits at the edge for vertex detection, to avoid false positives from deltas
      double GetEdgeMinHits(){ return fEdgeMinHits; }
      /// Get fractional minimum asymmetry for vertex detection, to avoid false positives from straight tracks
      double GetEdgeMaxSigma(){ return fEdgeMaxSigma; }
      /// Get fractional range to collect hits in for vertex detection
      double GetEdgeRange(){ return fEdgeRange; }
      /// Get threshold to classify turn as a vertex
      double GetEdgeThreshold(){ return fEdgeThreshold; }
      /// Get offset from center to classify turn as a vertex
      double GetEdgeOffset(){ return fEdgeOffset; }

      /// Get maximum distance for connecting hits when forming clusters
      int GetHVClusterExtrapolateDist(){ return fHVClusterExtrapolateDist; }
      /// Get maximum total distance for connecting hits when forming clusters
      int GetHVClusterExtrapolateLimit(){ return fHVClusterExtrapolateLimit; }
      /// Get distance for connecting hits when forming clusters
      int GetClusterConnectDist(){ return fClusterConnectDist; }
      /// Get maximum number of isolated horizontal or vertical cells at the edge of a HV cluster
      int GetHVEdgeDist(){ return fHVEdgeDist; }
      /// Get angle for discriminating between horizontal and vertical clusters
      float GetThresholdAngle(){ return fThresholdAngle; }
      /// Get minimum points along path segment required when establishing angle from dichotomy technique
      int GetDichotomyCutoff(){ return fDichotomyCutoff; }

      /// Get distance to check when cleaning up anomalous hits near the junction
      int GetAnomCheckDist(){ return fAnomCheckDist; }
      /// Get distance to use for control when cleaning up anomalous hits near the junction
      int GetAnomProjectDist(){ return fAnomProjectDist; }
      /// Get threshold when cleaning up anomalous hits near the junction
      double GetAnomMaxOffs(){ return fAnomMaxOffs; }

      /// Get maximum isolated clusters at the start of a path before they're merged into a nearby vertex
      int GetHVClusterMaxIso(){ return fHVClusterMaxIso; }

    TVector3 GetMinPos(){return TVector3(fMinX,fMinY,fMinZ);}

    TVector3 GetMaxPos(){return TVector3(fMaxX,fMaxY,fMaxZ);}

    double GetBField(){return fBField;}

    private:

      bool fHaveX;

      /// Distance between MM pad centres in y direction
      double fPadPitchY;
      // Distance between MM pad centres in z direction
      double fPadPitchZ;

      /// Default A* scale for one cell apart connection in x direction
      float fAStarXScale;
      /// Default A* scale for one cell apart connection in y direction
      float fAStarYScale;
      /// Default A* scale for one cell apart connection in z direction
      float fAStarZScale;
      /// Factor to weight A* heuristic by, to alter the performance and effectiveness of the algorithm
      float fAStarHeuristicFactor;

      /// Size of an individual x cell
      double fXCellSize;

      /// Minimum x id in event
      int fMinX;
      /// Minimum y id in event
      int fMinY;
      /// Minimum z id in event
      int fMinZ;
      /// Maximum x id in event
      int fMaxX;
      /// Maximum y id in event
      int fMaxY;
      /// Maximum z id in event
      int fMaxZ;
      /// Number of x ids in event
      int fSizeX;
      /// Number of y ids in event
      int fSizeY;
      /// Number of z ids in event
      int fSizeZ;

      /// Minimum useful size of pattern
      int fMinPatternPads;
      /// Minimum useful size of path
      int fMinPathClusters;

      /// Number of layers at edge of sub event to use when looking for edges
      int fEdgeLayers;

      /// How close cells can be to connect for pattern recognition and path finding, in x
      int fConnectDistX;
      /// How close cells can be to connect for pattern recognition and path finding, in y
      int fConnectDistY;
      /// How close cells can be to connect for pattern recognition and path finding, in z
      int fConnectDistZ;
      /// Scale for specifying rough structure for pattern recognition and path finding, in x
      int fStructDistX;
      /// Scale for specifying rough structure for pattern recognition and path finding, in y
      int fStructDistY;
      /// Scale for specifying rough structure for pattern recognition and path finding, in z
      int fStructDistZ;

      /// Distance for connecting hits through A*
      int fAStarConnectDist;
      /// Distance for associating hits with a path
      int fPathHitConnectDist;
      /// Distance for identifying extra hits by distance from an existing path
      int fExtraHitConnectDist;

      /// Maximum for merging overlapping path edge groups
      int fEdgeMergeStructDist;
      /// Maximum for merging overlapping vertex groups
      int fVertexMergeStructDist;
      /// Distance for finding vertices as points of divergence
      int fVertexFindStructDist;
      /// Distance for adding hits to vertices
      int fVertexHitStructDist;
      /// Distance for tagging a path as redundant when it comes close to a vertex
      int fVertexPathStructDist;
      /// Distance for merging two nearby clusters into a vetex
      int fClusterMergeStructDist;

      /// Maximum extent in MM pads for a path counted as an x-path
      int fXPathMaxPads;
      /// Minimum ratio between x end-start position and y-z end-start position for a path counted as an x-path
      double fXPathMinEndRatio;

      /// Minumum distance to run precheck in for vertex detection
      double fEdgePreDist;
      /// Threshold angle for precheck for vertex detection
      double fEdgePreAng;
      /// Minumum distance to collect hits in for vertex detection
      double fEdgeMinDist;
      /// Minumum distance to compute offset with for vertex detection
      double fEdgeMinDistLow;
      /// Fractional minimum hits at the edge for vertex detection, to avoid false positives from deltas
      double fEdgeMinHits;
      /// Fractional minimum asymmetry for vertex detection, to avoid false positives from straight tracks
      double fEdgeMaxSigma;
      /// Fractional range to collect hits in for vertex detection
      double fEdgeRange;
      /// Threshold to classify turn as a vertex
      double fEdgeThreshold;
      /// Offset from center to classify turn as a vertex
      double fEdgeOffset;

      /// Maximum distance for extrapolating hits when forming clusters
      int fHVClusterExtrapolateDist;
      /// Maximum total distance for extrapolating hits when forming clusters
      int fHVClusterExtrapolateLimit;
      /// Distance for connecting hits when forming clusters
      int fClusterConnectDist;
      /// Number of horizontal or vertial cells required in a row for categorising HV clusters
      int fMergeDist;
      /// Maximum number of isolated horizontal or vertical cells at the edge of a HV cluster
      int fHVEdgeDist;
      /// Distance in cells either side to sample when classifying HV clusters
      int fThresholdAngleRange;
      /// Angle for discriminating between horizontal and vertical clusters
      float fThresholdAngle;
      /// Minimum points along path required when establishing angle from dichotomy technique
      int fDichotomyCutoff;

      /// Minumum x-extent to break up a junction
      int fXSizeThreshold;
      /// Minumum path size to form from breaking up a junction
      int fPathSizeThreshold;
      /// Allow breaking long junctions between connecting paths as well as at edges facing away from them
      bool fBreakInMiddle;

      /// Distance to check when cleaning up anomalous hits near the junction
      int fAnomCheckDist;
      /// Distance to use for control when cleaning up anomalous hits near the junction
      int fAnomProjectDist;
      /// Threshold when cleaning up anomalous hits near the junction
      double fAnomMaxOffs;

      /// Maximum isolated clusters at the start of a path before they're merged into a nearby vertex
      int fHVClusterMaxIso;

    double fBField;

  };
}

#endif
