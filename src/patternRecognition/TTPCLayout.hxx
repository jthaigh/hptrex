#ifndef TTPCLAYOUT_HXX
#define TTPCLAYOUT_HXX

// c++
#include <iostream>
#include <map>
#include <algorithm>
#include <sstream>

// ROOT
#include <TVector3.h>

// nd280
#include <TOARuntimeParameters.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <THit.hxx>
#include <TGeomInfo.hxx>
#include "TTPCCalibration.hxx"
#include <TTPCGeom.hxx>
#include <TGeometryId.hxx>
#include "TTPCCalibration.hxx"

// TRex
#include <TTPCHitPad.hxx>

namespace ND{

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
  /// Structure for 2D cell with x and y id
  struct TTPCCell2D{
    int x;
    int y;
  };
  /// Structure for 3D cell with x, y and z id
  struct TTPCCell3D{
    int x;
    int y;
    int z;
  };
  /// Structure for 3D cell with x, y and z id and edge status
  struct TTPCCellInfo3D{
    int x;
    int y;
    int z;
    int edgeX;
    int edgeY;
    int edgeZ;
    int segX;
    int segY;
    int segZ;
  };
  /// Structure containing minima, maxima and ranges in x, y and z
  struct TTPCCellRanges3D{
    int minX;
    int minY;
    int minZ;

    int maxX;
    int maxY;
    int maxZ;

    int sizeX;
    int sizeY;
    int sizeZ;
  };
  /// Structure containing pad information
  struct TTPCPadStruct{
    int tpc;
    int half;
    int mm;
    int pad;
    int rawpad;
  };

  /// Class containing relevant information and methods for expressing 3D layout of hits in a way that makes sense to feature finding and path finding algorithms
  class TTPCLayout{
    public:
      /// Default constructor
      TTPCLayout();
      /// Default destructor
      ~TTPCLayout();

      // Set min and max times either side of the cathode
      void SetTimeRanges(double tNMin, double tNMax, double tPMin, double tPMax);
      /// Set cell id ranges in x, y and z
      void SetRanges(int minX,int maxX, int minY,int maxY, int minZ,int maxZ);
      /// Get cell id ranges in x, y and z
      void GetRanges(int& sizeX,int& minX,int& maxX, int& sizeY,int& minY,int& maxY, int& minZ,int& maxZ,int& sizeZ);
      /// Get 2D cell ranges in x view, y view or z view, at axis of 1, 2 or 3 respectively
      void GetRanges(int& sizeX,int& minX,int& maxX, int& sizeY,int& minY,int& maxY, int axis=1);

      /// Get cell id in x, y and z based on a hit position
      ND::TTPCCellInfo3D GetPadPosID(ND::THandle<ND::THit> hit, int tpcMask=-1);
      /// Get cell id in x, y and z based on a 3D position
      ND::TTPCCellInfo3D GetPadPosID(TVector3 pos, double time, int tpcMask=-1);

      /// Convert cell id in x, y and z to a single unique long integer
      long Mash(int x, int y, int z);
      /// Convert cell id in y and z to a single unique long integer
      long MashYZ(int y, int z);
      /// Convert cell id in x, y and z to a single unique long integer, returning -1 if bad ids are supplied
      long SafeMash(int x, int y, int z);
      /// Convert a single unique long integer into a cell id in x, y and z
      ND::TTPCCell3D UnMash(long id);

      /// Get converts a 3D position into MM pad information
      ND::TTPCPadStruct GlobalXYZToPos(TVector3 pos);

      /// Find appropriate distances for a given type of connection
      void GetTypeDistances(int& distX, int& distY, int& distZ, ND::TTPCConnection::Type type);

      // getters
      /// Get use of experimental edge detection for helping with deltas
      bool GetUseAltEdgeDetection(){ return fUseAltEdgeDetection; }
      /// Get use of experimental hit association to help grab unused hits
      int GetUseAltHitAssociation(){ return fUseAltHitAssociation; }
      /// Get use of charge cut for path finding
      int GetUsePatRecPathologyCut(){ return fUsePatRecPathologyCut; }

      /// Get absolute cut for minimum charge on cells going into pattern recognition
      double GetChargeCut(){ return fChargeCut; }
      /// Get cut on charge ratio for negative peak arriving before the main one
      double GetEarlyNegativePeakCut(){ return fEarlyNegativePeakCut; }
      /// Get cut on charge ratio for negative peak arriving after the main one
      double GetLateNegativePeakCut(){ return fLateNegativePeakCut; }
      /// Get absolute cut for maximum number of saturated hits in an ASIC before it's tagged
      double GetASICSaturationCut(){ return fASICSaturationCut; }
      /// Get absolute cut for maximum occupancy of a sub-region of an ASIC before its hits are tagged as dodgy
      double GetASICSubOccupancyCut(){ return fASICSubOccupancyCut; }
      /// Get absolute cut for maximum occupancy of an ASIC before its hits are tagged as dodgy
      double GetASICOccupancyCut(){ return fASICOccupancyCut; }
      /// Get distance for expanding from a dodgy ASIC to tag dodgy hits
      double GetASICSatExpansion(){ return fASICSatExpansion; }
      /// Get distance for expanding from a dodgy ASIC to tag dodgy hits
      double GetASICOccExpansion(){ return fASICOccExpansion; }
      /// Get empirical number of y pads in sub region of ASIC with overflowing charge
      double GetASICSplittingY(){ return fASICSplittingY; }
      /// Get empirical number of z pads in sub region of ASIC with overflowing charge
      double GetASICSplittingZ(){ return fASICSplittingZ; }

      /// Get drift speed
      double GetDriftSpeed(){ return fDriftSpeed; }

      /// Get minimum hit time in negative half
      double GetTNMin(){ return fTNMin; }
      /// Get maximum hit time in negative half
      double GetTNMax(){ return fTNMax; }
      /// Get minimum hit time in positive half
      double GetTPMin(){ return fTPMin; }
      /// Get maximum hit time in positive half
      double GetTPMax(){ return fTPMax; }
      /// Get number of bins in negative half
      int GetNegativeBins(){ return fTNegativeBins; }
      /// Get number of bins in positibe half
      int GetPositiveBins(){ return fTPositiveBins; }
      /// Get total number of time bins
      int GetTBins(){ return fTBins; }

      /// Get gap between MM pads
      double GetPadGap(){ return fPadGap; }
      /// Get distance between MM pad centres in y direction
      double GetPadPitchY(){ return fPadPitchY; }
      /// Get distance between MM pad centres in z direction
      double GetPadPitchZ(){ return fPadPitchZ; }

      /// Get wheter a track in this event crosses the central cathode
      bool GetXCathodeCross(){ return fXCathodeCross; }

      /// Get size of an individual x cell
      double GetXCellSize(){ return fXCellSize; }
      /// Get time bin size
      double GetTWidth(){ return fTWidth; }

      /// Get number of MM pads
      int GetMMPads(){ return fMMPads; }
      /// Get number of MM pads in y direction
      int GetYPads(){ return fYPads; }
      /// Get number of MM pads in z direction
      int GetZPads(){ return fZPads; }

      /// Get increase in x id Getrom crossing central cathode
      int GetXShiftFromC(){ return fXShiftFromC; }
      /// Get increase in y id Getrom crossing between two MM columns
      int GetYShiftFromXZ(){ return fYShiftFromXZ; }
      /// Get increase in y id Getrom moving to the next MM volume in y direction
      int GetYShiftFromMM(){ return fYShiftFromMM; }
      /// Get increase in z id Getrom moving to the next MM volume in z direction
      int GetZShiftFromMM(){ return fZShiftFromMM; }
      /// Get increase in z id Getrom moving to the next TPC
      int GetZShiftFromTPC(){ return fZShiftFromTPC; }

      /// Get minimum x id in event
      int GetMinX(){ return fMinX; }
      /// Get minimum y id in event
      int GetMinY(){ return fMinY; }
      /// Get minimum z id in event
      int GetMinZ(){ return fMinZ; }
      /// Get maximum x id in event
      int GetMaxX(){ return fMaxX; }
      /// Get maximum y id in event
      int GetMaxY(){ return fMaxY; }
      /// Get maximum z id in event
      int GetMaxZ(){ return fMaxZ; }
      /// Get number of x ids in event
      int GetSizeX(){ return fSizeX; }
      /// Get number of y ids in event
      int GetSizeY(){ return fSizeY; }
      /// Get number of z ids in event
      int GetSizeZ(){ return fSizeZ; }

      /// Get whether MM gaps can be crossed in x
      bool GetJumpX(){ return fJumpX; }
      /// Get whether MM gaps can be crossed in y
      bool GetJumpY(){ return fJumpY; }
      /// Get whether MM gaps can be crossed in z
      bool GetJumpZ(){ return fJumpZ; }

      /// Get offset added to x to effect a gap from crossing central cathode
      int GetGapOffsetX(){ return fGapOffsetX; }
      /// Get offset added to y to effect a gap from crossing between MM volumes in y direction 
      int GetGapOffsetY(){ return fGapOffsetY; }
      /// Get offset added to z to effect a gap from crossing between MM volumes in z direction
      int GetGapOffsetZ(){ return fGapOffsetZ; }
      /// Get multiplier for extra offsets applied to adjacent directions at a gap
      int GetGapOffsetAdjacent(){ return fGapOffsetAdjacent; }

      /// Get minimum number of pads constituting a useful pattern
      int GetMinPatternPads(){ return fMinPatternPads; }
      /// Get minimum number of clusters constituting a useful path
      int GetMinPathClusters(){ return fMinPathClusters; }

      /// Get rate at which delta search radius spreads out
      double GetDeltaSpreadRate(){ return fDeltaSpreadRate; }

      /// Get whether to use hits not on direct edge of an MM volume when looking for edges
      bool GetUseIndirectEdges(){ return fUseIndirectEdges; }
      /// Get number of layers at edge of sub event to search for track ends
      int GetEdgeLayers(){ return fEdgeLayers; }

      /// Default A* scale for one cell apart connection in x direction
      float fAStarXScale;
      /// Default A* scale for one cell apart connection in y direction
      float fAStarYScale;
      /// Default A* scale for one cell apart connection in z direction
      float fAStarZScale;
      /// Factor to weight A* heuristic by, to alter the performance and effectiveness of the algorithm
      float fAStarHeuristicFactor;
      /// Multiplicative penalty term to scale costs of connections breaking charge cut by
      float fAStarPathologyPenalty;
      /// Multiplicative penalty term to scale costs of connections breaking charge cut by, for associating hits when looking for edges only
      float fAStarAssociatePathologyPenalty;

      /// Default A* scale for one cell apart connection in x direction
      float GetAStarXScale(){ return fAStarXScale; }
      /// Default A* scale for one cell apart connection in y direction
      float GetAStarYScale(){ return fAStarYScale; }
      /// Default A* scale for one cell apart connection in z direction
      float GetAStarZScale(){ return fAStarZScale; }
      /// Factor to weight A* heuristic by, to alter the performance and effectiveness of the algorithm
      float GetAStarHeuristicFactor(){ return fAStarHeuristicFactor; }
      /// Get multiplicative penalty term to scale costs of connections breaking charge cut by
      float GetAStarPathologyPenalty(){ return fAStarPathologyPenalty; }
      /// Get multiplicative penalty term to scale costs of connections breaking charge cut by, for associating hits when looking for edges only
      float GetAStarAssociatePathologyPenalty(){ return fAStarAssociatePathologyPenalty; }

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

      /// Get A* cost for associating hits with a path
      float GetAltPathHitConnectDist(){ return fAltPathHitConnectDist; }
      /// Get A* cost for identifying extra hits by distance from an existing path
      float GetAltExtraHitConnectDist(){ return fAltExtraHitConnectDist; }
      /// Get A* cost for removing redundant edges by distance from an existing path
      float GetAltEdgeHitConnectDist(){ return fAltEdgeHitConnectDist; }
      /// Get A* cost for removing redundant vertices by distance from an existing vertex
      float GetAltVertexHitConnectDist(){ return fAltVertexHitConnectDist; }

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
      /// Get number of horizontal or vertial cells required in a row for categorising HV clusters
      int GetMergeDist(){ return fMergeDist; }
      /// Get maximum number of isolated horizontal or vertical cells at the edge of a HV cluster
      int GetHVEdgeDist(){ return fHVEdgeDist; }
      /// Get distance in cells either side to sample when classifying HV clusters
      int GetThresholdAngleRange(){ return fThresholdAngleRange; }
      /// Get angle for discriminating between horizontal and vertical clusters
      float GetThresholdAngle(){ return fThresholdAngle; }
      /// Get minimum points along path segment required when establishing angle from dichotomy technique
      int GetDichotomyCutoff(){ return fDichotomyCutoff; }

      /// Get minumum x-extent to break up a junction
      int GetXSizeThreshold(){ return fXSizeThreshold; }
      /// Get minumum path size to form from breaking up a junction
      int GetPathSizeThreshold(){ return fPathSizeThreshold; }
      /// Get whether to allow breaking long junctions between connecting paths as well as at edges facing away from them
      bool GetBreakInMiddle(){ return fBreakInMiddle; }

      /// Get distance to check when cleaning up anomalous hits near the junction
      int GetAnomCheckDist(){ return fAnomCheckDist; }
      /// Get distance to use for control when cleaning up anomalous hits near the junction
      int GetAnomProjectDist(){ return fAnomProjectDist; }
      /// Get threshold when cleaning up anomalous hits near the junction
      double GetAnomMaxOffs(){ return fAnomMaxOffs; }

      /// Get maximum isolated clusters at the start of a path before they're merged into a nearby vertex
      int GetHVClusterMaxIso(){ return fHVClusterMaxIso; }

      /// Get minimum fraction of non-detla hits for avoiding delta classification
      float GetNonDeltaFraction(){ return fNonDeltaFraction; }
      /// Get minimum absolute non-detla hits for avoiding delta classification
      int GetNonDelta(){ return fNonDelta; }

    private:
      /// Use experimental edge detection for helping with deltas 
      bool fUseAltEdgeDetection;
      /// Use experimental hit association to help grab unused hits
      int fUseAltHitAssociation;
      /// Use charge cut for path finding
      int fUsePatRecPathologyCut;

      /// Absolute cut for minimum charge on cells going into pattern recognition
      double fChargeCut;
      /// Cut on charge ratio for negative peak arriving before the main one
      double fEarlyNegativePeakCut;
      /// Cut on charge ratio for negative peak arriving after the main one
      double fLateNegativePeakCut;
      /// Absolute cut for maximum number of saturated hits in an ASIC before it's tagged
      double fASICSaturationCut;
      /// Absolute cut for maximum occupancy of a sub-region of an ASIC before its hits are tagged as dodgy
      double fASICSubOccupancyCut;
      /// Absolute cut for maximum occupancy of an ASIC before its hits are tagged as dodgy
      double fASICOccupancyCut;
      /// Distance for expanding from a dodgy ASIC to tag dodgy hits
      double fASICSatExpansion;
      /// Distance for expanding from a dodgy ASIC to tag dodgy hits
      double fASICOccExpansion;
      /// Empirical number of y pads in sub region of ASIC with overflowing charge
      double fASICSplittingY;
      /// Empirical number of z pads in sub region of ASIC with overflowing charge
      double fASICSplittingZ;

      /// Drift speed
      double fDriftSpeed;

      /// Minimum hit time in negative half
      double fTNMin;
      /// Maximum hit time in negative half
      double fTNMax;
      /// Minimum hit time in positive half
      double fTPMin;
      /// Maximum hit time in positive half
      double fTPMax;
      /// Number of bins in negative half
      int fTNegativeBins;
      /// Number of bins in positibe half
      int fTPositiveBins;
      /// Total number of time bins
      int fTBins;

      /// Gap between MM pads
      double fPadGap;
      /// Distance between MM pad centres in y direction
      double fPadPitchY;
      /// Distance between MM pad centres in z direction
      double fPadPitchZ;

      /// Wheter a track in this event crosses the central cathode
      bool fXCathodeCross;
      /// Size of an individual x cell
      double fXCellSize;
      /// Time bin size
      double fTWidth;

      /// Number of MM pads
      int fMMPads;
      /// Number of MM pads in y direction
      int fYPads;
      /// Number of MM pads in z direction
      int fZPads;

      /// Increase in x id from crossing central cathode
      int fXShiftFromC;
      /// Increase in y id from crossing between two MM columns
      int fYShiftFromXZ;
      /// Increase in y id from moving to the next MM volume in y direction
      int fYShiftFromMM;
      /// Increase in z id from moving to the next MM volume in z direction
      int fZShiftFromMM;
      /// Increase in z id from moving to the next TPC
      int fZShiftFromTPC;

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

      /// Whether MM gaps can be crossed in x
      bool fJumpX;
      /// Whether MM gaps can be crossed in y
      bool fJumpY;
      /// Whether MM gaps can be crossed in z
      bool fJumpZ;

      /// Offset added to x to effect a gap from crossing central cathode
      int fGapOffsetX;
      /// Offset added to y to effect a gap from crossing between MM volumes in y direction 
      int fGapOffsetY;
      /// Offset added to z to effect a gap from crossing between MM volumes in z direction
      int fGapOffsetZ;
      /// Multiplier for extra offsets applied to adjacent directions at a gap
      int fGapOffsetAdjacent;

      /// Minimum useful size of pattern
      int fMinPatternPads;
      /// Minimum useful size of path
      int fMinPathClusters;

      /// Rate at which delta search radius spreads out
      double fDeltaSpreadRate;

      /// Whether to use hits not on direct edge of an MM volume when looking for edges
      bool fUseIndirectEdges;
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

      /// A* cost for associating hits with a path
      float fAltPathHitConnectDist;
      /// A* cost for identifying extra hits by distance from an existing path
      float fAltExtraHitConnectDist;
      /// A* cost for removing redundant edges by distance from an existing path
      float fAltEdgeHitConnectDist;
      /// A* cost for removing redundant vertices by distance from an existing vertex
      float fAltVertexHitConnectDist;

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

      /// Minimum fraction of non-detla hits for avoiding delta classification
      float fNonDeltaFraction;
      /// Minimum absolute non-detla hits for avoiding delta classification
      int fNonDelta;
  };
}

#endif
