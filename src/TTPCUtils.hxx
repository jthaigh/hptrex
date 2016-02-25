#ifndef TTPCUtils_hxx_seen
#define TTPCUtils_hxx_seen

#include "TTPCPath.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCT0.hxx"
#include <TG4Trajectory.hxx>

namespace TTPCUtils {
  /// Sorting function order the Paths from longest to shortest.
  bool SortShortestToLongest( ND::THandle< ND::TTPCPath > first, ND::THandle< ND::TTPCPath > second );

  /// Extract the detector bit and convert into a string of type "TPCx".
  std::string DetectorToString( ND::THandle<ND::TReconBase> track );

  /// Act as the copy constructor that TReconVertex doesn't have.
  ND::THandle<ND::TReconVertex> CopyCreateTReconVertex(ND::TReconVertex *Original);

  /// Return the matching chi2 between the track state and the given cluster.
  /// Also provides the breakdown of the chi2 for the 3 directions.
  double State2CluChi2(State helixState, ND::THandle<ND::TTPCHVCluster> Cluster, double Chi2[3],int NDOF[3]);

  /// Return the matching chi2 between the track state and the given cluster.
  double State2CluChi2(State helixState, ND::THandle<ND::TTPCHVCluster> Cluster);

  bool SafeSort( double first, double second );

  double GetXWithT0(double oldX, double T0, ND::THandle<ND::THit> RefHit);

  State InvertDirection( State inputState );

  /// Find which ends of the tracks are the closest to one another.
  /// This is very important when matching and merging tracks.
  void FindClosestEnds(ND::THandle<ND::TReconBase> PathA, ND::THandle<ND::TReconBase> PathB, unsigned int &EndA, unsigned int &EndB);

  /// Simple code to decide crudely the sense based on the first two clusters
  int SenseFromTwoClusters(ND::THandle<ND::TTPCPath> Path2, ND::THandle<ND::TTPCHVCluster> ChosenClu);

  /// Merge 2 paths while taking care of the the cluster ordering
  ND::THandle<ND::TTPCPath> MergePaths(ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPath> PathB);

  /// Merge 2 paths using the Kalman filter fit. This should be used only for cathode crossers a prori.
  ND::THandle<ND::TReconBase> MergeAndFitObjectsWithRecPack(ND::THandle<ND::TReconBase> T1, ND::THandle<ND::TReconBase> T2);

  ND::THandle<ND::TG4Trajectory> FindTrueTrajectories( ND::THandle<ND::THitSelection> recoHits, double &complete, double &clean);

  bool TrueStateNear3Dpoint( ND::THandle<ND::THitSelection> recoHits, TVector3 pos3d, State &trueState);

  void PrintTrackInfo(ND::THandle<ND::TReconBase> track);

  void PrintConstituentMap(ND::THandle<ND::TReconBase> track, int level);

  int GetIdFromOAEvent(ND::THandle<ND::TReconBase> track);

  void ExtractPIDsFromPatterns(ND::THandle<ND::TReconObjectContainer> PatternContainer, ND::THandle<ND::TReconObjectContainer> PIDcontainer);

  void HVClustersPrintout(ND::THandle<ND::THitSelection> inputClu, bool Extended = false);
};


#endif
