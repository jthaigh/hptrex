#ifndef TTPCUtils_hxx_seen
#define TTPCUtils_hxx_seen

#include "TTPCPath.hxx"
#include "TTPCHVCluster.hxx"

namespace TTPCUtils {

  /// Act as the copy constructor that TReconVertex doesn't have.
  ND::THandle<ND::TReconVertex> CopyCreateTReconVertex(ND::TReconVertex *Original);

  /// Return the matching chi2 between the track state and the given cluster.
  /// Also provides the breakdown of the chi2 for the 3 directions.
  double State2CluChi2(State helixState, ND::THandle<ND::TTPCHVCluster> Cluster, double Chi2[3],int NDOF[3]);

  /// Return the matching chi2 between the track state and the given cluster.
  double State2CluChi2(State helixState, ND::THandle<ND::TTPCHVCluster> Cluster);

  bool SafeSort( double first, double second );

  /// Find which ends of the tracks are the closest to one another.
  /// This is very important when matching and merging tracks.
  void FindClosestEnds(ND::THandle<ND::TReconBase> PathA, ND::THandle<ND::TReconBase> PathB, unsigned int &EndA, unsigned int &EndB);

  /// Simple code to decide crudely the sense based on the first two clusters
  int SenseFromTwoClusters(ND::THandle<ND::TTPCPath> Path2, ND::THandle<ND::TTPCHVCluster> ChosenClu);

  /// Merge 2 paths while taking care of the the cluster ordering
  ND::THandle<ND::TTPCPath> MergePaths(ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPath> PathB);

  /// Merge 2 paths using the Kalman filter fit. This should be used only for cathode crossers a prori.
  ND::THandle<ND::TReconBase> MergeAndFitObjectsWithRecPack(ND::THandle<ND::TReconBase> T1, ND::THandle<ND::TReconBase> T2);
};


#endif
