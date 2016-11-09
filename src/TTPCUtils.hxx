#ifndef TTPCUtils_hxx_seen
#define TTPCUtils_hxx_seen

#include "TTRExPath.hxx"
#include "TTPCHVCluster.hxx"

namespace TTPCUtils {

  bool SafeSort( double first, double second );

  /// Find which ends of the tracks are the closest to one another.
  /// This is very important when matching and merging tracks.
  void FindClosestEnds(ND::THandle<ND::TReconBase> PathA, ND::THandle<ND::TReconBase> PathB, unsigned int &EndA, unsigned int &EndB);

  /// Simple code to decide crudely the sense based on the first two clusters
  int SenseFromTwoClusters(ND::THandle<ND::TTPCPath> Path2, ND::THandle<ND::TTPCHVCluster> ChosenClu);

  /// Merge 2 paths while taking care of the the cluster ordering
  ND::THandle<ND::TTPCPath> MergePaths(ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPath> PathB);

};


#endif
