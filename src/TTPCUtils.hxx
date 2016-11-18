#ifndef TTPCUtils_hxx_seen
#define TTPCUtils_hxx_seen

#include "TTRExPath.hxx"
#include "TTRExHVCluster.hxx"

namespace TTPCUtils {

  bool Curvature_to_MomentumAndCharge(const TVector3& pos, const TVector3& dir, double curv, double& p, double& q);

  bool MomentumAndCharge_to_Curvature(const TVector3& pos, const TVector3& dir, double p, double q, double& curv);

  bool SafeSort( double first, double second );

  /// Find which ends of the tracks are the closest to one another.
  /// This is very important when matching and merging tracks.
  void FindClosestEnds(trex::TTRExPath& PathA, trex::TTRExPath& PathB, unsigned int &EndA, unsigned int &EndB);

  /// Simple code to decide crudely the sense based on the first two clusters
  int SenseFromTwoClusters(trex::TTRExPath& Path2, trex::TTRExHVCluster& ChosenClu);

  /// Merge 2 paths while taking care of the the cluster ordering
  trex::TTRExPath* MergePaths(trex::TTRExPath& PathA, trex::TTRExPath& PathB);

  void ReverseStateSenseAndCharge(std::vector<double>& propagState);

};


#endif
