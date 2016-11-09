#ifndef TTPCRecPackUtils_hxx_seen
#define TTPCRecPackUtils_hxx_seen

#include "TTPCHVCluster.hxx"
#include <TRecPackManager.hxx>

namespace TTPCRecPackUtils {

  /// Parameters that we want to locally modify in TREx
  struct SavedParam {
    int length_sign;
    bool allow_zero;
    bool eloss_on;
    bool ms_on;
    bool eloss_fluct_on;
    bool unique_surface;
  };

  void InitRecPackManager(void);
  
  /// Modify the configuration of the RecPack propagation to suit propagation inside the TPCs.
  /// Return values that we are going to modify to not mess up the rest of the reconstruction.
  SavedParam InitForPropagationInTPC();

  /// Propagate the given state to the given cluster.
  /// The T0 must have been provided to the cluster to get the right X position.
  /// The given state is modified to make it easy to propagate
  /// along the path, always propagating from one cluster to the next.
  bool PropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState, double &Length);

  /// Convenient method calling PropagateToHVCluster with its length argument.
  bool PropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState);

  /// Propagate the given state to the given cluster.
  /// The T0 must have been provided to the cluster to get the right X position.
  /// The given state is modified to make it easy to propagate
  /// along the path, always propagating from one cluster to the next.
  bool FullPropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState, double &Length);

  /// Convenient method calling PropagateToHVCluster with its length argument.
  bool FullPropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState);

  /// Reset the configuration of the RecPack propagation to what
  /// it was before calling InitForPropagationInTPC().
  void ResetAfterPropagationInTPC(SavedParam &Params);
};



#endif

