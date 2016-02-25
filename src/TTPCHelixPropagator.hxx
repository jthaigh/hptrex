#ifndef TTPCHelixPropagator_hxx_seen
#define TTPCHelixPropagator_hxx_seen

#include <TND280Event.hxx>
#include <TRecPackManager.hxx>

#include "TTPCHVCluster.hxx"

namespace ND {
  class TTPCHelixPropagator;
}

class ND::TTPCHelixPropagator {
public:

  virtual ~TTPCHelixPropagator() {};

  /// Get a pointer to the singleton instance of the helix propagator.
  static TTPCHelixPropagator& Get(void);

  /// Just set all the internal variables to 0.0.
  /// It is good practice to call this when you are
  /// done propagating on the given helix.
  void Reset();


  bool InitHelixPosDirQoP(State &RPState, bool FirstCluIsVertical);
  bool InitHelixPosDirQoP(double *Param, bool FirstCluIsVertical);

  void ReloadHelixPosTanCurv(double *Param);
  bool ReloadHelixPosDirQoP(State &RPState);
  bool ReloadHelixPosDirQoP(double *Param);

  /// Just return the quadrant where first cluster is.
  int GetQuadrant();
  /// Just return the sense of the track.
  int GetSense();

  /*!
    Propagate the given state to the given cluster.
    The T0 must have been provided to the cluster to get the right X position.
    This propagation works only in steps, going from cluster to cluster
    and saving the last result for the next propagation.
    DO NOT propagate for example from the first to the last cluster because
    depending on the topology of the track, you may get the wrong position
    on the helix.
  */
  bool PropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster);

  /// Return the results in position, direction (3D vector) and curvature.
  void GetHelixPosDirCurv(double *Result);

  /// Return the results in position, tangent x, tangent y or z (depends on the quadrant) and curvature.
  void GetHelixPosTanCurv(double *Result);

  /// Return the results in position, direction (3D vector) and charge over momentum.
  void GetHelixPosDirQoP(double *Result);

  /// Convert vector and covariance from PosDirCurv representation to PosDirQoP.
  void PosTanCurvToPosDirQoP(EVector &ptcVect, EMatrix &ptcCova, EVector &pdqpVect, EMatrix &pdqpCova );


private:
  TTPCHelixPropagator();

  static TTPCHelixPropagator* _helixPropagator;
  /// Calculates fPhiZY and fPhiQuad for the first point.
  void CalculatePhisAndCenters();
  /// Determine the quadrant at the first point.
  void FindQuadrant();

  double fX0;
  double fY0;
  double fZ0;
  double fDirX0;
  double fDirY0;
  double fDirZ0;
  double fPhiZY;
  double fPhiQuad;
  double fRho0;
  // The center of the circular projection in YZ.
  double fYc;
  double fZc;
  /*!
    Defines the quadrants
    
      \   2  /
       \    /
        \  /
      3  \/  1
         /\
        /  \
       /    \
      /   4  \
  
    Quadrants 1 and 3 are used when propagating horizontal clusters.
    Quadrants 2 and 4 are used when propagating vertical clusters.
    The phi angle is local to the quadrant which means that
    it is defined to always be well within the interval ]-pi/2;pi/2[.
    This way we always avoid problems when calculating
    for example DeltaPhi when phi approaches 2pi.
  */
  unsigned int fQuadrant;
  /// Keep the orientation of the first cluster.
  /// This is essential to the quadrant determination.
  unsigned int fFirstCluIsVert;

  /// This is equivalent to the sense in the RecPack states.
  /// Keep the sense from the Init call. Crucial for ReloadHelixPosTanCurv.
  int fInitSense;
};


namespace ND{
  TTPCHelixPropagator& helixPropagator();
};

#endif
