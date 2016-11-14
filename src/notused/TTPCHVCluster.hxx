#ifndef TTPCHVCluster_hxx_seen
#define TTPCHVCluster_hxx_seen

#include <THit.hxx>
#include <TComboHit.hxx>
#include <recpack/Surface.h>

namespace ND {
  class TTPCHVCluster;
}


/// Represents a horizontal or vertical cluster of pads.
/// This is the object needed to calculate the likelihood.
class ND::TTPCHVCluster : public THit {
public:
  TTPCHVCluster();
  TTPCHVCluster(ND::THitSelection& hits, bool isAVerticalCluster);
  virtual ~TTPCHVCluster();

  //////// A lot of code is copy-pasted from TComboHit because
  //////// TComboHit had too many private variables to make inheritance
  //////// practical
  /// Get the geometry id for the volume that contains this hit.  This
  /// should be handled carefully with a TComboHit since the volume
  /// containing this hit might not be enumerated.
  ND::TGeometryId GetGeomId(void) const;

  /// Return the calibrated "charge" for the hit.
  double GetCharge(void) const;
  double GetMaxCharge(void) const;
  double GetIntCharge(void) const; // Integrated

  /// Return the calibrated "time" for the hit.
  double GetTime(void) const;

  /// Return the largest difference in drift distance between the waveforms.
  double GetDeltaDrift(void) const;

  /// The center of the volume associated with this hit.
  const TVector3& GetPosition(void) const;

  /// The calibrated center of the volume associated with this hit
  /// which means taking into account the T0 and the field corrections.
  const TVector3& GetCalibPosition(void) const;

  /// Return true if this hit has useful X information.
  bool IsXHit(void) const;

  /// Return true if this hit has useful Y information.
  bool IsYHit(void) const;

  /// Return true if this hit has useful Z information.
  bool IsZHit(void) const;

  /// Return the "spread" of the hit position.  This is the extent in the X,
  /// Y, and Z directions.  For instance, a P0D bar is about
  /// (3cm)x(2m)x(1.5cm), so the spread is (1.3cm)x(1m)x(0.75cm)
  const TVector3& GetSpread(void) const;

  /// Return the "uncertainty" of the hit position.  
  const TVector3& GetUncertainty(void) const;

  /// Return the "uncertainty" for the time measurement.
  double GetTimeUncertainty(void) const;

  /// Get the parent geometry object.
  TGeoNode* GetParentNode(int i) const;

  /// Get a constant version of hits associated with this combination hit.
  const THitSelection& GetHits(void) const;

  /// Add a hit to the combo hits.
  void AddHit(ND::THandle<ND::THit>& hit);

  /// Copy a hit selection into the combo hit.
  void AddHitSelection(THitSelection& hits);

  /// Allow access to the non-constant version of the hits.
  void OpenHits();

  /// Recalculate internal state variables in this hit list of
  /// associated hits.
  void CloseHits() const;

  /// Define the cluster as horizontal
  void SetIsHorizontal();
  /// Define the cluster as vertical
  void SetIsVertical();

  /// Returns true if this is a horizontal cluster
  bool IsHorizontal();
  /// Returns true if this is a vertical cluster
  bool IsVertical();

  /// Define the cluster as an end node
  void SetEndNode();
  /// Clear the EndNode status
  void ClearEndNode();

  /// Returns true if this cluster is an end node
  bool IsEndNode();

  /// Get the field correction for this measurement. This is normally computed based on the seeding
  double GetFieldCorrection(void);
  double GetFieldCorrection(void) const;
  /// Add the field correction of this measurement. This is normally computed based on the seeding
  void SetFieldCorrection(double val);

  /// Is this pad against the vertical edge of the MM ?
  bool IsAtVertEdge();
  /// Is this pad against the horizontal edge of the MM ?
  bool IsAtHoriEdge();

  /// Does any pad have suspicious timing sucha as large difference between the peak time
  /// and the fitted time ?
  bool HasSuspiciousPadTiming();

  /// Returns number of peaks in the waveforms with the maximum number of peaks in the cluster.
  unsigned int GetMaxNPeaks();
  /// Returns number of saturated waveforms in the cluster.
  unsigned int GetNSaturated();


  int AllFitted();

  /// Short cut to access the X position of the cluster
  double X();
  /// Short cut to access the Y position of the cluster
  double Y();
  /// Short cut to access the Z position of the cluster
  double Z();
  /// Returns the "calibrated" X position, i.e. takes into account T0
  double CalibX();
  /// Returns the "calibrated" Y position, i.e. adds DeltaY correction
  double CalibY();
  /// Returns the "calibrated" Z position, i.e. adds DeltaY correction
  double CalibZ();
  /// Provide the calibrated position using the provided T0 and DeltaY and DeltaZ
  TVector3 GetCalibratedPosition();
  /// Returns the drift distance for this cluster.
  double GetDriftDistance();
  /// Store the T0 for this cluster and calculate the "calibrated" X position
  void SetT0(double T0);
  /// Store the DeltaY and DeltaZ from field corrections for example
  /// and calculate the "calibrated" Y and Z positions.
  void SetDeltaYZ(double dY, double dZ);
  /// Is the T0 available ?
  bool hasT0();
  /// Are DeltaY and DeltaZ corrections available ?
  bool hasDeltaYZ();
  /// Access DeltaY directly
  double GetDeltaY();
  /// Access DeltaZ directly
  double GetDeltaZ();

  /// Defines if this cluster should be used for the seeding based
  /// on the nature of the cluster and its content.
  void SetOkForSeed(bool OkOrNot) {fOkForSeed = OkOrNot;};
  /// Reports is this cluster should be used for the seeding based
  /// on the nature of the cluster and its content.
  bool isOkForSeed() {return fOkForSeed;};

  /// Defines if this cluster should be used for the likelihood fit based
  /// on the nature of the cluster and its content.
  void SetOkForFit(bool OkOrNot) {fOkForFit = OkOrNot;};
  /// Reports is this cluster should be used for the likelihood fit based
  /// on the nature of the cluster and its content.
  bool isOkForFit() {return fOkForFit;};

  /// Defines if this cluster is turned off or on.
  /// This is separate from the isOkForSeed and isOkForFit which check
  /// if the cluster is of good quality for an algorithm.
  /// The "usability" is about if we want the algorithms to use it or not,
  /// for example in the case of spiralling tracks where we want to fit only
  /// the first loop, so we turn off the other clusters.
  void SetUsable(bool UsableOrNot) {fUsable = UsableOrNot;};
  /// Reports is this cluster has been turned off or on.
  bool isUsable() {return fUsable;};

  /// Define and store the propagation surface for this cluster
  /// in order to speed up propagations.
  void PreparePropagSurf(const std::string surfname);
  /// Return the propagation surface to perform a propagation.
  Surface& GetPropagSurf();
  /// Call this once done with all propagations to clean up
  /// the RecPack geometry from useless surfaces.
  void RemovePropagSurf();
void PrintSurface() {std::cout<<fPropagSurfName<<std::endl;};

  /// Convert into a pure oaEvent for compatibility with the rest of the ND280 software.
  ND::THandle<ND::TComboHit> ConvertToOAEvent();


private:
  /// Runtime parameters
  bool fUseExtrapolatedCharge;      //! Don't save in ROOT.
  bool fUseFittedCharge;            //! Don't save in ROOT.
  bool fUseExtrapolatedIntCharge;   //! Don't save in ROOT.
  double fHeightToAreaFactor;       //! Don't save in ROOT.

  /// Mark the cache as invalid.
  mutable bool fHitsOpen; //! Don't save in ROOT.

  /// The hits that contribute to this combined hit.
  THitSelection fHits;

  /// The current volume and parents that contain this cluster.  This has
  /// the constraint that it must contain at least two nodes, the volume
  /// itself, and the parent volume.  If the size is less than 2, the vector
  /// is considered empty.
  mutable std::vector<TGeoNode*> fParentNodes; //! Don't save in ROOT.

  /// The geometry id of the vlume that contains the hit.  This is the
  /// volume of the averaged position, and may not contain any of the
  /// constituient hits.
  mutable TGeometryId fGeomId; //! Don't save in ROOT.

  /// The measured "charge" for this hit.
  mutable double fCharge; //! Don't save in ROOT.

  mutable double fIntCharge; //! Don't save in ROOT.

  /// The measured "time" for this hit.
  mutable double fTime; //! Don't save in ROOT.

  /// The largest difference in difference between the waveforms.
  mutable double fDeltaDrift; //! Don't save in ROOT.

  /// The position of the hit.
  mutable TVector3 fPosition; //! Don't save in ROOT.

  /// Drift distance calculated from the given T0.
  mutable double fDriftDistance; //! Don't save in ROOT.
  /// The position of the hit taking into account T0 and E/F field corrections to Y and Z.
  mutable TVector3 fCalibratedPosition; //! Don't save in ROOT.
  mutable double fT0; //! Don't save in ROOT.
  mutable double fDeltaYZ[2]; //! Don't save in ROOT.

  /// The spread of the hit position.
  mutable TVector3 fSpread; //! Don't save in ROOT.

  /// The uncertainty of the hit position.
  mutable TVector3 fUncertainty; //! Don't save in ROOT.

  /// The uncertainty in the hit time.
  mutable double fTimeUncertainty; //! Don't save in ROOT.

  /// Flags for what type of position information is available.
  mutable bool fIsXHit; //! Don't save in ROOT.
  mutable bool fIsYHit; //! Don't save in ROOT.
  mutable bool fIsZHit; //! Don't save in ROOT.

  /// True if this is a vertical cluster, false if this is horizontal
  bool fIsVertical;

  /// True if a node should be created for this cluster 
  bool fEndNode;

  /// Field distortion correction in Y or Z, whichever is a free coordinate.
  /// This is the correction that has to be added to the predicted coordinate
  /// of the track at the center of the corresponding pad.
  mutable double fFieldCorrection; 

  /// Maximum charge in a single time sample in the cluster.
  mutable double fMaxCharge;
  /// Max number of peaks per waveform for all the pads.
  mutable int fMaxNPeaks;
  /// Number of saturated waveforms
  mutable int fNSaturated;

  /// Is the cluster considered too close to the edge ?
  mutable bool fAtVertEdge;
  mutable bool fAtHoriEdge;

  mutable bool fSuspectPadTiming;

  /// Are all the waveforms were fitted ?
  mutable bool fAllFitted;

  /// Should this cluster be used for seeding ?
  bool fOkForSeed;
  /// Should this cluster be used for fitting ?
  bool fOkForFit;
  /// Define which part of the track should be used.
  /// This is particular for spiraling delta-rays.
  bool fUsable;

  /// RecPack propagation surface.
  Surface fPropagSurf;
  std::string fPropagSurfName;

protected:

  ClassDef(TTPCHVCluster,1);
};


#endif
