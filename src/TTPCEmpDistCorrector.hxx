#ifndef TTPCEmpDistCorrector_h
#define TTPCEmpDistCorrector_h

#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <utility>  //for pair
#include <vector>

///    \brief Provides services to the TPC fitting methods to correct for extra distortions
///    \author Christian Hansen, christian.hansen@cern.ch
///
///    TTPCEmpDistCorrector
///   
///    This class uses a 3D distortion map, achieved by using laser map and cosmics, 
///    to provide a service to correct for deviations of reconstructed TPC hits due to 
///    non perfect nominal fields and gas inhomogeneouses. This class should normally
///    not be called directly, but one can turn on/off Distortion Correction from 
///    TTPCClusterCorrection.


class TTPCEmpDistCorrector {
  
public:
  TTPCEmpDistCorrector();
  ~TTPCEmpDistCorrector();

  ///Need to set the complete file name for where the distortion map is 
  ///saved (since for now it is simply an ASCII file, should soon be in DB)
  bool SetDistMap(TString distMap);
  bool SetAllDistMap();
  bool DistMapIsSet() {return fDistMapIsSet;};

  ///Takes a 3D position as input. This should normaly be a reconstructed point.
  ///As an output it gives a new position for the reconstructed point corrected
  ///for e.g. nonperfect nominal fields.
  ///To be used in the "inverse approach", i.e. by the chi-square fit.
  bool GetDistortionCorrectedPoint(const TVector3& recPointWithHomogenousFields, TVector3 &corrPos);

  ///Takes a 3D position as input. This should normaly be a point on a helix
  ///defined by a set of track parameters suggested by Minuit.
  ///The other argument is a TLorentzVector and should be a point on the pad plane.
  ///The distortions belonging to the the track point will be applied on the pad plane point
  ///(i.e. only the pad plane point's y, z and t will be changed, x is always on pad plane).
  ///To be used by the "direct approach", i.e. by the maximum likelihood fit. 
  ///Since the likelihood reconstruction gets a very good seed track before the fit this only need to be 
  ///called once per each track point (not for every small change in the track parameters 
  ///during the optimization)???? OR????
  bool ApplyDistOnPadPoint(const TVector3& trackPoint, TLorentzVector &padPoint);
  
private: 
  bool fDistMapIsSet;
  bool fFourHits;
  double fMaxDriftInX;
  double fGlobalXOfPadPlane;
  std::vector<std::pair<TVector3, TLorentzVector > > fLaserHitsAndDistortions;
  void ReadInFromDB();
  bool GetClosestHit(const TVector3& recPos, std::pair<TVector3, TLorentzVector >& closestHit);
  bool GetClosestFourHits(const TVector3& recPos, std::vector<std::pair<TVector3,TLorentzVector> >& closestFourHits) ;
  bool Correct(const TVector3& recPos, const TVector3& fullDriftDistCorr, TVector3& corrPos);
  bool CorrectWeightedAverage(const TVector3& recPos,std::vector<std::pair<TVector3,TLorentzVector> >& closestFourHits , TVector3& corrPos);
  bool ApplyWeightedAverage(const TVector3& trackPos,std::vector<std::pair<TVector3,TLorentzVector> >& closestFourHits , TLorentzVector& padPlanePoint) ;
  bool Apply(const TVector3& trackPos, const TLorentzVector& fullDriftDist, TLorentzVector& padPlanePoint);
};
  
#endif // #ifdef TTPCEmpDistCorrector_hxx
