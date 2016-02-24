#ifndef TTPCFieldCorrector_h
#define TTPCFieldCorrector_h

#include "TVector3.h"
#include "TLorentzVector.h"
#include <TDriftManager.hxx>

///    \brief Provides services to the TPC fitting methods to correct for deviations due to inhomogenous fields
///    \author Christian Hansen, christian.hansen@cern.ch
///
///    TTPCFieldCorrector
///   
///    This class uses one of the description of the fields available (for now select in
///    oaCalibrationDatabase) to provide a service to correct for deviations of reconstructed 
///    TPC hits due to inhomogenous fields. 
///
///    Field Corrections are given by the use of nominal field maps, i.e. the best known 
///    representations of the true magnetic and electric fields (see magnetCalib). 
///    Drift electrons' paths through the drift volumes are simulated by solving RK4 
///    equations for steps with adaptive step sizes through the inhomogenetical fields.
///    This corrects the error imposed when assuming perfect E and B fields.
///    This class should normally not be called directly, but one can turn on/off 
///    Field Correction from TDisplacmentCorrector.

class TTPCFieldCorrector {
  
public:
  TTPCFieldCorrector(bool isMC);
  ~TTPCFieldCorrector();

  ///CORRECTION WITH DIRECT DRIFT
  ///Takes a 4D position as input. This should normaly be a point (x,y,z,t) on a helix
  ///given from a set of track parameters for which the likelihood for the 
  ///observed charge distribution will be calculated. 
  ///As an output this function gives the hit on the padplane (i.e. x' = x_padplane,y',z',t') 
  ///after drift under the effect of E and B fields. To be used by the "direct approach",
  ///i.e. by the maximum likelihood fit. Since the likelihood reconstruction gets a 
  ///very good seed track before the fit this only need to be called once 
  ///per each track point (not for every small change in the track parameters 
  ///during the optimization).
  bool GetHitOnPadPlane(const TLorentzVector& trackPoint, TLorentzVector& padPoint);

  //CORRECTION WITH REVERSE DRIFT (will take a point in pad plane and drift it back to the cathode)-> For Laser

  bool GetHitOnCathode(const TVector3& padPos, TVector3 &catPos);


  ///CORRECTION WITH REVERSE DRIFT (will have to be tested that the B field has a 
  ///                               one-to-one mapping before this is used; 
  ///                               tested for NOMAD and SIMFIELD for a set of track points.) 
  ///Takes a 3D position as input. This should normaly be a reconstructed point,
  ///reconstructed assuming constant drift velocity from the data base and
  ///assumed perfect (only x components) B and E fields.       
  ///As an output it gives a new position for the reconstructed point, using the            
  ///inhomogenetic fields. To be used in the "inverse approach", i.e. by the           
  ///chi-square fit.                                                                 
  bool GetFieldCorrectedPoint(const TVector3& recPointWithHomogenousFields, TVector3 &origPos);

  ///Same as above but the input has not assumed any used drift velocity, but                   
  ///gives the drift time instead                                                   
  bool GetFieldCorrectedPoint(const TVector3& padPlanePoint, double time, TVector3 &origPos);

  /// Sets the member variable fDriftVelocity.
  /// Added so the value originally set from runtime parameters (see constructor)
  /// can be replaced with a value from the calibration database.
  void SetDriftVelocity(double driftVelocity) {
    fDriftVelocity = driftVelocity;
  };

  /// This is the access point to turn on or off the EField distortion
  void SetDoEFieldDistortion(bool isEFieldDistOn) {
    fDriftManager.SetDoEFieldDistortion(isEFieldDistOn);
  };
  
  /// If you use the EField model = 0, you can create TEFieldModel
  /// derived object and assigned them to each TPC drift volume.
  /// tpc =  0, 1, 2
  /// rp = 0, 1
  void SetEFieldModel( unsigned int tpc, unsigned int rp, ND::TEFieldModel *Model){
    fDriftManager.SetEFieldModel(tpc, rp ,Model);
  };
  
private: 
  int fBFieldType;
  double fBFieldStrengthX;
  double fBFieldStrengthY;
  double fBFieldStrengthZ;

  double fCathode;
  int fUseAdaptiveStepSize;

  /// Value of the DriftVelocity from the data base or parameters file.
  double fDriftVelocity;

  double fMaxDriftGlobalX;

  ND::TDriftManager fDriftManager;

};
  
#endif // #ifdef TTPCFieldCorrector_hxx
