#ifndef TTPCClusterCorrection_h
#define TTPCClusterCorrection_h

#include <THandle.hxx>
#include <TVector3.h>

#include "TTPCPath.hxx"
#include "TTPCFieldCorrector.hxx"
#include "TTPCEmpDistCorrector.hxx"


///    \brief Provides services to the TPC fitting methods to correct for spatial distortions
///    \author Christian Hansen, christian.hansen@cern.ch
///
///    TTPCClusterCorrection
///   
///    This class uses TTPCFieldCorrector and TTPCEmpDistCorrector to provide a sum of
///    both corrections as a service to fit methods. It is possible to turn on/off either one
///    of these two sub corrections with SetDoFieldCorrection and SetDoDistortionCorrection.
///
///    Two different fit methods will call the Displacement Corrector; the Likelihood Fit and
///    the Single Point Fit (chi square), and the Displacement Corrector works differently for
///    the two different fit methods. 
///    
///    The Single Point Fit gives a reconstructed point to the Displacement Corrector.
///    This point was reconstructed assuming perfect electrical and magnetic fields and
///    no distortions. The Displacement Corrector tells the Single Point fit where that 
///    point actually belongs taking field and distortion corrections into account.
///    
///    The Likelihood Fit works differently and therefore should call different functions
///    of the Displacement Corrector. The Likelihood Fit is supposed to minimize the 
///    difference between the pad-charge distribution obtained from a track and the 
///    pad-charge distribution seen on pads. So TMinuit is used to suggest a set of track
///    parameters. Points along the helix defined by the track parameters are given to
///    the Displacement Corrector which tells where on the pad place (y, z) and when (x)
///    a drift electron arrives, taking inhomogenous fields and other distortions into 
///    account. Displacement Corrector probably only needs to be called once per track
///    since the seed track for Likelihood fit is a good approximation.
///
///    For more detailed instructions how to use, see mainpage.dox.

class TTPCClusterCorrection {
  
public:
  TTPCClusterCorrection();
  ~TTPCClusterCorrection();

  ///Takes a 3D position as input. This should normaly be a reconstructed point,
  ///reconstructed assuming constant drift velocity from the data base and 
  ///assumed perfect (only x components) B and E fields. 
  ///As an output it gives a new position for the reconstructed point, using
  ///both the TTPCFieldCorrector and TTPCEmpDistCorrector 
  ///To be used in the "inverse approach", i.e. by the chi-square fit.
  bool GetCorrectedPoint(const TVector3& recPointWithHomogenousFields, TVector3& corrPos);

  //CORRECTION WITH REVERSE DRIFT (will take a point in pad plane and drift it back to the cathode)-> For Laser
  bool GetHitOnCathode(const TVector3& padPos, TVector3& corrPos);

  ///Takes a 4D point as input. This should normaly be a point (x,y,z,t) on a helix
  ///given from a set of track parameters (suggested by Minuit for the optimization) 
  ///for which the likelihood for the observed charge distribution will be calculated. 
  ///As an output this function gives the hit on the padplane (i.e. x' = x_padplane,y',z',t') 
  ///when corrections from TTPCFieldCorrector and TTPCEmpDistCorrector has been taken
  ///into account. To be used by the "direct approach", i.e. by the maximum likelihood fit. 
  ///Since the likelihood reconstruction gets a very good seed track before the fit this 
  ///only need to be called once per each track point (not for every small change in the 
  ///track parameters during the optimization).
  bool GetHitOnPadPlane(const TLorentzVector& trackPoint, TLorentzVector& padPoint);

  ///Load the configuration wanted
  void LoadConfiguration(unsigned int cf);
  void ResetDefaultConfiguration() {LoadConfiguration(0);};
  unsigned int GetNbConfigurations() {return fConfList.size();};
  std::string GetConfigurationName() {return fConfig->GetName();};

  void Apply(ND::THandle<ND::TTPCPath>);
  
  bool SetDistMap(TString distMap);
  bool SetAllDistMap();

  bool isFieldCorrectorOn() {return fConfig->isFieldCorrectorOn();};
  bool isDistortionCorrectorOn() {return fConfig->isDistortionCorrectorOn();};
  bool isAnyCorrectionOn() {return fConfig->isAnyCorrectionOn();};  

private: 
  TTPCFieldCorrector* fFieldCorrector;
  TTPCEmpDistCorrector* fDistCorrector;

  /// Switch to correct the seed position or the cluster position for each cluster.
  bool fCorrectSeedPosition;

  // Only one set of maps can be initialized at once so
  // all configurations must use the same
  // Flag to initialize all distortion maps. Default is disable.
  bool fAllDistMap;
  // Name of the distortion map
  TString fDistMapName;

  class CorrectorConfig {
  public:
    CorrectorConfig(){
      fDoFieldCorrection = false;
      fDoDistortionCorrection = false;
      fDoEFieldDistortion = false;
    }
    ~CorrectorConfig() {};

    void SetName(std::string ConfName) {fName = ConfName;};
    void SetDoFieldCorrection(bool doFieldCorr) {fDoFieldCorrection = doFieldCorr;};
    void SetDoDistortionCorrection(bool doDistCorr) {fDoDistortionCorrection = doDistCorr;};
    void SetDoEFieldDistortion(bool doEFieldDist) {fDoEFieldDistortion = doEFieldDist;};

    std::string GetName()         {return fName;};
    bool isFieldCorrectorOn()     {return fDoFieldCorrection;};
    bool isDistortionCorrectorOn(){return fDoDistortionCorrection;};
    bool isEFieldDistortionOn()   {return fDoEFieldDistortion;};
    bool isAnyCorrectionOn()      {return ( fDoFieldCorrection || fDoDistortionCorrection || fDoEFieldDistortion);};  
    
  private:
    // A string to name the configuration
    std::string fName;
    // Flag to enable the correction. Default is disabled.
    bool fDoFieldCorrection;
    // Flag to enable the distortion correction. Default is disabled.
    bool fDoDistortionCorrection;
    // Flag to enable the EField correction. Default is disabled.
    // This is in simulate mode for MC, i.e. this will add the effect
    // of the EField distortion on a perfect MC.
    // Of course it will correct on data.
    bool fDoEFieldDistortion;
  };

  std::vector<CorrectorConfig*> fConfList;
  CorrectorConfig *fConfig;

};
  
#endif // #ifdef TTPCClusterCorrection_hxx
