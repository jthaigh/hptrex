#include "TTPCClusterCorrection.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCRecPackUtils.hxx"

#include <TOARuntimeParameters.hxx>
#include <HEPUnits.hxx>
#include <TGeomInfo.hxx>
#include <TGeoManager.h>
#include <iostream>

using namespace std;

// *****************************************************************************
TTPCClusterCorrection::TTPCClusterCorrection() {

///// Create default configuration and set 

  fFieldCorrector = NULL;
  fDistCorrector = NULL;

  fConfig = new CorrectorConfig();
  fConfList.push_back(fConfig);

  fCorrectSeedPosition = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.CorrectSeedPosition");

  fAllDistMap = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.AllDistMap");

  // Distortion Map Name
  TString fDistMapName = ND::TOARuntimeParameters::Get().GetParameterS("trexRecon.Reco.Correction.DistMapName");
  TString fullFileName = getenv("TREXRECONROOT") + TString("/parameters/") + fDistMapName ; 

  // It should be already the case but let's be extra careful
  fConfig = fConfList[0];
  fConfig->SetName("Main");
  // Use this boolean to define complicated modes
  bool doEFieldDistortion;    
  if( !ND::tpcCalibration().IsMC()){
    fConfig->SetDoFieldCorrection(ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.BFieldOnData"));
    fConfig->SetDoDistortionCorrection(ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.EmpiricalDistortionOnData"));
    doEFieldDistortion    = (bool) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.EFieldOnData");
  } else {
    fConfig->SetDoFieldCorrection(ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.BFieldOnMC"));
    fConfig->SetDoDistortionCorrection(ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.EmpiricalDistortionOnMC"));
    doEFieldDistortion    = (bool) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.EFieldOnMC");
  }
  // The Field corrector is absolutely necessary for the EField distortion at the moment
  if (doEFieldDistortion && fConfig->isFieldCorrectorOn())
    fConfig->SetDoEFieldDistortion(true);

  // Create the other configurations

  // ===> Refit with only distortion correction on
  bool RefitDistCorrOnData = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.RefitEmpDistOnData");
  bool RefitDistCorrOnMC = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.RefitEmpDistOnMC");
  if ( ((!ND::tpcCalibration().IsMC()) && RefitDistCorrOnData ) || ( ND::tpcCalibration().IsMC() && RefitDistCorrOnMC ) ){
    CorrectorConfig *thisConfig = new CorrectorConfig();
    fConfList.push_back(thisConfig);
    thisConfig->SetName("DistortionCorr");
    
    thisConfig->SetDoFieldCorrection(false);
    thisConfig->SetDoDistortionCorrection(true);
    thisConfig->SetDoEFieldDistortion(false);
  }

  // ===> Refit with BField and distortion correction on
  bool RefitBFieldAndDistCorrOnData = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.RefitBFieldAndEmpDistOnData");
  bool RefitBFieldAndDistCorrOnMC = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.RefitBFieldAndEmpDistOnMC");
  if ( ((!ND::tpcCalibration().IsMC()) && RefitBFieldAndDistCorrOnData ) || ( ND::tpcCalibration().IsMC() && RefitBFieldAndDistCorrOnMC ) ){
    CorrectorConfig *thisConfig = new CorrectorConfig();
    fConfList.push_back(thisConfig);
    thisConfig->SetName("BFieldAndDistCorr");
    
    thisConfig->SetDoFieldCorrection(true);
    thisConfig->SetDoDistortionCorrection(true);
    thisConfig->SetDoEFieldDistortion(false);
  }

  // ===> Refit with EField distortion on for MC.
  bool RefitEFieldDistOnMC = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.RefitEFieldOnMC");
  if ( ND::tpcCalibration().IsMC() && RefitEFieldDistOnMC ){
    CorrectorConfig *thisConfig = new CorrectorConfig();
    fConfList.push_back(thisConfig);
    thisConfig->SetName("EFieldDistortion");
    
    thisConfig->SetDoFieldCorrection(false);
    thisConfig->SetDoDistortionCorrection(false);
    thisConfig->SetDoEFieldDistortion(true);
  }
  // ===> Refit with EField distortion off for Data.
  bool RefitEFieldDistOffData = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Correction.RefitEFieldOffData");
  if ( (!ND::tpcCalibration().IsMC()) && RefitEFieldDistOffData ){
    CorrectorConfig *thisConfig = new CorrectorConfig();
    fConfList.push_back(thisConfig);
    thisConfig->SetName("EFieldDistortion");
    
    thisConfig->SetDoFieldCorrection(true);
    thisConfig->SetDoDistortionCorrection(false);
    thisConfig->SetDoEFieldDistortion(false);
  }

  bool needFieldCorrector = false;
  bool needDistortionCorrection = false;
  bool needEFieldDistortion = false;
  for ( unsigned int cf = 0; cf < fConfList.size(); cf++ ){
    needFieldCorrector = needFieldCorrector || fConfList[cf]->isFieldCorrectorOn();
    needDistortionCorrection = needDistortionCorrection || fConfList[cf]->isDistortionCorrectorOn();
    needEFieldDistortion = needEFieldDistortion || fConfList[cf]->isEFieldDistortionOn();
  }

  // Just to be very careful, don't even initialize the pointers if we don't do any field or distortion correction.
  if ( needFieldCorrector || needDistortionCorrection){
    fFieldCorrector = new TTPCFieldCorrector(ND::tpcCalibration().IsMC());
    fDistCorrector = new TTPCEmpDistCorrector();
  }

  if( needFieldCorrector ) {
    std::cout << "trexRecon: Field corrector enabled " << std::endl; 
  }


  if ( needDistortionCorrection ) {

    if(fAllDistMap){
      if (!SetAllDistMap()) {
        std::cout<<"ERROR in trexRecon TTPCClusterCorrection::SetAllDistMap, distortion map could not be set! Exit!"<<std::endl;
        exit(1);
      }
      else std::cout<<"trexRecon TTPCClusterCorrection::SetAllDistMap distortion maps opened!"<<std::endl;    

    }else {
      if (!SetDistMap(fullFileName)) {
        std::cout<<"ERROR in trexRecon TTPCClusterCorrection::SetDistMap, distortion map could not be set! Exit!"<<std::endl;
        exit(1);
      }
      else std::cout<<"trexRecon TTPCClusterCorrection::SetDistMap," <<  fullFileName <<" distortion map opened!"<<std::endl;
    }
  }

  if ( needEFieldDistortion ) {
    std::cout << "trexRecon: EField Distortion enabled " << std::endl; 
  }

  std::cout<<"trexRecon::TTPCClusterCorrection configurations"<<std::endl;
  std::cout<<setw(20)<<"Configuration";
  std::cout<<" |  BField Corr.  |  Dist. Corr.  |  EField Dist."<<std::endl;
  for ( unsigned int cf = 0; cf < fConfList.size(); cf++ ){
    std::cout<< setw(20) <<fConfList[cf]->GetName()<<setw(10)<<
    fConfList[cf]->isFieldCorrectorOn()<<setw(15)<<
    fConfList[cf]->isDistortionCorrectorOn()<<setw(15)<<
    fConfList[cf]->isEFieldDistortionOn()<<std::endl;
  }

  // This automatically takes care of special cases like switching the EField On if necessary
  ResetDefaultConfiguration();
  
}



// *****************************************************************************
TTPCClusterCorrection::~TTPCClusterCorrection() {
  if (!fFieldCorrector)
    delete fFieldCorrector;
  if (!fDistCorrector)
    delete fDistCorrector;
}

// *****************************************************************************
void TTPCClusterCorrection::LoadConfiguration(unsigned int cf){
  if (cf >= fConfList.size())
    std::cerr<<"TTPCClusterCorrection::LoadConfiguration(): ERROR. Configuration number "<<cf<<" requested doesn't exist ! It will crash now"<<std::endl;
  fConfig = fConfList[cf];
  // to be sure it crashes right away if the wrong configuration was requested
  fConfig->isFieldCorrectorOn();
  // Special treatment for the EField
  if (fFieldCorrector != NULL)
    fFieldCorrector->SetDoEFieldDistortion(fConfig->isEFieldDistortionOn());
}


// *****************************************************************************
void TTPCClusterCorrection::Apply(ND::THandle<ND::TTPCPath> thePath){
  State propagState = thePath->GetFrontSeedState();

  ND::THandle<ND::THitSelection> HVclu = thePath->GetHits();
  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  ND::rpman("TREx").model_svc().model().intersector().set_length_sign(0);
  for (ND::THitSelection::const_iterator Hit = HVclu->begin(); Hit != HVclu->end(); Hit++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = (*Hit);
    double x_raw, y_raw, z_raw;

    bool RPpropagOk = false;
    double propagLength = 99999999.;
    // Perform the propagation of the seed.
    if ( fCorrectSeedPosition){
      RPpropagOk = TTPCRecPackUtils::PropagateToHVCluster(Cluster, propagState, propagLength);
      x_raw = propagState.vector()[0];
      y_raw = propagState.vector()[1];
      z_raw = propagState.vector()[2];
    }

    // RPpropagOk will be false if fCorrectSeedPosition is false and if the propagation failed.
    // In both cases we will use the cluster position instead.
    if ( (!RPpropagOk) || propagLength > 500. ) {
      x_raw = Cluster->CalibX();
      y_raw = Cluster->Y();
      z_raw = Cluster->Z();
    }

    TVector3 padPoint(x_raw, y_raw, z_raw);
    TVector3 trackPoint;
    GetCorrectedPoint(padPoint,trackPoint); 
    // WARNING: The code in tpcRecon has this sign convention for the DeltaY and DeltaZ.
    // I don't understand it ... why are Y and Z opposite signs ???
    Cluster->SetDeltaYZ( (y_raw - trackPoint.Y()), (trackPoint.Z() - z_raw) );
  }

  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);

}



// *****************************************************************************
bool TTPCClusterCorrection::SetDistMap(TString distMap) {
  return fDistCorrector->SetDistMap(distMap);
}

// *****************************************************************************
bool TTPCClusterCorrection::SetAllDistMap() {
  return fDistCorrector->SetAllDistMap();
}


// *****************************************************************************
bool TTPCClusterCorrection::GetCorrectedPoint(const TVector3& recPos, TVector3& corrPos) {
  TVector3 prevPos = recPos;
  corrPos = recPos;
  bool DidSmthng = false;

  if (fConfig->isFieldCorrectorOn()) {
    if (!fFieldCorrector->GetFieldCorrectedPoint(prevPos,corrPos)){
      prevPos = corrPos;
      DidSmthng = true;
    }
  }
  if (fConfig->isDistortionCorrectorOn()) {
    if (!fDistCorrector->DistMapIsSet()) cout<<"ERROR: Distortion map not set!"<<endl;
    if (!fDistCorrector->GetDistortionCorrectedPoint(prevPos,corrPos)){
      prevPos = corrPos;
      DidSmthng = true;
    }
  }

  return DidSmthng;
}


// *****************************************************************************
//CORRECTION WITH REVERSE DRIFT (will take a point in pad plane and drift it back to the cathode)-> For Laser
bool TTPCClusterCorrection::GetHitOnCathode(const TVector3& padPos, TVector3& corrPos) {
  TVector3 prevPos = padPos;
  corrPos = padPos;
  bool DidSmthng = false;

  if (fConfig->isDistortionCorrectorOn()) {
    if (!fDistCorrector->DistMapIsSet()) cout<<"ERROR: Distortion map not set!"<<endl;
    if (!fDistCorrector->GetDistortionCorrectedPoint(prevPos,corrPos)){
      prevPos = corrPos;
      DidSmthng = true;
    }
  }
  if (fConfig->isFieldCorrectorOn()) {
    if (!fFieldCorrector->GetHitOnCathode(prevPos,corrPos)){
      prevPos = corrPos;
      DidSmthng = true;
    }
  }

  return DidSmthng;
}


// *****************************************************************************
bool TTPCClusterCorrection::GetHitOnPadPlane(const TLorentzVector& trackPoint, TLorentzVector& padPoint) {
  if ((fConfig->isFieldCorrectorOn())&&(fConfig->isDistortionCorrectorOn())) {
    //From the track point see where an electron end up on the pad plane taking (only) the nominal fields into account
    if (!fFieldCorrector->GetHitOnPadPlane(trackPoint, padPoint)) return false;
    //From UV Laser tests and Cosmics tests we know emperically how much wrong the drift simulation with 
    //nominal fields gives (saved in a distortion map). Apply these extra distortions on the
    //field corrected pad point, using the original track point to know what distortions to apply
    if (!fDistCorrector->DistMapIsSet()) cout<<"ERROR: Distortion map not set!"<<endl;
    if (!fDistCorrector->ApplyDistOnPadPoint(trackPoint.Vect(), padPoint)) return false;

    return true;
  } else if (fConfig->isFieldCorrectorOn()) {
    if (!fFieldCorrector->GetHitOnPadPlane(trackPoint, padPoint)) return false;
    return true;
  } else if (fConfig->isDistortionCorrectorOn()) {
    TLorentzVector distortedPoint(0,0,0,0);
    if (!fDistCorrector->DistMapIsSet()) cout<<"ERROR: Distortion map not set!"<<endl;
    if (!fDistCorrector->ApplyDistOnPadPoint(trackPoint.Vect(),distortedPoint)) return false;
    //Just to get the pad plane x coordinate
    //Get a reference to the singleton instance of TPC geometry information
    // Use the TPC pad manager to get the pad geometry
    ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();
    // Get the pad plane x coordinate
    double padPlaneX = tpcGeom.GetMaxDriftDistance() + tpcGeom.GetCathodeWidth() / 2.;
    if (trackPoint.X() < 0) padPlaneX = -padPlaneX;
    //Set the pad point that is not field corrected but only distortion corrected
    padPoint = TLorentzVector(TVector3(padPlaneX,distortedPoint.Y(),distortedPoint.Z()),trackPoint.T());
    cout<<"WARNING in TTPCClusterCorrection::GetHitOnPadPlane: Distortion correction is done but not ";
    cout<<"Field correction, so there is not time (i.e. X) correction!"<<endl;
    return true;
  }
  return false;
}
