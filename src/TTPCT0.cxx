#include "TTPCT0.hxx"
#include "TTPCCalibration.hxx"

#include <THandle.hxx>
#include <TGeomInfo.hxx>
#include <TRealDatum.hxx>
#include <TIntegerDatum.hxx>

//*****************************************************************************
TTPCT0::TTPCT0(double Default){
  Init();
  fT0Range[0] = Default;
  fT0Range[1] = Default;
  fT0 = Default;
}

//*****************************************************************************
TTPCT0::TTPCT0(TTPCT0 &T0Copy){
  Init();

  double tmpMin, tmpMax;
  T0Copy.GetT0Range(tmpMin, tmpMax);
  fT0Range[0] = tmpMin;
  fT0Range[1] = tmpMax;
  fT0 = T0Copy.GetT0();
}

//*****************************************************************************
void TTPCT0::Init(){
  fClosestHitFnd = false;
  fMatchedHitFnd = false;
  fSource = kNoT0src;
  fChi2 = 0xABCDEF;

  int sdim = 7;
  EVector v1 = EVector(sdim,0);
  EMatrix C1 = EMatrix(sdim,sdim,0);
  fRes = HyperVector(v1,C1);
}

//*****************************************************************************
void TTPCT0::LoadMatchedHit(ND::THandle<ND::THit> hitRes, double &chi2, HyperVector &Residuals){
  LoadClosestHit(hitRes, chi2);
  fRes = Residuals;
  fClosestHitFnd = false;
  fMatchedHitFnd = true;
}

//*****************************************************************************
void TTPCT0::LoadClosestHit(ND::THandle<ND::THit> hitRes, double &chi2){
  fHit = hitRes;
  fChi2 = chi2;
  fClosestHitFnd = true;
  fMatchedHitFnd = false;
  fT0 = fHit->GetTime();

  // Code to determine the source FGD based on the hit that provided the t0
  if ((fHit->GetGeomId().GetSubsystemId() == ND::GeomId::Def::kFGD)) {
    if (ND::TGeomInfo::Get().FGD().IsInFGD1(fHit->GetPosition())) {
      fSource = kFGD1T0src;
    } else {
      fSource = kFGD2T0src;
    }
  } else if ((fHit->GetGeomId().GetSubsystemId() == ND::GeomId::Def::kDSECal)) {
    fSource = kDsECALT0src;
  } else if ((fHit->GetGeomId().GetSubsystemId() == ND::GeomId::Def::kTECal)) {

    ND::GeomId::Def::DetectorId det;
    ND::GeomId::Def::ECal::MagnetClams clam;
    ND::GeomId::Def::ECal::ECalModules mod;
    ND::TGeomInfo::ECAL().GetModule(fHit).GetIdentity(det, clam, mod);
    if (mod == ND::GeomId::Def::ECal::kTopModule) {
      fSource = kTopBrECALT0src;
    }
    else if (mod == ND::GeomId::Def::ECal::kBottomModule) {
      fSource = kBotBrECALT0src;
    }
    else { // kSideModule
      if (clam == ND::GeomId::Def::ECal::kNegXClam) {
        fSource = kNegXBrECALT0src;
      }
      else {
        fSource = kPosXBrECALT0src;
      }
    }
    
  } else if ((fHit->GetGeomId().GetSubsystemId() == ND::GeomId::Def::kP0D)) {
    fSource = kP0DT0src;
  } else if ((fHit->GetGeomId().GetSubsystemId() == ND::GeomId::Def::kSMRD)) {
    fSource = kSMRDT0src; 
  }

}

//*****************************************************************************
void TTPCT0::LoadCathodeHit(ND::THandle<ND::THit> hitRes, double T0){
  fHit = hitRes;
  fClosestHitFnd = false;
  fMatchedHitFnd = false;
  fT0 = T0;

  fSource = kCathodeT0src;
}


//*****************************************************************************
void TTPCT0::LoadTimeRange(double firstT0, double lastT0){
  if ( firstT0 < lastT0 ){
    fT0Range[0] = firstT0;
    fT0Range[1] = lastT0;
  } else {
    fT0Range[0] = lastT0;
    fT0Range[1] = firstT0;
  }
}

//*****************************************************************************
void TTPCT0::GetT0Range(double &lower, double &upper){
  lower = fT0Range[0];
  upper = fT0Range[1];
}

//*****************************************************************************
bool TTPCT0::HitTimeWithinT0Range(ND::THandle<ND::THit> candHit){
  // If no T0 range is available, just accept all hits.
  if ( fT0Range[0] == fT0Range[1])
    return true;

  // Add 200ns just to avoid removing time too close to the edge.
  // This is much more than the resolution of the various detectors.
  if ( (fT0Range[0]-200.) < candHit->GetTime() && candHit->GetTime() < (fT0Range[1]+200.) )
    return true;
  return false;
}


//*****************************************************************************
double TTPCT0::GetHitCharge(){
  if (fHit){
    return fHit->GetCharge();
  } else {
    return -0xABCDEF;
  }
}


//*****************************************************************************
void TTPCT0::FillTRealData(ND::THandle<ND::TReconBase> TRecB){
  ND::TIntegerDatum* t0source = new ND::TIntegerDatum("T0Source", fSource);
  TRecB->AddDatum(t0source);

  // When we don't have a T0, we save the possible range of T0s.
  // The analyzer can decide to use it or not.
  if ( fSource == kNoT0src){
    ND::TRealDatum* T0Range = new ND::TRealDatum("T0Range", fT0Range[0]);
    T0Range->push_back(fT0Range[1]);
    TRecB->AddDatum(T0Range);
    // Get first cluster or hitpad for drift sense
    ND::THandle<ND::THit> Hit = *(TRecB->GetHits()->begin());
    double dSense = ND::TGeomInfo::TPC().GetDriftSense(Hit->GetGeomId());
    
    double diffX = dSense * fabs(fT0Range[0] - fT0) * ND::tpcCalibration().GetDriftVelocity();
    ND::TRealDatum* DeltaX = new ND::TRealDatum("T0Range_DeltaX", diffX);
    diffX = dSense * fabs(fT0Range[1] - fT0) * ND::tpcCalibration().GetDriftVelocity();
    DeltaX->push_back(diffX);
    TRecB->AddDatum(DeltaX);

  }
}


//*****************************************************************************
std::string ConvertT0idxToName(int source){
  switch (source)
  {
    case kNoT0src:         return "No T0";
    case kFGD1T0src:       return "FGD1";
    case kFGD2T0src:       return "FGD2";
    case kTopBrECALT0src:  return "Top BrECal";
    case kBotBrECALT0src:  return "Bottom BrECal";
    case kNegXBrECALT0src: return "Negative X BrECal";
    case kPosXBrECALT0src: return "Positive X BrECal";
    case kDsECALT0src:     return "DsECal";
    case kP0DT0src:        return "P0D";
    case kSMRDT0src:       return "SMRD";
    case kCathodeT0src:    return "Cathode";
    default:            return "T0 INDEX IS OUT OF RANGE !";
             break;
  }
}
