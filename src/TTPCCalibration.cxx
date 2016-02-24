#include "TTPCCalibration.hxx"
#include <TOARuntimeParameters.hxx>
#include <TGeomInfo.hxx>
#include <TTPC_Drift_Velocity_Table.hxx>
#include <TEventFolder.hxx>
#include <TTPCHitPad.hxx>
#include <TTPCHVCluster.hxx>

/// The static member pointer to the singleton.
ND::TTPCCalibration* ND::TTPCCalibration::_tpcCalibration = NULL;

//*****************************************************************************
ND::TTPCCalibration& ND::tpcCalibration(){
//*****************************************************************************

  return ND::TTPCCalibration::Get();

}

//*****************************************************************************
ND::TTPCCalibration& ND::TTPCCalibration::Get(void) {
//*****************************************************************************
  if (!_tpcCalibration) {
    _tpcCalibration = new ND::TTPCCalibration();
  }

  return *_tpcCalibration;
}

//*****************************************************************************
ND::TTPCCalibration::TTPCCalibration(){
//*****************************************************************************
  const ND::TND280Event& Event = *(ND::TEventFolder::GetCurrentEvent());
  const ND::TND280Context& ND280cont = Event.GetContext();
  fEvtIsMC = ND280cont.IsMC();

  fLastEvent = -1;
  fLastTime = -1;

  fDriftVelocity = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCgas.DriftSpeed");;
  fUseMeasDriftSpeed = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.TPCgas.UseMeasDriftSpeed");

  fTimeOffset_MC = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.TimeOffsetMC");
  fTimeOffset_Data = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.TimeOffset");

  fDefaultCosmT0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.DefaultCosmT0");
  fDefaultSpillT0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.DefaultSpillT0");
  std::cout<<"TTPCT0Finder Default Spill T0 = "<<fDefaultSpillT0<<" Default Cosmic T0 = "<<fDefaultCosmT0<<std::endl;

  fMaxDrift = ND::TGeomInfo::Get().TPC().GetMaxDriftDistance(0,0);
}


//*****************************************************************************
void ND::TTPCCalibration::ReadCalibration(const ND::TND280Event& Event) {
//*****************************************************************************
  const ND::TND280Context& ND280cont = Event.GetContext();

  fPathIdOffset = 1;
  fJunctionIdOffset = 1;
  fPatternIdOffset = 1;

  fTimeOffset = fTimeOffset_Data;

  // Obtain calibration constant for drift velocity if event is not from MC:
  if (fEvtIsMC) {
    fTimeOffset = fTimeOffset_MC;
    return;
  }

  // This event was already calibrated. Don't call the database again.
  if ( fLastEvent == ND280cont.GetEvent() && fLastTime == ND280cont.GetTimeStamp() )
    return;

  fLastEvent = ND280cont.GetEvent();
  fLastTime = ND280cont.GetTimeStamp();

  // TPC calibrated drift velocity
  if (fUseMeasDriftSpeed) {
    ND::TResultSetHandle<ND::TTPC_Drift_Velocity_Table> rs(ND280cont);
    if (rs.GetNumRows() != 0) {

      fDriftVelocity = rs.GetRow(0)->GetDriftVelocity();
      // Convert from calib table units of cm/us to mm/ns
      fDriftVelocity = fDriftVelocity / 100.0;

    }
  }

  // First set a default value from the parameters file
  fDefaultT0 = fDefaultSpillT0;
  Int_t tbits = Event.GetHeader().GetTriggerBits();
  if ( (tbits & ND::Trigger::kTFBCosmic) || (tbits & ND::Trigger::kFGDCosmic) ) fDefaultT0 = fDefaultCosmT0;

}



//*****************************************************************************
bool ND::TTPCCalibration::IsMC(){
//*****************************************************************************
  return fEvtIsMC;
}



//*****************************************************************************
double ND::TTPCCalibration::GetDriftVelocity(){
//*****************************************************************************
  return fDriftVelocity;
}


//*****************************************************************************
double ND::TTPCCalibration::GetTimeOffset(){
//*****************************************************************************
  return fTimeOffset;
}


//*****************************************************************************
unsigned int ND::TTPCCalibration::GetPathId(){
//*****************************************************************************
  return fPathIdOffset++;
}

//*****************************************************************************
unsigned int ND::TTPCCalibration::GetJunctionId(){
//*****************************************************************************
  return fJunctionIdOffset++;
}

//*****************************************************************************
unsigned int ND::TTPCCalibration::GetPatternId(){
//*****************************************************************************
  return fPatternIdOffset++;
}


//*****************************************************************************
void ND::TTPCCalibration::OffsetIds(){
//*****************************************************************************
  // Ids start at 1 and therefore at 1001
  fPatternIdOffset  = TREXSTDOUTIDOFFSET+1;
  fPathIdOffset     = TREXSTDOUTIDOFFSET+1;
  fJunctionIdOffset = TREXSTDOUTIDOFFSET+1;
}


//*****************************************************************************
void ND::TTPCCalibration::SetDefaultT0(double DftT0){
//*****************************************************************************
  fDefaultT0 = DftT0;
}

//*****************************************************************************
double ND::TTPCCalibration::GetDefaultT0(){
//*****************************************************************************
  return fDefaultT0;
}
