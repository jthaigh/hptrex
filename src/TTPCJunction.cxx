#include "TTPCJunction.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCUtils.hxx"
#include "TTPCCalibration.hxx"

#include <iostream>
#include <TrackingUtils.hxx>
#include <TMultiHit.hxx>
#include <TIntegerDatum.hxx>
#include <TRealDatum.hxx>

// *********************************************************************************
ClassImp(ND::TTPCJunction);
ND::TTPCJunction::TTPCJunction() : TReconVertex(){
  fId = 0;
  fT0.SetDefaultT0(0.0);
}


// *********************************************************************************
ND::TTPCJunction::TTPCJunction(const TLorentzVector &PosTime) : TReconVertex(){
  ND::THandle<ND::TVertexState> tstate = this->GetState();
  tstate->SetPosition(PosTime);
  fId = 0;
  fT0.SetDefaultT0(0.0);
}


// *********************************************************************************
ND::TTPCJunction::~TTPCJunction() { }


// *********************************************************************************
void ND::TTPCJunction::SetId(unsigned int theId){
  fId = theId;
}


// *********************************************************************************
unsigned int ND::TTPCJunction::GetId(){
  return fId;
}


// *********************************************************************************
void ND::TTPCJunction::InitialSetup(){
  if(! GetHits()){
    std::cout<<"TTPCJunction::InitialSetup: WARNING. No hits available. Setup aborted."<<std::endl;
    return;
  }
  ND::THandle<ND::TTPCHitPad> hPad = *(GetHits()->begin());
  if(! hPad){
    std::cout<<"TTPCJunction::InitialSetup: WARNING. The first THit is not a TTPCHVCluster. Setup aborted."<<std::endl;
    return;
  }
  AddDetector(TrackingUtils::GetDetector(hPad));

}


// *********************************************************************************
void ND::TTPCJunction::SetT0(TTPCT0 &T0){
  fT0 = T0;
  ND::THandle<ND::TVertexState> tstate = this->GetState();
  TLorentzVector PosTime = tstate->GetPosition();
  ND::THandle<ND::TTPCHitPad> hitPad = *(GetHits()->begin());
  double Sense = hitPad->DriftSense();
  double RPx = hitPad->GetPosition().X();
  // How far are we from the read out plane ?
  double DriftDistance = (PosTime.T() - fT0.GetT0() - ND::tpcCalibration().GetTimeOffset()) * ND::tpcCalibration().GetDriftVelocity();
  PosTime.SetX( RPx - (Sense * DriftDistance));
  PosTime.SetT( T0.GetT0());
  tstate->SetPosition(PosTime); 
}

// *********************************************************************************
bool ND::TTPCJunction::HasT0(){
  return (fT0.GetSource() != kNoT0src);
}

// *********************************************************************************
TTRExT0Source ND::TTPCJunction::GetT0Source(){
  return fT0.GetSource();
}

// *********************************************************************************
double ND::TTPCJunction::GetT0(){
  return fT0.GetT0();
}


// *********************************************************************************
double ND::TTPCJunction::GetCalibX(ND::THandle<ND::THit> rawHit){
  ND::THandle<ND::TTPCHitPad> Hit = rawHit;
  double Sense = Hit->DriftSense();
  double RPx = Hit->GetPosition().X();
  // How far are we from the read out plane ?
  double DriftDistance = (Hit->GetTime() - fT0.GetT0() - ND::tpcCalibration().GetTimeOffset()) * ND::tpcCalibration().GetDriftVelocity();
  return( RPx - (Sense * DriftDistance));
  
}


// *********************************************************************************
bool ND::TTPCJunction::IsPathConnected(unsigned int WantedPathId){
  for (ND::TReconObjectContainer::iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
    ND::THandle< ND::TTPCPath > Path = (*ptc);
    if (WantedPathId == Path->GetId() )
      return true;
  }
  return false;
}



// *********************************************************************************
ND::THandle<ND::TReconVertex> ND::TTPCJunction::ConvertToOAEvent() {
  // TReconVertex doesn't have a copy constructor and the TReconBase copy constructor
  // doesn't copy the nodes and the state so we must do it by hand.
  ND::THandle<ND::TReconVertex> tvertex = TTPCUtils::CopyCreateTReconVertex( this);

  // The waveforms are not automatically converted so we have to do it by hand.
  ND::THandle<ND::THitSelection> oldHits = GetHits();
  ND::THitSelection* newHits = new ND::THitSelection();
  if(! oldHits){
    // TODO: Proper exception
    throw;
    return tvertex;
  }
  for (ND::THitSelection::const_iterator tmpHit = oldHits->begin(); tmpHit != oldHits->end(); tmpHit++) {
    ND::THandle<ND::TTPCHitPad> oHit = *tmpHit;
    ND::THandle<ND::TMultiHit> nHit = oHit->ConvertToOAEvent();
    newHits->push_back(nHit);
  }
  tvertex->AddHits(newHits);

  ND::TIntegerDatum* junctionId = new ND::TIntegerDatum("JunctionId",fId);
  tvertex->AddDatum(junctionId);

  fT0.FillTRealData(tvertex);

  return tvertex;
}

