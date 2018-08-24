#include "TTPCHitPad.hxx"

//ClassImp(trex::TTPCHitPad);


trex::TTPCHitPad::TTPCHitPad(double eDep, TLorentzVector pos4): TrueTrack(0), fTrueTrackID(-1){

  fChargeFit=eDep;
  fPosition=pos4.Vect();
  fTimeFit=pos4.T();
}
