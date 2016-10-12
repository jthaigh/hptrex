
#ifndef TTRUEHIT_HXX              
#define TTRUEHIT_HXX       

#include "TLorentzVector.h"

class TTrueHit {

public:

  TTrueHit() : TrueEdep(0),TruePos4(0,0,0,0),pdg(0),TrueTrackID(0),charge(0) {}

  TTrueHit(double edep, TLorentzVector pos, int Pdg, int trackId, int Charge){
    TrueEdep=edep;
    TruePos4=pos;
    pdg=Pdg;
    TrueTrackID=trackId;
    charge=Charge;
  }
  
  bool operator <(TTrueHit a){

    if(TrueTrackID < a.TrueTrackID){
      return true;}
    else{return false;}

  }

  double TrueEdep;
  TLorentzVector TruePos4;
  int pdg;
  int TrueTrackID;
  int charge;

};

#endif
