#ifndef TTPCHitPad_hxx_seen
#define TTPCHitPad_hxx_seen

#include <iostream>
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTrueTrack.hxx"

using namespace std;

namespace trex {
  class TTPCHitPad;
}

/// Contains all the information relative to a waveform from a TPC pad.
/// This may contains more than one peak if the peaks are too closed
/// to be split into individual waveforms.
class trex::TTPCHitPad {
public:

  TTPCHitPad() :fChargeFit(0),fTimeFit(0),fPosition(0,0,0), TrueTrack(0), fTrueTrackID(-1) {};

  ~TTPCHitPad(){}

  TTPCHitPad(double eDep, TLorentzVector pos4);

  double GetCharge(){return fChargeFit;}

  double GetTime(){return fTimeFit;}

  /// Short cut to access the Y position of the pad
  double Y(){return fPosition.Y();};
  /// Short cut to access the Z position of the pad
  double Z(){return fPosition.Z();};

  TVector3 GetPosition(){return fPosition;}
  
  //Each HitPad points at a TrueTrack object that hold track-level truth information
  void SetTrueTrack(trex::TTrueTrack* track) {TrueTrack = track;fTrueTrackID=track->GetTrackID();}
  trex::TTrueTrack* GetTrueTrack() {return TrueTrack;}

  int GetTrueTrackID() { return fTrueTrackID;}

  char Print(){

    std::cout << "This Hit is at position: " << Y() << " : " << Z() << " and Charge: " << GetCharge()  << std::endl;
    char out = 'd';
    return out;
  }
  

private:

  /// Amplitude for the peaks, as produced by the analytic gaussian fit (UNFIT if not fit)
  double fChargeFit;

  /// Fitted times for the peaks, as produced by the analytic gaussian fit (UNFIT if not fit)
  double fTimeFit;
  TVector3 fPosition;

  TTrueTrack* TrueTrack; //!
  int fTrueTrackID;
  
};

#endif
