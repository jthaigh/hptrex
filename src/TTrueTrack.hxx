#ifndef TTRUETRACK_HXX
#define TTRUETRACK_HXX

#include "TVector3.h"

namespace trex{
  
  class TTrueTrack {
    
  public:
    
    //Use constructor-like setter. Values should not be modified once set. 
    TTrueTrack() {}

    //Setter for values that come straight from the input file
    void SetEntries(int pdg, int trackNum, int id,/* int pOrPi,*/ int parent, TVector3 initial, TVector3 final, double mom, int NParticles){
      fTrackPDG = pdg;
      fTrackNumber = trackNum;
      fTrackID = id;
     // fTrackProOrPi = pOrPi;
      fTrackParentID = parent;
      fTrackInitialPos = initial;
      fTrackFinalPos = final;
      fMomentum = mom;
      fNParticles = NParticles;
    }


    void SetNumberOfHits(int hitNum) {fNumberOfHits = hitNum;}
    
    int GetTrackPDG(){return fTrackPDG;}
    int GetTrackNumber(){return fTrackNumber;}
    int GetTrackID(){return fTrackID;}
   // int GetTrackProOrPi(){return fTrackProOrPi;}
    int GetNParticles(){return fNParticles;}
    int GetTrackParentID(){return fTrackParentID;}
    TVector3 GetTrackInitialPos(){return fTrackInitialPos;}
    TVector3 GetTrackFinalPos(){return fTrackFinalPos;}
    double GetMomentum(){return fMomentum;}
    int GetNumberOfHits(){return fNumberOfHits;}

  private:
    
    int fTrackPDG;
    int fTrackNumber;
    int fTrackID;
 //   int fTrackProOrPi;
    int fNParticles;
    int fTrackParentID;
    TVector3 fTrackInitialPos;
    TVector3 fTrackFinalPos;
    double fMomentum;
    int fNumberOfHits;
   
  };
  
}
#endif
