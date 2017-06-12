#ifndef TTRUETRACK_HXX
#define TTRUETRACK_HXX

#include "TVector3.h"

namespace trex{
  
  class TTrueTrack {
    
  public:
    
    //Use constructor-like setter. Values should not be modified once set. 
    TTrueTrack() {}

    void SetEntries(int pdg, int id, int pOrPi, int parent, TVector3 initial, TVector3 final){
      TrackPDG = pdg;
      TrackID = id;
      TrackProOrPi = pOrPi;
      TrackParentID = parent;
      TrackInitialPos = initial;
      TrackFinalPos = final;
    }
    
    int GetTrackPDG() {return TrackPDG;}
    int GetTrackID() {return TrackID;}
    int GetTrackProOrPi() {return TrackProOrPi;}
    int GetTrackParentID() {return TrackParentID;}
    TVector3 GetTrackInitialPos() {return TrackInitialPos;}
    TVector3 GetTrackFinalPos() {return TrackFinalPos;}
    
    
  private:
    
    int TrackPDG;
    int TrackID;
    int TrackProOrPi;
    int TrackParentID;
    TVector3 TrackInitialPos;
    TVector3 TrackFinalPos;
    
  };
  
}
#endif
