#ifndef TTRExPIDAlgorithm_HXX
#define TTRExPIDAlgorithm_HXX

#include "TTRExPattern.hxx"
#include "TTRExPath.hxx" 
#include "TTRExHVCluster.hxx"
#include "TTrueTrack.hxx"

namespace trex{
  
  class TTRExPIDAlgorithm{
    
  public:
    
    TTRExPIDAlgorithm(); 

    void Process(std::vector<trex::TTRExPattern>& allPatterns);
    
    
    void dEdx(trex::TTRExPath& path);
    void ChargeSum(trex::TTRExPath& path);
    void TrackLength(trex::TTRExPath& path);
    
    void TrueTrackLength(trex::TTRExPath& path);
    void TrackCleanliness(trex::TTRExPath& path);
    void PID(trex::TTRExPath& path);

  private:
    
    int fTrueEventMultiplicity;
    int fRecoEventMultiplicity;
    //double fPIDcut=3.5e-6;
    
  };
    
}

#endif 
