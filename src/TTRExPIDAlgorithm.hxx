#ifndef TTRExPIDAlgorithm_HXX
#define TTRExPIDAlgorithm_HXX

#include "TTRExPattern.hxx"
#include "TTRExPath.hxx" 
#include "TTRExHVCluster.hxx"

namespace trex{
  
  class TTRExPIDAlgorithm{
    
  public:
    
    TTRExPIDAlgorithm(); 

    void Process(std::vector<trex::TTRExPattern>& allPatterns);
    
    
    void dEdx(trex::TTRExPath& path);
    void ChargeSum(trex::TTRExPath& path);
    void TrackLength(trex::TTRExPath& path);
    
    
  private:
    
    
    
  };
    
}

#endif 
