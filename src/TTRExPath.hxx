#ifndef TTRExPath_HXX
#define TTRExPath_HXX



#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"

namespace trex{
  
  class TTRExPath {
    
  public:
    
    TTRExPath() : fClusters(0) {};
    
    TTRExPath(std::vector<trex::TTRExHVCluster> clusters){
      fClusters = clusters;
    }
    
    void SetClusters(std::vector<trex::TTRExHVCluster> clusters){
      fClusters = clusters;
    }
    
    std::vector<trex::TTRExHVCluster>&  GetClusters(){
      return fClusters;
    }


    void Print() {
      for(int i=0; i<fClusters.size(); ++i){
	fClusters[i].Print();
      }
    }
      
    
  private:
    
    std::vector<trex::TTRExHVCluster>  fClusters;
    //INCLUDE MORE TRACKING AND FIT VARIABLES HERE
    
  };
}






#endif
