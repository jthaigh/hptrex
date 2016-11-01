#ifndef TTRExHVCluster_HXX
#define TTRExHVCluster_HXX

//c++
#include <vector>

//TREx
#include "TTPCHitPad.hxx" 


namespace trex{
  
  
  class TTRExHVCluster{

  public:

    TTRExHVCluster() : fcHitPtrs(0), fIsVertical(0) {};

    //TTRExHVCluster(bool HV, std::vector<trex::TTPCHitPad> cHits) {

    //fIsVertical = HV;
    //fcHits = cHits;
      
    // }

    
    TTRExHVCluster(bool HV, std::vector<trex::TTPCHitPad*> cHits){
      
      fIsVertical = HV;
      fcHitPtrs = cHits;

    }
    
    void SetHits(std::vector<trex::TTPCHitPad*> cHits){
      fcHitPtrs = cHits;
    }

    void SetIsVertical(bool HV){
      fIsVertical = HV;
    }

    bool IsVertical() {

      return fIsVertical;
    }
    
    bool IsHorizontal() {

      return (not fIsVertical);
    }

    std::vector<trex::TTPCHitPad*> GetClusterHits(){
      
      return fcHitPtrs;
    }

    
    int GetNHits(){
      
      return fcHitPtrs.size();
      
    }
    
    //GetOutputCluster() {
      
    //std::vector<trex::TTPCHitPad> hits;
      
    //for(auto id=fcHits.begin(); id!=fcHits.end(), ++id){

    //hits.push_back(*id);
    // }

    //return std::move(hits);

    //}


    void Print() {

      string orientation;
      if(fIsVertical){
	orientation = "vertical";
      }else{
	orientation = "horizontal";
      }
      
      std::cout << "This is a " << orientation << " Cluster containing " << this->GetNHits() << " Hits." << std::endl;
    }


  private:

    std::vector<trex::TTPCHitPad*> fcHitPtrs;
    //std::vector<trex::TTPCHitPad> fcHits;
    bool fIsVertical;
    
    
  };
  
  
}
  
#endif
