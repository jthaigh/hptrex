#ifndef TTRExHVCluster_HXX
#define TTRExHVCluster_HXX

//c++
#include <vector>

//TREx
#include <TTPCHitPad.hxx> 


namespace trex{
  
  
  class TTRExHVCluster{

  public:

    TTRExHVCluster() : fcHits(0), fIsVertical(0) {};

    TTRExHVCluster(bool HV, std::vector<trex::TTPCHitPad> cHits) {

      fIsVertical = HV;
      fcHits = cHits;
      
    }

    
    TRExHVCluster(bool HV, std::vector<trex::TTPCHitPad*> cHits){
      
      fIsVertical = HV;
      fcHits = cHits;

    }
    
    SetHits(std::vector<trex::TTPCHitPad*> cHits){
      fcHits = cHits;
    }

    SetIsVertical(bool HV){
      fIsVertical = HV;
    }

    IsVertical() {

      return fIsVertical;
    }
    
    IsHorizontal() {

      return (not fIsVertical);
    }

    GetCluster(){
      
      return fcHits;
    }
    
    GetOutputCluster() {
      
      std::vector<trex::TTPCHitPad> hits;
      
      for(auto id=fcHits.begin(); id!=fcHits.end(), ++id){

	hits.push_back(*id);
      }

      return std::move(hits);

    }




  private:

    std::vector<trex::TTPCHitPad*> fcHitPtrs;
    //std::vector<trex::TTPCHitPad> fcHits;
    bool fIsVertical;
    
    
  }
  
