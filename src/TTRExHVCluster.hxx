#ifndef TTRExHVCluster_HXX
#define TTRExHVCluster_HXX

//c++
#include <vector>
#include <TVector3>

//TREx
#include "TTPCHitPad.hxx" 


namespace trex{
  
  
  class TTRExHVCluster{

  public:

    TTRExHVCluster() : fcHitPtrs(0), fIsVertical(0), fOkForSeed(0), fPosition(0), fCharge(0) {};

    //TTRExHVCluster(bool HV, std::vector<trex::TTPCHitPad> cHits) {

    //fIsVertical = HV;
    //fcHits = cHits;
      
    // }

    
    TTRExHVCluster(bool HV, std::vector<trex::TTPCHitPad*> cHits){
      
      fIsVertical = HV;
      fcHitPtrs = cHits;
      fOkForSeed = true;
    }
    
    void SetHits(std::vector<trex::TTPCHitPad*> cHits){
      fcHitPtrs = cHits;
    }

    void SetIsVertical(bool HV){
      fIsVertical = HV;
    }

    void SetOkForSeed(bool ok){

      fOkForSeed = ok;
    }

    void SetPosition(TVector3 pos){
    
      fPosition = pos;  
    }

    void SetCharge(double charge){
      
      fCharge = charge;
    }

    bool IsVertical() {

      return fIsVertical;
    }
    
    bool IsHorizontal() {

      return (not fIsVertical);
    }

    bool IsOkForSeed(){
      
      return fOkForSeed;
    }

    std::vector<trex::TTPCHitPad*> GetClusterHits(){
      
      return fcHitPtrs;
    }

    
    int GetNHits(){
      
      return fcHitPtrs.size();
      
    }

    double X() {
      return fPosition.X();
    }

    double Y() {
      return fPosition.Y();
    }

    double Z() {
      return fPosition.Z();
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

    TVector3 fPosition;
    double fCharge;
    std::vector<trex::TTPCHitPad*> fcHitPtrs;
    //std::vector<trex::TTPCHitPad> fcHits;
    bool fIsVertical;
    bool fOkForSeed;
    
  };
  
  
}
  
#endif
