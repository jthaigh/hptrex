#ifndef TTRExHVCluster_HXX
#define TTRExHVCluster_HXX

//c++
#include <vector>

//root
#include <TVector3.h>

//TREx
#include "TTPCHitPad.hxx" 

namespace trex{
  
  
  class TTRExHVCluster{

  public:

    TTRExHVCluster() : fcHitPtrs(0), fIsVertical(0), fPosition(), fCharge(0), fOkForSeed(1), fOkForFit(1), fIsUsable(1){};
    
    TTRExHVCluster(bool HV, std::vector<trex::TTPCHitPad*> cHits){
      
      fIsVertical = HV;
      fcHitPtrs = cHits;
      fOkForSeed = true;
    }

    double X(){return fPosition.X();}

    double Y(){return fPosition.Y();}

    double Z(){return fPosition.Z();}
    
    TVector3 GetPosition(){return fPosition;}

    bool isOkForFit(){return fOkForFit;}

    void SetOkForFit(bool isOk){fOkForFit=isOk;}

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

    void SetUsable(bool isUsable){fIsUsable=isUsable;}

    bool isUsable(){return fIsUsable;}

    //MDH TODO: Implement these
    double GetCharge(){return fCharge;}
    
    double GetDeltaDrift(){return 0.;}

    double GetDeltaY(){return 0.;}

    double GetDeltaZ(){return 0.;}

    double GetDriftDistance(){return 0.;}

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

    //MDH TODO: Implement this
    //This will calculate mean cluster position
    void CloseHits(){}

    std::vector<trex::TTPCHitPad*> fcHitPtrs;
    bool fIsVertical;
    TVector3 fPosition;
    double fCharge;
    //std::vector<trex::TTPCHitPad> fcHits;
    bool fOkForSeed;
    bool fOkForFit;
    bool fIsUsable;

  };
  
  
}
  
#endif
