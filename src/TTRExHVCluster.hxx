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


    TTRExHVCluster() : fcHitPtrs(0), fIsVertical(0), fPosition(), fCharge(0), fOkForSeed(1), fOkForFit(1), fIsUsable(1), fEndNode(0){};

    
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
      CloseHits();
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

    bool isOkForSeed(){
      
      return fOkForSeed;
    }

    std::vector<trex::TTPCHitPad*>& GetClusterHits(){

      return fcHitPtrs;
    }

    
    int GetNHits(){
      
      return fcHitPtrs.size();
      
    }
    
    void SetEndNode(void){
      fEndNode = true;
    }

    void ClearEndNode(void){
      fEndNode = false;
    }

    bool IsEndNode(void){
      return fEndNode;
    }
    
    void SetUsable(bool isUsable){fIsUsable=isUsable;}

    bool isUsable(){return fIsUsable;}

    double GetCharge(){return fCharge;}
    
    double GetDeltaDrift(){return fDeltaDrift;}

    //MDH TODO: Need to know position of readout plance to do this properly
    //If transverse diffusion is not in simulation, need to remove code which
    //uses this.
    double GetDriftDistance(){return fPosition.X();}

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

    void CloseHits();

    std::vector<trex::TTPCHitPad*> fcHitPtrs;
    bool fIsVertical;
    TVector3 fPosition;
    double fCharge;
    bool fOkForSeed;
    bool fOkForFit;
    bool fIsUsable;
    bool fEndNode;
    double fDeltaDrift;

  };
  
  
}
  
#endif
