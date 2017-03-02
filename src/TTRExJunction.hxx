#ifndef TTRExJunction_HXX
#define TTRExJunction_HXX

#include "TTPCHitPad.hxx"
#include "TTRExPath.hxx"


namespace trex{

  class TTRExPath;
  
  class TTRExJunction{
    
  public:
    
    
    TTRExJunction(){
      //std::cout<<"Created a junction!"<<std::endl;
    }

    TTRExJunction(const std::vector<trex::TTPCHitPad*> &Hits){
      SetHits(Hits);
      //std::cout<<"Created a junction!"<<std::endl;
    }
    
    TTRExJunction(const TVector3 &Position){
      fPosition=Position;
      //std::cout<<"Created a junction!"<<std::endl;
    }
    
    virtual ~TTRExJunction(){};
    
    void SetHits(const std::vector<TTPCHitPad*>& theHits);
    void SetId(unsigned int theId){fId=theId;}
    unsigned int GetId(){return fId;}
        
    std::vector<trex::TTPCHitPad*>& GetHits(){
      return fHits;
    }

    TVector3 GetPosition(){return fPosition;}
    
    /// Get number of paths associated with this junction
    unsigned int GetNPaths(){
      return fConnectedPaths.size();
    }
        
    void Print(){
      for(unsigned int i=0; i<fHits.size(); ++i){
	fHits[i]->Print();
      }   
    }

    std::vector< trex::TTRExPath* >& GetConnectedPaths(){
      return fConnectedPaths;
    }
    
    /// Check if a path with this Id is connected to this junction
    bool IsPathConnected(unsigned int WantedPathId);
    
    void SetConnectedPaths(std::vector<trex::TTRExPath*> &paths);
  
    void AddConnectedPath(trex::TTRExPath* path);

    
  private:
    unsigned int fId;
    std::vector<trex::TTRExPath*> fConnectedPaths;
    std::vector<unsigned int> fConnectedPathsId;
    std::vector<trex::TTPCHitPad*>  fHits;
    TVector3 fPosition;
  };
}




#endif
