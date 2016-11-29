#ifndef TTRExJunction_HXX
#define TTRExJunction_HXX

#include "TTPCHitPad.hxx"
#include "TTRExPath.hxx"


namespace trex{

  class TTRExPath;
  
  class TTRExJunction{
    
  public:
    
    
    TTRExJunction(){};
    TTRExJunction(const std::vector<trex::TTPCHitPad*> &Hits){
      fHits = Hits;
    }
    
    void SetConnectedPaths(std::vector<trex::TTRExPath*> &paths);
  
    void AddConnectedPath(trex::TTRExPath* path);
    
    TTRExJunction(const TVector3 &Position){};
    
    virtual ~TTRExJunction(){};
    
    void SetHits(std::vector<TTPCHitPad*>& theHits){fHits=theHits;}
    void SetId(unsigned int theId){fId=theId;}
    unsigned int GetId(){return fId;}
        
    std::vector<trex::TTPCHitPad*> GetHits(){
      return fHits;
    }
    
    //MDH TODO: Implement these methods properly
    /// Get number of paths associated with this junction
    unsigned int GetNPaths(){
      return fConnectedPaths.size();
    }
    
    /// Check if a path with this Id is connected to this junction
    bool IsPathConnected(unsigned int WantedPathId);
        
    void Print(){
      for(unsigned int i=0; i<fHits.size(); ++i){
	fHits[i]->Print();
      }   
    }
    
    
  private:
    unsigned int fId;
    std::vector<trex::TTRExPath*> fConnectedPaths;
    std::vector<unsigned int> fConnectedPathsId;
    std::vector<trex::TTPCHitPad*>  fHits;
    
  };
}




#endif
