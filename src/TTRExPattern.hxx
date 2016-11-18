#ifndef TTRExPattern_HXX
#define TTRExPattern_HXX

#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"
#include "TTRExPath.hxx"
#include "TTRExJunction.hxx"

namespace trex{
  
  class TTRExPattern;
  
  
  class TTRExEvent {
    
  public: 
    
    TTRExEvent() : fPatterns(0) {};
    
    TTRExEvent(std::vector<trex::TTRExPattern> patterns){
      fPatterns=patterns;
    };

    std::vector<trex::TTRExPattern> GetPatterns(){
      return fPatterns;
    }

    void SetPatterns(std::vector<trex::TTRExPattern> patterns){
      fPatterns = patterns;
    }


    void Clear() {
      
      fPatterns.clear();

    }

    
  private:
    
    std::vector<trex::TTRExPattern> fPatterns;
    
    
  };
    
    
    
    class TTRExPattern {
      
    public:
      
      TTRExPattern() : fPaths(0),fJunctions(0),fPathsToJunctionMap(0), fIsUsable(true){};
      
      TTRExPattern(std::vector<trex::TTRExPath>& paths, std::vector<trex::TTRExJunction>& junctions, std::vector< std::vector<unsigned int> > map){
	  fPaths = paths;
	  fJunctions = junctions;
	  fPathsToJunctionMap = map;
	};


      std::vector<trex::TTRExPath>& GetPaths(){
	return fPaths;
      }

      std::vector<trex::TTRExJunction>& GetJunctions(){
	return fJunctions;
      }

      std::vector< std::vector<unsigned int> > GetMap(){
	return fPathsToJunctionMap;
      }


      void SetPaths(std::vector<trex::TTRExPath> paths){
	fPaths = paths;
      }

      void SetJunctions(std::vector<trex::TTRExJunction> junctions){
	fJunctions = junctions;
      }

      void SetMap(std::vector< std::vector<unsigned int> > map){
	fPathsToJunctionMap = map;
      }

      void Clear() {
	fPaths.clear();
	fJunctions.clear();
	fPathsToJunctionMap.clear();
      }
      

      void PrintSize() {

	int pathSize = fPaths.size();
	int junctSize = fJunctions.size();

       
	std::cout << "This Pattern holds " << pathSize << " Paths and " << junctSize << " Junctions" << std::endl;

      }

      void Print() {
	
	int pathSize = fPaths.size();
	int junctSize = fJunctions.size();

	std::cout << "This Pattern has the following content: " << std::endl;
	
	for(int i=0; i<pathSize; ++i){
	  std::cout << " Path # " << i << " contains the following Hits: " << std::endl;
	  fPaths[i].Print();	  
	}
	
	for(int i=0; i<junctSize; ++i){
	  std::cout << " Junction # " << i << " contains the following Hits: " << std::endl;
	  const std::vector<trex::TTPCHitPad*> jHits;
          for(auto iHit=jHits.begin(); iHit!=jHits.end(); ++iHit){
            (**iHit).Print();
          }
	  }
	
      }

      unsigned int GetId(){return fId;}

      void SetId(unsigned int id){fId=id;}
      
      void SetUsable(bool isUsable){fIsUsable=isUsable;}

      bool IsUsable(){return fIsUsable;}

    private:
      
      unsigned int fId;
      std::vector<trex::TTRExPath> fPaths;
      std::vector<trex::TTRExJunction> fJunctions;
      std::vector< std::vector<unsigned int> > fPathsToJunctionMap;
      bool fIsUsable;

    };
	
	
}

      
      
#endif
