#ifndef TTRExPattern_HXX
#define TTRExPattern_HXX

#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"
#include "TTRExPath.hxx"

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
      
      TTRExPattern() : fPaths(0),fJunctions(0),fPathsToJunctionMap(0){};
      
      TTRExPattern(std::vector<trex::TTRExPath> paths, std::vector<std::vector<trex::TTPCHitPad> > junctions, std::vector< std::vector<unsigned int> > map){
	  fPaths = paths;
	  fJunctions = junctions;
	  fPathsToJunctionMap = map;
	};


      std::vector<trex::TTRExPath>& GetPaths(){
	return fPaths;
      }

      std::vector<std::vector<trex::TTPCHitPad> >& GetJunctions(){
	return fJunctions;
      }

      std::vector< std::vector<unsigned int> > GetMap(){
	return fPathsToJunctionMap;
      }


      void SetPaths(std::vector<trex::TTRExPath> paths){
	fPaths = paths;
      }

      void SetJunctions(std::vector<std::vector<trex::TTPCHitPad> > junctions){
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
          for(int j=0; j<fJunctions[i].size(); ++j){
            fJunctions[i][j].Print();
          }
        }
	
      }
      
      
    private:
      
      std::vector<trex::TTRExPath> fPaths;
      std::vector<std::vector<trex::TTPCHitPad> > fJunctions;
      std::vector< std::vector<unsigned int> > fPathsToJunctionMap;
	
    };
	
	
}

      
      
#endif
