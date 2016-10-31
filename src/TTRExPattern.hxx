#ifndef TTRExPattern_HXX
#define TTRExPattern_HXX

#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"


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
      
      TTRExPattern() : fPaths(0),fJunctions(0){};
      
      TTRExPattern(std::vector<std::vector<trex::TTRExHVCluster> > paths, std::vector<std::vector<trex::TTPCHitPad> > junctions){
	  fPaths = paths;
	  fJunctions = junctions;
	};


      std::vector<std::vector<trex::TTRExHVCluster> > GetPaths(){
	return fPaths;
      }

      std::vector<std::vector<trex::TTPCHitPad> > GetJunctions(){
	return fJunctions;
      }

      void SetPaths(std::vector<std::vector<trex::TTRExHVCluster> > paths){
	fPaths = paths;
      }

      void SetJunctions(std::vector<std::vector<trex::TTPCHitPad> > junctions){
	fJunctions = junctions;
      }

      void Clear() {
	fPaths.clear();
	fJunctions.clear();
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
	  for(int j=0; j<fPaths[i].size(); ++j){
	    fPaths[i][j].Print(); 
	  }
	}
	
	
	for(int i=0; i<junctSize; ++i){
	  std::cout << " Junction # " << i << " contains the following Hits: " << std::endl;
          for(int j=0; j<fJunctions[i].size(); ++j){
            fJunctions[i][j].Print();
          }
        }
	
      }

	
      private:
	
	std::vector<std::vector<trex::TTRExHVCluster> > fPaths;
	std::vector<std::vector<trex::TTPCHitPad> > fJunctions;
	
	
    };
	
	
}

      
      
#endif
