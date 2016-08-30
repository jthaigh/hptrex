#ifndef TTRExPattern_HXX
#define TTRExPattern_HXX

#include "TTPCHitPad.hxx"

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
      
      TTRExPattern(std::vector<std::vector<trex::TTPCHitPad> > paths, std::vector<std::vector<trex::TTPCHitPad> > junctions){
	  fPaths = paths;
	  fJunctions = junctions;
	};


      std::vector<std::vector<trex::TTPCHitPad> > GetPaths(){
	return fPaths;
      }

      std::vector<std::vector<trex::TTPCHitPad> > GetJunctions(){
	return fJunctions;
      }

      void SetPaths(std::vector<std::vector<trex::TTPCHitPad> > paths){
	fPaths = paths;
      }

      void SetJunctions(std::vector<std::vector<trex::TTPCHitPad> > junctions){
	fJunctions = junctions;
      }

      void Clear() {
	fPaths.clear();
	fJunctions.clear();
      }
      
	
      private:
	
	std::vector<std::vector<trex::TTPCHitPad> > fPaths;
	std::vector<std::vector<trex::TTPCHitPad> > fJunctions;
	
	
    };
	
	
}

      
      
#endif
