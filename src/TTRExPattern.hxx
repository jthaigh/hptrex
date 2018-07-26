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
    
    TTRExEvent(std::vector<trex::TTRExPattern>& patterns){
      fPatterns=patterns;
    };

    std::vector<trex::TTRExPattern>& GetPatterns(){
      return fPatterns;
    }

    void Clear() {
      
      fPatterns.clear();

    }

    
  private:
    
    std::vector<trex::TTRExPattern> fPatterns;
        
  };
    
    
    
    class TTRExPattern {
      
    public:
      
      TTRExPattern() : fPaths(0),fJunctions(0),fJunctionsToPathsMap(0), fIsUsable(true){};
      
      TTRExPattern(std::vector<trex::TTRExPath>& paths, std::vector<trex::TTRExJunction>& junctions, std::vector< std::vector<unsigned int> > map){
	  fPaths = paths;
	  fJunctions = junctions;
	  fJunctionsToPathsMap = map;
	};


      std::vector<trex::TTRExPath>& GetPaths(){
	return fPaths;
      }

      std::vector<trex::TTRExJunction>& GetJunctions(){
	return fJunctions;
      }

      std::vector< std::vector<unsigned int> >& GetMap(){
	return fJunctionsToPathsMap;
      }


      void SetPaths(std::vector<trex::TTRExPath> paths){
	fPaths = paths;
      }

      void SetJunctions(std::vector<trex::TTRExJunction> junctions){
	fJunctions = junctions;
      }

      void SetMap(std::vector< std::vector<unsigned int> > map){
	fJunctionsToPathsMap = map;
      }

      void Clear() {
	fPaths.clear();
	fJunctions.clear();
	fJunctionsToPathsMap.clear();
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
      std::vector< std::vector<unsigned int> > fJunctionsToPathsMap;
      bool fIsUsable;

    };

  struct WritablePattern{
    std::vector<std::vector<trex::TTPCHitPad> > Paths;
    std::vector<std::vector<trex::TTPCHitPad> > Junctions;
    std::vector<std::vector<unsigned int> > JunctionsToPathsMap;
    std::vector<double> dEdx;
    std::vector<double> TrackLength;
    std::vector<double> ChargeSum;
    std::vector<int> PID;
    std::vector<double> TrueTrackLength;
    std::vector<double> TrackCompleteness;
    std::vector<double> TrackCleanliness;
    std::vector<int> NumberOfTrueHitsFound;
    std::vector<int> PDG;
    std::vector<TVector3> InitialPosition;
    std::vector<TVector3> FinalPosition;
    std::vector<double> Momentum;
    std::vector<int> TrackNumber;
    std::vector<int> TrackID;
    std::vector<int> ParentID;
   // std::vector<int> ProOrPi;
    std::vector<int> TrueNumberOfHits;
    std::vector<int> NParticles;
  };

  struct WritableEvent{
    std::vector<WritablePattern> patterns;
    
    //Add event level truth information here
    //Event Multiplicity
    int RecoMultiplicity=0;
    //True Event multiplicity
    int TrueMultiplicity=0;
    

    void FillFromEvent(trex::TTRExEvent& evt){
      patterns.clear();
      RecoMultiplicity=0;
     
      for(auto iPat=evt.GetPatterns().begin();iPat!=evt.GetPatterns().end();++iPat){
	
	patterns.emplace_back();

	std::cout << "This pattern contains "<< iPat->GetPaths().size() << " Paths" << std::endl;
	patterns.back().JunctionsToPathsMap=iPat->GetMap();
	RecoMultiplicity+=iPat->GetPaths().size();
	
	for(auto iPath=iPat->GetPaths().begin();iPath!=iPat->GetPaths().end();++iPath){
	  
	  patterns.back().Paths.emplace_back();
	  
	  patterns.back().dEdx.push_back(iPath->GetdEdx());
	  patterns.back().TrackLength.push_back(iPath->GetTrackLength());
	  patterns.back().ChargeSum.push_back(iPath->GetChargeSum());
	  patterns.back().PID.push_back(iPath->GetPID());
	  patterns.back().TrueTrackLength.push_back(iPath->GetTrueTrackLength());
	  patterns.back().TrackCompleteness.push_back(iPath->GetTrackCompleteness());
	  patterns.back().TrackCleanliness.push_back(iPath->GetTrackCleanliness());
	  patterns.back().NumberOfTrueHitsFound.push_back(iPath->GetNumberOfTrueHitsFound());
	  patterns.back().PDG.push_back(iPath->GetPDG());
	  patterns.back().InitialPosition.push_back(iPath->GetInitialPosition());
	  patterns.back().FinalPosition.push_back(iPath->GetFinalPosition());
	  patterns.back().Momentum.push_back(iPath->GetMomentum());
	  patterns.back().TrackNumber.push_back(iPath->GetTrackNumber());
	  patterns.back().TrackID.push_back(iPath->GetTrackID());
	  patterns.back().ParentID.push_back(iPath->GetParentID());
	//  patterns.back().ProOrPi.push_back(iPath->GetProOrPi());
	  patterns.back().TrueNumberOfHits.push_back(iPath->GetTrueNumberOfHits());
	  patterns.back().NParticles.push_back(iPath->GetNParticles());
	  //Add more PID and Truth Variables here
	  
	  
	  for(auto iCluster=iPath->GetClusters().begin();iCluster!=iPath->GetClusters().end();++iCluster){
	    
	    for(auto iHit=(*iCluster)->GetClusterHits().begin();iHit!=(*iCluster)->GetClusterHits().end();++iHit){
	      patterns.back().Paths.back().emplace_back(**iHit);
	    }
	  } 
	}
	std::cout << "check3" << std::endl;
	  
	for(auto iJunct=iPat->GetJunctions().begin();iJunct!=iPat->GetJunctions().end();++iJunct){
	    
	  patterns.back().Junctions.emplace_back();
	  
	  for(auto iHit=iJunct->GetHits().begin();iHit!=iJunct->GetHits().end();++iHit){
	    patterns.back().Junctions.back().emplace_back(**iHit);
	  }	    
	}
	
      }
    }

    void SetTrueMultiplicity(int multi){
      TrueMultiplicity=multi;
    }

  };	
  
}
    
      
      
#endif
