#ifndef TTRExPath_HXX
#define TTRExPath_HXX


#include "TTRExJunction.hxx"
#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"

namespace trex{
  
  struct TTPCLogLikelihood{
    double Total;
    double X;
    double HV;
  };

  class TTRExJunction;
  
  class TTPCPathFitResults {
  public:
    TTPCPathFitResults(){
      Sigma = 0.0;
      eSigma = 0.0;
      Curvature = 0.0;
      eCurvature = 0.0;
      LogLikelihood.Total = 0.0;
      LogLikelihood.X = 0.0;
      LogLikelihood.HV = 0.0;
      fitSteps = 0;
      IsFitReliable = false;
    }
    std::vector<double> FitState;
    bool IsFitReliable;
    double Sigma;
    double eSigma;
    double Curvature;
    double eCurvature;
    TTPCLogLikelihood LogLikelihood;
    unsigned int fitSteps;
  };
  
  
  class TTRExPath {
    
  public:
    
    TTRExPath() : fClusters(0), fHasChi2Fit(false), fHasRunFit(false), fHasLikelihoodFit(false),fHasFitState(false){};
    
    TTRExPath(std::vector<trex::TTRExHVCluster*> clusters) : fHasChi2Fit(false), fHasRunFit(false), fHasLikelihoodFit(false),fHasFitState(false){
      fClusters = clusters;
    }    
    
    void SetClusters(std::vector<trex::TTRExHVCluster*> clusters){
      fClusters = clusters;
    }
        
    void SetConnectedJunctions(std::vector<trex::TTRExJunction*> &juncts);
    
    void AddConnectedJunction(trex::TTRExJunction* junct);
    
    unsigned int GetId(){return fId;}
    
    void SetId(unsigned int id){fId=id;}
    
    
    std::vector<trex::TTRExHVCluster*>&  GetClusters(){
      return fClusters;
    }
    
    
    void Print() {
      for(unsigned int i=0; i<fClusters.size(); ++i){
	fClusters[i]->Print();
      }
    }
    
    void SaveFitState(TTPCPathFitResults& results){fHasFitState=true;fFitState=results;}
    
    bool HasFitState(){return fHasFitState;}
    
    int GetConnectedEnd(unsigned int junctionId){return 0;}
    
    std::vector<double> GetFrontFitState(){return std::vector<double>(0);}
    std::vector<double> GetBackFitState(){return std::vector<double>(0);}
    
    std::vector<double> GetFrontSeedState(){return std::vector<double>(0);}
    std::vector<double> GetBackSeedState(){return std::vector<double>(0);}
    
    bool HasSeedState(){return false;}
    
    void SaveSeedStates(std::vector<double>& frontSeedState, std::vector<double>& backSeedState){
      fFrontSeedState=frontSeedState;
      fBackSeedState=backSeedState;
    }
    
    void SaveMatchedPath(unsigned int id,trex::TTPCLogLikelihood likelihood){}
    void SaveMatchedPattern(unsigned int patternId, unsigned int pathId, trex::TTPCLogLikelihood likelihood){}
    
    bool IsFrontConnected(){return false;}
    bool IsBackConnected(){return false;}
    
    unsigned int GetNMatchedPattern(){return 0;}
    
    unsigned int GetMatchPatternId(unsigned int n){return 0;}
    
    unsigned int GetPatternMatchPathId(unsigned int n){return 0;}
    
    unsigned int NbEndsFreeToMatch(){return 0;}
    
    unsigned int GetMatchPathIdIndex(unsigned int n){return 0;}
    
    double GetPathMatchLikelihood(unsigned int n){return 0.;}
    
    bool IsEndFreeToMatch(int end){return false;}
    
    double GetPatternMatchLikelihood(unsigned int n){return 0.;}
    
    double GetLogLikelihood(){return 0.;}
    
    void SetEndNotFreeToMatch(int end){}
    
    bool HasChi2Fit(){return fHasChi2Fit;}
    
    bool HasRunFit(){return fHasRunFit;}
    
    bool HasLikelihoodFit(){return fHasLikelihoodFit;}
    
    void SetHasChi2Fit(bool hasChi2Fit){fHasChi2Fit=hasChi2Fit;}
    
    void SetHasRunFit(bool hasRunFit){fHasRunFit=hasRunFit;}
    
    void SetHasLikelihoodFit(bool hasLikelihoodFit){fHasLikelihoodFit=hasLikelihoodFit;}
    
  private:
    
    unsigned int fId;
    
    std::vector<trex::TTRExHVCluster*>  fClusters;
    //INCLUDE MORE TRACKING AND FIT VARIABLES HERE
    
    trex::TTPCPathFitResults fFitState;
    
    bool fHasChi2Fit;
    
    bool fHasRunFit;
    
    bool fHasLikelihoodFit;

    bool fHasFitState;
    
    std::vector<double> fFrontSeedState;
    
    std::vector<double> fBackSeedState;
    
    std::vector<trex::TTRExJunction*> fConnectedJunctions;
    std::vector<unsigned int> fConnectedJunctionsId;
  
  };
}






#endif
