#ifndef TTRExPath_HXX
#define TTRExPath_HXX


#include "TTRExJunction.hxx"
#include "TTPCHitPad.hxx"
#include "TTRExHVCluster.hxx"
#include "TTPCHelixPropagator.hxx"

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
  
  /// Structure to store the likelihood when matching paths.                                                                                  
  struct PathMatchInfo {
  public:
    int PathId;
    TTPCLogLikelihood MatchLikelihood;
  };

  /// Structure to store the likelihood when matching patterns.                                                                               
  struct PatternMatchInfo {
  public:
    int PatternId;
    int PathId;
    TTPCLogLikelihood MatchLikelihood;
  };
  
  
  class TTRExPath {
    
  public:
    
    TTRExPath() : fClusters(0), fHasChi2Fit(false), fHasRunFit(false), fHasLikelihoodFit(false),fHasFitState(false), fId(0),
		  fFrontIsConnected(false), fBackIsConnected(false){};
    
    TTRExPath(std::vector<trex::TTRExHVCluster*> clusters) : fHasChi2Fit(false), fHasRunFit(false), fHasLikelihoodFit(false),fHasFitState(false), fId(0),fFrontIsConnected(false), fBackIsConnected(false){
      fClusters = clusters;
    }
    
    
    void SetClusters(std::vector<trex::TTRExHVCluster*> clusters){
      fClusters = clusters;
    }
    
    std::vector<trex::TTRExJunction*> GetConnectedJunctions(){
      return fConnectedJunctions;
    }
    
    std::vector<unsigned int> GetConnectedJunctionsId(){
      return fConnectedJunctionsId;
    }
    
    //Adds an individual connected junction
    void AddConnectedJunction(trex::TTRExJunction* junct);
    
    unsigned int GetId(){return fId;}
    
    void SetId(unsigned int id){fId=id;}
    
    std::vector<trex::TTRExHVCluster*>&  GetClusters(){
      return fClusters;
    }
    
    //
    //    std::vector<trex::TTRExHVCluster*>& GetHits(){
    //  return GetClusters();
    //}
    
    
    //print method for debugging
    void Print() {
      for(unsigned int i=0; i<fClusters.size(); ++i){
	fClusters[i]->Print();
      }
    }
    
    
    //MDH TODO: Implement these
    
    //Convert the end of the path into end nodes
    void SetEndClustersToNodes();
    
    
    //Store the fit result. This method should probably doing a lot more!
    void SaveFitState(TTPCPathFitResults& results);
    
    bool HasFitState(){return fHasFitState;}
    
    
    void SaveFitState(std::vector<double> inState); 
    
    
    //bool HasFitState();
    
    bool HasReliableFitState();
    
    //P.D.: Need to check this works with how ids get assigned
    int GetConnectedEnd(unsigned int junctionId);
    
    std::vector<double> GetFitState();
    
    TTPCPathFitResults GetFitResults(){
      return fFitState;
    }
    
    std::vector<double> GetFrontFitState();
    std::vector<double> GetBackFitState();
    
    std::vector<double> GetFrontSeedState();
    std::vector<double> GetBackSeedState();
    
    //originally a call to the ReconBase
    bool HasSeedState(){return fHasChi2Fit;}
    
    //save the front and back seed states
    void SaveSeedStates(std::vector<double>& frontSeedState, std::vector<double>& backSeedState);
    
    
    void SaveMatchedPath(unsigned int mPatternId,trex::TTPCLogLikelihood matchLklhd);
    void SaveMatchedPattern(unsigned int mPatternId, unsigned int mPathId, trex::TTPCLogLikelihood matchLklhd);
    
    bool IsFrontConnected(){return fFrontIsConnected;}
    bool IsBackConnected(){return fBackIsConnected;}
    
    
    unsigned int GetNMatchedPattern();
    
    unsigned int GetMatchPatternId(unsigned int n);
    
    unsigned int GetPatternMatchPathId(unsigned int i);
    
    unsigned int NbEndsFreeToMatch();
    
    unsigned int GetNMatchedPath();
    
    unsigned int GetMatchPathIdIndex(unsigned int pathId);
    
    double GetPathMatchLikelihood(unsigned int n);
    
    bool IsEndFreeToMatch(unsigned int end);
    
    double GetPatternMatchLikelihood(unsigned int n);
    
    double GetLogLikelihood();
    
    
    void SetEndNotFreeToMatch(int end);
    
    
    bool HasChi2Fit(){return fHasChi2Fit;}
    
    bool HasRunFit(){return fHasRunFit;}
    
    bool HasLikelihoodFit(){return fHasLikelihoodFit;}
    
    void SetHasChi2Fit(bool hasChi2Fit){fHasChi2Fit=hasChi2Fit;}
    
    void SetHasRunFit(bool hasRunFit){fHasRunFit=hasRunFit;}
    
    void SetHasLikelihoodFit(bool hasLikelihoodFit){fHasLikelihoodFit=hasLikelihoodFit;}
    
    void SetdEdx(double dEdx){fdEdx=dEdx;}
    double GetdEdx(){return fdEdx;}
    
    void SetChargeSum(double charge){fChargeSum=charge;}
    double GetChargeSum(){return fChargeSum;}
    
    void SetTrackLength(double length){fTrackLength=length;}
    double GetTrackLength(){return fTrackLength;}
    
  private:
    
    unsigned int fId = 0;
    
    std::vector<trex::TTRExHVCluster*>  fClusters;
    
    //INCLUDE MORE TRACKING AND FIT VARIABLES HERE
    
    trex::TTPCPathFitResults fFitState;
    
    std::vector<double> fFrontFitState;
    std::vector<double> fBackFitState;
    
    bool fHasChi2Fit;
    
    bool fHasRunFit;
    
    bool fHasLikelihoodFit;
    
    bool fHasFitState;
    
    std::vector<double> fFrontSeedState;
    std::vector<double> fBackSeedState;
    
    bool fFrontIsConnected;
    bool fBackIsConnected;
    
    bool fEndFreeToMatch[2];
    
    std::vector<trex::TTRExJunction*> fConnectedJunctions;
    std::vector<unsigned int> fConnectedJunctionsId;
    
    /// List of connected paths and their matching chi2.                 
    std::vector<PathMatchInfo> fPathsMatched;
    
    /// List of connected paths and their matching chi2.                  
    std::vector<PatternMatchInfo> fPatternsMatched;
    
    unsigned int fPID;
    double fdEdx;
    double fTrackLength;
    double fChargeSum;
    
  };
}

#endif
