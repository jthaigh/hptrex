#ifndef TTPCPath_hxx_seen
#define TTPCPath_hxx_seen

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconPID.hxx>
#include <TReconNode.hxx>
#include <recpack/State.h>

#include "TTPCJunction.hxx"
#include "TTPCT0.hxx"
#include "TrackingUtils.hxx"

namespace ND {
  class TTPCPath;
}

struct TTPCLogLikelihood{
    double Total;
    double X;
    double HV;
};

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
    State FitState;
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

class ND::TTPCPath : public TReconTrack {
public:
  TTPCPath();
  TTPCPath(unsigned int Id);
  virtual ~TTPCPath();

  void SetId(unsigned int theId);
  unsigned int GetId();

  /// Setup things like the detector bit for the path.
  /// This should be called right after AddHits.
  void InitialSetup(unsigned int PatternId);

  TVector3 GetFirstPosition();
  TVector3 GetLastPosition();

  double GetFirstTime();
  double GetLastTime();

  void AddJunctionId(unsigned int Id);
  void ClearJunctionIds();
  std::vector<unsigned int> GetJunctionIds();
  unsigned int GetNJunctionIds();
  int GetConnectedEnd(unsigned int JunctionId);

  /// Creates the propagation surfaces.
  /// This must be called before any propagation to an HVcluster.
  void PrepareForPropagation();
  /// Call RemovePropagSurf for each cluster.
  void CleanUp();

  void SetFrontConnection(bool Connect);
  void SetBackConnection(bool Connect);
  bool IsFrontConnected();
  bool IsBackConnected();
  void SaveSeedStates(const State &frontState, const State &backState);
  bool HasSeedState();
  State GetFrontSeedState();
  State GetBackSeedState();

  void SetEndNotFreeToMatch(unsigned int End);
  bool IsEndFreeToMatch(unsigned int End);
  unsigned int NbEndsFreeToMatch();

  void SetCathodeCrosser(bool CathCross);
  bool IsCathodeCrosser();
  
  void SetLength(double length);
  double GetLength();

  void SaveFitState( State inState);
  void SaveFitState( TTPCPathFitResults &fitRes);
  TTPCPathFitResults GetFitResults();
  double GetLogLikelihood();
  bool HasFitState();
  bool HasReliableFitState();
  State GetFitState();
  State GetFrontFitState();
  State GetBackFitState();


  void SetT0(TTPCT0 &T0);
  bool HasT0();
  TTRExT0Source GetT0Source();
  double GetT0();
  const TTPCT0& GetTTPCT0();

  void SaveMatchedPath(int mPathId, TTPCLogLikelihood matchLklhd);
  int GetNMatchedPath();
  double GetPathMatchLikelihood(int i);
  int GetMatchPathId(int i);
  int GetMatchPathIdIndex(int pathId);

  void SaveMatchedPattern(int mPatternId, int mPathId, TTPCLogLikelihood matchLklhd);
  int GetNMatchedPattern();
  double GetPatternMatchLikelihood(int i);
  int GetMatchPatternId(int i);
  int GetPatternMatchPathId(int i);
  int GetMatchPatternIdIndex(int patternId);

  void SetEndClustersToNodes();

  void SetPID(ND::TReconPID::ParticleId pid, double weight);
  void SaveInRealDatum(std::string name, double value);

  ND::THandle<ND::TReconBase> ConvertToOAEvent();

  /// Get whether path is in x direction and requiring special treatment
  bool GetIsXPath(){ return fIsXPath; }
  /// Set whether path is in x direction and requiring special treatment
  void SetIsXPath(bool isXPath){ fIsXPath = isXPath; }

private:
  void Init();

  void PrepareTReconTrack( State HelixStart, State HelixEnd, ND::THandle<ND::TReconTrack> Track);

  ND::THandle<ND::TReconBase> GetOAEventObj();

  TTPCT0 fT0;

  /// Whether path is in x direction and requiring special treatment
  bool fIsXPath;

  State fFrontSeedState;
  State fBackSeedState;

  bool fFrontIsConnected;
  bool fBackIsConnected;

  bool fEndFreeToMatch[2];

  State fFrontFitState;
  State fBackFitState;
  TTPCPathFitResults fFitResults;

  TTPCTrackType fTrackType;

  double fLength;
  double fChi2;
  double fNDOF;

  unsigned int fId;
  unsigned int fPatternId;
  std::vector<unsigned int> fJunctionId;

  ND::TReconPID::ParticleId fPID;
  double fPIDweight;

  /// List of connected paths and their matching chi2.
  std::vector<PathMatchInfo> fPathsMatched;

  /// List of connected paths and their matching chi2.
  std::vector<PatternMatchInfo> fPatternsMatched;

  /// Name of the variables and their values.
  /// Each name will be saved in a TRealDatum with the corresponding value.
  std::map<std::string,double> fInRealDatum;

  ND::THandle<ND::TReconBase> fOAEventObj;

  ClassDef(TTPCPath,1);
};


#endif
