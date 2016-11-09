#include "TTPCPath.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCUtils.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCCalibration.hxx"

#include <TOARuntimeParameters.hxx>
#include <TIntegerDatum.hxx>
#include <TRealDatum.hxx>
#include <TrackingUtils.hxx>
#include <TRecPackManager.hxx>
#include <RecPackConverters.hxx>

#include <sstream>
#include <cmath>

// *********************************************************************************
ClassImp(ND::TTPCPath);
ND::TTPCPath::TTPCPath() : TReconTrack(){
  fId = 0;
  fPatternId = 0;
  Init();

}


// *********************************************************************************
ND::TTPCPath::TTPCPath(unsigned int Id) : TReconTrack(){
  fId = Id;
  fPatternId = 0;
  Init();

}


// *********************************************************************************
void ND::TTPCPath::Init(){
  fIsXPath = false;

  fLength = 0.0;
  fChi2 = 0.0;
  fNDOF = 0.0;

  fPID = ND::TReconPID::kNotSet;

  fFrontIsConnected = false;
  fBackIsConnected = false;

  fEndFreeToMatch[0] = true;
  fEndFreeToMatch[1] = true;

  fTrackType = kUNKNOWNTRACK;

}


// *********************************************************************************
ND::TTPCPath::~TTPCPath() { }


// *********************************************************************************
void ND::TTPCPath::SetId(unsigned int theId){
  fId = theId;
}


// *********************************************************************************
unsigned int ND::TTPCPath::GetId(){
  return fId;
}


// *********************************************************************************
void ND::TTPCPath::InitialSetup(unsigned int PatternId){
  if(! GetHits()){
    std::cout<<"TTPCPath::InitialSetup: WARNING. No hits available. Setup aborted."<<std::endl;
    return;
  }
  ND::THandle<ND::TTPCHVCluster> cluster = *(GetHits()->begin());
  if(! cluster){
    std::cout<<"TTPCPath::InitialSetup: WARNING. The first THit is not a TTPCHVCluster. Setup aborted."<<std::endl;
    return;
  }
  if(! cluster->GetHits().size()){
    std::cout<<"TTPCPath::InitialSetup: WARNING. No hits available in the first cluster. Setup aborted."<<std::endl;
    return;
  }
  AddDetector(TrackingUtils::GetDetector(*(cluster->GetHits().begin())));
  fPatternId = PatternId;

  PrepareForPropagation();
}



// *********************************************************************************
// TODO: use an exception instead
TVector3 ND::TTPCPath::GetFirstPosition(){
  if (this->size() == 0){
    std::cerr<<"TTPCPath ERROR: No THits in TTPCPath !"<<std::endl;
    return TVector3(0.,0.,0.);
  }
  ND::THitSelection::const_iterator it = this->GetHits()->begin();
  return (*it)->GetPosition();
}


// *********************************************************************************
// TODO: use an exception instead
TVector3 ND::TTPCPath::GetLastPosition(){
  if (this->size() == 0){
    std::cerr<<"TTPCPath ERROR: No THits in TTPCPath !"<<std::endl;
    return TVector3(0.,0.,0.);
  }
  ND::THitSelection::const_reverse_iterator it = this->GetHits()->rbegin();
  return (*it)->GetPosition();
}


// *********************************************************************************
// TODO: use an exception instead
double ND::TTPCPath::GetFirstTime(){
  if (this->size() == 0){
    std::cerr<<"TTPCPath ERROR: No THits in TTPCPath !"<<std::endl;
    return 0.;
  }
  ND::THitSelection::const_iterator it = this->GetHits()->begin();
  return (*it)->GetTime();
}


// *********************************************************************************
// TODO: use an exception instead
double ND::TTPCPath::GetLastTime(){
  if (this->size() == 0){
    std::cerr<<"TTPCPath ERROR: No THits in TTPCPath !"<<std::endl;
    return 0.;
  }
  ND::THitSelection::const_reverse_iterator it = this->GetHits()->rbegin();
  return (*it)->GetTime();
}


// *********************************************************************************
void ND::TTPCPath::AddJunctionId(unsigned int Id) {
  fJunctionId.push_back(Id);
}


// *********************************************************************************
void ND::TTPCPath::ClearJunctionIds() {
  fJunctionId.clear();
}


// *********************************************************************************
std::vector<unsigned int> ND::TTPCPath::GetJunctionIds() {
  return fJunctionId;
}


// *********************************************************************************
unsigned int ND::TTPCPath::GetNJunctionIds() {
  return fJunctionId.size();
}


// *********************************************************************************
int ND::TTPCPath::GetConnectedEnd(unsigned int JunctionId) {
  if (fJunctionId.size() == 1){
    if ( fFrontIsConnected)
      return -1;
    else
      return 1;
  }
  if ( (*(fJunctionId.begin()) == JunctionId))
    return -1;
  else if ( (*(fJunctionId.rbegin()) == JunctionId))
    return 1;
  else
    return 0;
}


// *********************************************************************************
void ND::TTPCPath::PrepareForPropagation(){
  unsigned countClu = 0;
  std::stringstream ss;
  for (ND::THitSelection::const_iterator itClu = GetHits()->begin() ; itClu != GetHits()->end(); itClu++, countClu++) {
    // Use the Path ID to make sure that the cluster name is unique.
    // It doesn't matter if this cluster ends up in a merged path later
    // because the name is stored in the cluster.
    ss << "path"<< fId << "clu" << countClu;
    ND::THandle<ND::TTPCHVCluster> Clu = *itClu;
    Clu->PreparePropagSurf(ss.str());
    ss.str(std::string());
  }
}


// *********************************************************************************
void ND::TTPCPath::CleanUp(){
  for (ND::THitSelection::const_iterator itClu = GetHits()->begin() ; itClu != GetHits()->end(); itClu++) {
    ND::THandle<ND::TTPCHVCluster> Clu = *itClu;
    Clu->RemovePropagSurf();
  }
}


// *********************************************************************************
void ND::TTPCPath::SetFrontConnection(bool Connect){
  fFrontIsConnected = Connect;
  if (fFrontIsConnected)
    fEndFreeToMatch[0] = false;
  else
    fEndFreeToMatch[0] = true;
}

// *********************************************************************************
void ND::TTPCPath::SetBackConnection(bool Connect){
  fBackIsConnected = Connect;
  if (fBackIsConnected)
    fEndFreeToMatch[1] = false;
  else
    fEndFreeToMatch[1] = true;
}


// *********************************************************************************
bool ND::TTPCPath::IsFrontConnected(){
  return fFrontIsConnected;
}

// *********************************************************************************
bool ND::TTPCPath::IsBackConnected(){
  return fBackIsConnected;
}

// *********************************************************************************
void ND::TTPCPath::SetEndNotFreeToMatch(unsigned int End){
  // TODO: Proper exception needed
  if (End > 1)
    throw;
  fEndFreeToMatch[End] = false;
}

// *********************************************************************************
bool ND::TTPCPath::IsEndFreeToMatch(unsigned int End){
  // TODO: Proper exception needed
  if (End > 1)
    throw;
  return fEndFreeToMatch[End];
}

// *********************************************************************************
unsigned int ND::TTPCPath::NbEndsFreeToMatch(){
  unsigned int count = 0;
  for ( int i = 0; i < 2; i++)
    if (fEndFreeToMatch[i])
      count++;
  return count;
}

// *********************************************************************************
void ND::TTPCPath::SaveSeedStates(const State &frontState, const State &backState){
  SetStatus(ND::TReconBase::kChi2Fit);
  fFrontSeedState = frontState;
  fBackSeedState = backState;
}

// *********************************************************************************
bool ND::TTPCPath::HasSeedState(){
  return ( CheckStatus(ND::TReconBase::kChi2Fit));
}

// *********************************************************************************
State ND::TTPCPath::GetFrontSeedState(){
  if ( CheckStatus(ND::TReconBase::kChi2Fit)){
    return fFrontSeedState;
  }
  return State();
}


// *********************************************************************************
State ND::TTPCPath::GetBackSeedState(){
  if ( CheckStatus(ND::TReconBase::kChi2Fit)){
    return fBackSeedState;
  }
  return State();
}


// *********************************************************************************
void ND::TTPCPath::SetCathodeCrosser(bool CathCross){
  if (CathCross)
    fTrackType = kCATHCROSSTRACK;
  else
    fTrackType = kSTANDARDTRACK;
}


// *********************************************************************************
bool ND::TTPCPath::IsCathodeCrosser(){
  return (fTrackType == kCATHCROSSTRACK);
}


// *********************************************************************************
void ND::TTPCPath::SetLength(double length){
  fLength = length;
}

// *********************************************************************************
double ND::TTPCPath::GetLength(){
  return fLength;
}


// *********************************************************************************
void ND::TTPCPath::SaveFitState( State inState){
  TTPCPathFitResults tmpFitRes;
  tmpFitRes.FitState = inState;
  tmpFitRes.IsFitReliable = true;
  SaveFitState(tmpFitRes);
}

// *********************************************************************************
void ND::TTPCPath::SaveFitState( TTPCPathFitResults &fitRes){
  // Create the front and back states
  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  // Just to be sure that we don't miss the first cluster if the state was defined
  // at the second or third.
  ND::rpman("TREx").model_svc().model().intersector().set_length_sign(0);

  // This is the running state at each cluster
  State localState = fitRes.FitState;
  fChi2 = 0.0;

  ND::THitSelection::const_iterator tmpTRExClu = GetHits()->begin();
  bool frontStateDone = false;
  State LastState;
  for ( ; tmpTRExClu != GetHits()->end(); tmpTRExClu++) {
    //------- Propagate the helix to each cluster --------
    ND::THandle<ND::TTPCHVCluster> CurrentTRExClu = *tmpTRExClu;
    // If propagation fails, let's try our luck on the next cluster.
    if (!TTPCRecPackUtils::PropagateToHVCluster(CurrentTRExClu, localState))
      continue;

    if (!frontStateDone){
      frontStateDone = true;
      fFrontFitState = localState;
    }
    fBackFitState = localState;
  }

  // If the propagation of the front state failed, this is probably a bad fit od a delta ray.
  if (!frontStateDone)
    return;    

  SetStatus(ND::TReconBase::kLikelihoodFit);
  ClearStatus(ND::TReconBase::kSuccess);
  if (fitRes.IsFitReliable)
    SetStatus(ND::TReconBase::kSuccess);
  fFitResults = fitRes;

  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);
}

// *********************************************************************************
TTPCPathFitResults ND::TTPCPath::GetFitResults(){
  return fFitResults;
}

// *********************************************************************************
double ND::TTPCPath::GetLogLikelihood(){
  return fFitResults.LogLikelihood.Total;
}

// *********************************************************************************
bool ND::TTPCPath::HasFitState(){
  return ( CheckStatus(ND::TReconBase::kLikelihoodFit) &&
       CheckStatus(ND::TReconBase::kSuccess));
}

// *********************************************************************************
bool ND::TTPCPath::HasReliableFitState(){
  return ( CheckStatus(ND::TReconBase::kLikelihoodFit) &&
       CheckStatus(ND::TReconBase::kSuccess) && fFitResults.IsFitReliable);
}

// *********************************************************************************
State ND::TTPCPath::GetFitState(){
  if ( CheckStatus(ND::TReconBase::kLikelihoodFit) &&
       CheckStatus(ND::TReconBase::kSuccess)){
    return fFrontFitState;
  }
  return State();
}

// *********************************************************************************
State ND::TTPCPath::GetFrontFitState(){
  if ( CheckStatus(ND::TReconBase::kLikelihoodFit) &&
       CheckStatus(ND::TReconBase::kSuccess)){
    return fFrontFitState;
  }
  return State();
}


// *********************************************************************************
State ND::TTPCPath::GetBackFitState(){
  if ( CheckStatus(ND::TReconBase::kLikelihoodFit) &&
       CheckStatus(ND::TReconBase::kSuccess)){
    return fBackFitState;
  }
  return State();
}


// *********************************************************************************
void ND::TTPCPath::SetT0(TTPCT0 &T0){
  fT0 = T0;
  // Save the old position of the first cluster
  ND::THitSelection::const_iterator tmpTRExClu = GetHits()->begin();
  ND::THandle<ND::TTPCHVCluster> TRExClu = *tmpTRExClu;
  double OldX = TRExClu->CalibX();
  for ( ; tmpTRExClu != GetHits()->end(); tmpTRExClu++) {
    TRExClu = *tmpTRExClu;
    TRExClu->SetT0(fT0.GetT0());
  }

  // By how much has X changed ?
  tmpTRExClu = GetHits()->begin();
  TRExClu = *tmpTRExClu;
  double DeltaX = TRExClu->CalibX() - OldX;

  // Readjust the state X positions when we set the new T0.
  if (HasSeedState()){
    EVector vect = fFrontSeedState.vector();
    EMatrix cova = fFrontSeedState.matrix();
    vect[0] += DeltaX;
    fFrontSeedState.set_hv(HyperVector(vect,cova));
    vect = fBackSeedState.vector();
    cova = fBackSeedState.matrix();
    vect[0] += DeltaX;
    fBackSeedState.set_hv(HyperVector(vect,cova));
  }
  // Readjust the seed X position now that we have the T0
  if (HasFitState()){
    EVector vect = fFrontFitState.vector();
    EMatrix cova = fFrontFitState.matrix();
    vect[0] += DeltaX;
    fFrontFitState.set_hv(HyperVector(vect,cova));
    vect = fBackFitState.vector();
    cova = fBackFitState.matrix();
    vect[0] += DeltaX;
    fBackFitState.set_hv(HyperVector(vect,cova));
  }
}

// *********************************************************************************
bool ND::TTPCPath::HasT0(){
  return (fT0.GetSource() != kNoT0src);
}

// *********************************************************************************
TTRExT0Source ND::TTPCPath::GetT0Source(){
  return fT0.GetSource();
}

// *********************************************************************************
double ND::TTPCPath::GetT0(){
  return fT0.GetT0();
}

// *********************************************************************************
const TTPCT0& ND::TTPCPath::GetTTPCT0(){
  return fT0;
}


// *********************************************************************************
void ND::TTPCPath::SaveMatchedPath(int mPathId, TTPCLogLikelihood matchLklhd){
  PathMatchInfo newMatch;
  newMatch.PathId = mPathId;
  newMatch.MatchLikelihood = matchLklhd;
  fPathsMatched.push_back(newMatch);
}


// *********************************************************************************
int ND::TTPCPath::GetNMatchedPath(){
  return fPathsMatched.size();
}


// *********************************************************************************
double ND::TTPCPath::GetPathMatchLikelihood(int i){
  return fPathsMatched[i].MatchLikelihood.Total;
}


// *********************************************************************************
int ND::TTPCPath::GetMatchPathId(int i){
  return fPathsMatched[i].PathId;
}

// *********************************************************************************
int ND::TTPCPath::GetMatchPathIdIndex(int pathId){
  int idx = 0;
  for( std::vector<PathMatchInfo>::iterator pmi = fPathsMatched.begin(); pmi != fPathsMatched.end(); pmi++){
    if ( (*pmi).PathId == pathId)
      return idx;
    idx++;
  }
  return -1;
}



// *********************************************************************************
void ND::TTPCPath::SaveMatchedPattern(int mPatternId, int mPathId, TTPCLogLikelihood matchLklhd){
  PatternMatchInfo newMatch;
  newMatch.PatternId = mPatternId;
  newMatch.PathId = mPathId;
  newMatch.MatchLikelihood = matchLklhd;
  fPatternsMatched.push_back(newMatch);
}


// *********************************************************************************
int ND::TTPCPath::GetNMatchedPattern(){
  return fPatternsMatched.size();
}


// *********************************************************************************
double ND::TTPCPath::GetPatternMatchLikelihood(int i){
  return fPatternsMatched[i].MatchLikelihood.Total;
}


// *********************************************************************************
int ND::TTPCPath::GetMatchPatternId(int i){
  return fPatternsMatched[i].PatternId;
}

// *********************************************************************************
int ND::TTPCPath::GetPatternMatchPathId(int i){
  return fPatternsMatched[i].PathId;
}

// *********************************************************************************
int ND::TTPCPath::GetMatchPatternIdIndex(int patternId){
  int idx = 0;
  for( std::vector<PatternMatchInfo>::iterator pmi = fPatternsMatched.begin(); pmi != fPatternsMatched.end(); pmi++){
    if ( (*pmi).PatternId == patternId)
      return idx;
    idx++;
  }
  return -1;
}



// *********************************************************************************
void ND::TTPCPath::SetEndClustersToNodes(){
  ND::THitSelection::const_iterator tmpClu = GetHits()->begin();
  ND::THandle<ND::TTPCHVCluster> Clu = *tmpClu;
  Clu->SetEndNode();
  ND::THitSelection::const_reverse_iterator tmpClu2 = GetHits()->rbegin();
  Clu = *tmpClu2;
  Clu->SetEndNode();
}



// *********************************************************************************
void ND::TTPCPath::SetPID(ND::TReconPID::ParticleId pid, double weight){
  fPID = pid;
  fPIDweight = weight;
}


// *********************************************************************************
void ND::TTPCPath::SaveInRealDatum(std::string name, double value){
  fInRealDatum[name] = value;
}


// *********************************************************************************
ND::THandle<ND::TReconBase> ND::TTPCPath::GetOAEventObj(){
  return fOAEventObj;
}


// *********************************************************************************
void ND::TTPCPath::PrepareTReconTrack( State HelixStart, State HelixEnd, ND::THandle<ND::TReconTrack> Track) {

  // TODO: Compare cluster list sizes as safety check

  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  ND::rpman("TREx").model_svc().model().intersector().set_length_sign(0);

  // This is the running state at each cluster
  double stepLength = 0.0;
  fLength = 0.0;
  fChi2 = 0.0;
  fNDOF = 0.0;
  double Chi2_dir[3];
  int NDOF_dir[3]; 
 
  for (int i = 0; i < 3; i++){
    Chi2_dir[i] = 0.0;
    NDOF_dir[i] = 0;
  }

  double Chi2X = 0.0;
  double Chi2Y = 0.0;
  double Chi2Z = 0.0;

  int NDOFX = 0;
  int NDOFY = 0;
  int NDOFZ = 0;

  State localState;
  // Reverse or not depending on the first and last state.
  // Go downstream, i.e. path state with positive Z direction.
  // Give priority to negative tracks for curving back tracks.
  if(HelixStart.name(RP::representation) != RP::pos_dir_curv)
    RP::rep().convert(HelixStart, RP::pos_dir_curv);
  if(HelixEnd.name(RP::representation) != RP::pos_dir_curv)
    RP::rep().convert(HelixEnd, RP::pos_dir_curv);

  double StartdZ = HelixStart.vector()[5];
  double EnddZ   = HelixEnd.vector()[5];
  double QoP     = HelixStart.vector()[6];
  ND::THandle<ND::THitSelection> TRExClusters;

  if ( ((StartdZ * EnddZ) > 0.0 && StartdZ > 0.0) ||
       ((StartdZ * EnddZ) < 0.0 && QoP < 0.0) ){
//  if(1){
    localState = HelixStart;
    TRExClusters = GetHits();
  } else {
    // Reverse the cluster ordering !!!
    TRExClusters = ND::THandle<ND::THitSelection>(new ND::THitSelection());
    for (ND::THitSelection::const_reverse_iterator Clu = GetHits()->rbegin(); Clu != GetHits()->rend(); Clu++) {
      TRExClusters->push_back(*Clu);
    }
    // Reverse also the end state !!!
    localState = HelixEnd;
    ND::tman().ReverseStateSenseAndCharge(localState);
  }

  fTrackType = kSTANDARDTRACK;
  if ( (StartdZ * EnddZ) < 0.0 )
    fTrackType = kCURVBACKTRACK;

  // The clusters are not automatically converted so we have to do it by hand.
  // Do it early so we can use the new clusters to fill the nodes. 
  ND::THandle<ND::THitSelection> oaEvtClusters = Track->GetHits();
  // Fill these once we have figured out the ordering
  oaEvtClusters->clear();
  for (ND::THitSelection::const_iterator tmpClu = TRExClusters->begin(); tmpClu != TRExClusters->end(); tmpClu++) {
    ND::THandle<ND::TTPCHVCluster> oClu = *tmpClu; 
    ND::THandle<ND::TComboHit> nClu = oClu->ConvertToOAEvent(); 
    oaEvtClusters->push_back(nClu);
  }

  // Get Track Node container
  ND::TReconNodeContainer& trkNodes =  Track->GetNodes();

  ND::THitSelection::const_iterator tmpTRExClu = TRExClusters->begin();
  ND::THitSelection::const_iterator tmpoaEvtClu = oaEvtClusters->begin();
  ND::THandle<ND::TTPCHVCluster> FirstTRExClu = *(TRExClusters->begin());
  ND::THandle<ND::TTPCHVCluster> LastTRExClu = *(TRExClusters->rbegin());
  for ( ; tmpTRExClu != TRExClusters->end(); tmpTRExClu++, tmpoaEvtClu++) {
    //------- Propagate the helix to each cluster --------
    ND::THandle<ND::TTPCHVCluster> CurrentTRExClu = *tmpTRExClu;
    // If propagation fails, let's try our luck on the next cluster.
    stepLength = 0.0;
    if (!TTPCRecPackUtils::FullPropagateToHVCluster(CurrentTRExClu, localState, stepLength))
      continue;
    // TODO: What to do if the first and or last nodes don't work ?

    fLength += stepLength;

    if (CurrentTRExClu->isOkForFit()){
      fChi2 += TTPCUtils::State2CluChi2(localState, CurrentTRExClu, Chi2_dir, NDOF_dir);
      fNDOF += 2.0;
  
      Chi2X += Chi2_dir[0];
      Chi2Y += Chi2_dir[1];
      Chi2Z += Chi2_dir[2];
  
      NDOFX += NDOF_dir[0];
      NDOFY += NDOF_dir[1];
      NDOFZ += NDOF_dir[2];
    }
    // We save only the end clusters as nodes
    if( !(CurrentTRExClu->IsEndNode()))
      continue;

    // Create the node
    ND::THandle<ND::TReconNode> recoNode(new ND::TReconNode);

    // Prepare the node object
    ND::THandle<ND::TReconCluster> recoObject(new ND::TReconCluster);
    ND::THandle<ND::TComboHit> combo = *tmpoaEvtClu;
    recoObject->FillFromHits("tpcnodeHits", combo->GetHits());

    // Prepare the node state
    double sigma_Q = 0.1; //place holder
    localState.set_hv("EDeposit",HyperVector( recoObject->GetEDeposit(), sigma_Q ));    // the EDeposit
    // Always save the T0 because we have a default set in TTPCT0Finder even if no T0 was found.
    localState.set_hv("time",HyperVector( GetT0() , 1. ));    // the track T0

    ND::THandle<ND::TReconState> recoState(new ND::TTrackState);
    ND::converter().State_to_TReconState(localState,(*recoState));

    // Save the node object
    ND::THandle<ND::TReconBase> recoBaseObject = recoObject;
    recoNode->SetObject(recoBaseObject);

    if ( CurrentTRExClu->IsVertical() ){
      recoState->SetFixed(3); // The z plane is fixed.
    } else {
      recoState->SetFixed(2); // The y plane is fixed.
    }

    // Save the node state
    ND::THandle<ND::TTrackState>  castState = recoState;
    recoNode->SetState(recoState);

    trkNodes.push_back(recoNode); 
  }

  fNDOF -= 7;
  NDOFX -= 2;
  NDOFY -= 3;
  NDOFZ -= 2;
   
  //Saving individual Chi-squares & dof in Data for output Tree

  SaveInRealDatum("TrkChi2X",Chi2X);
  SaveInRealDatum("TrkChi2Y",Chi2Y);
  SaveInRealDatum("TrkChi2Z",Chi2Z);
  
  SaveInRealDatum("TrkNDOFX",NDOFX);
  SaveInRealDatum("TrkNDOFY",NDOFY);
  SaveInRealDatum("TrkNDOFZ",NDOFZ);
  
  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);

  Track->SetQuality(fChi2);
  Track->SetNDOF(fNDOF);
  // Define the overall state of the track using the first node state.
  ND::THandle<ND::TTrackState> ttrackState = Track->GetState();
  if( trkNodes.size() ) {
    ND::THandle<ND::TTrackState> nodeState = trkNodes.front()->GetState(); 
    if( nodeState && ttrackState ) {
      *ttrackState = *nodeState;
    }
  }


}



// *********************************************************************************
ND::THandle<ND::TReconBase> ND::TTPCPath::ConvertToOAEvent() {

  // Only Gas Output paths that were not merged or 
  if (fOAEventObj && fId < TREXSTDOUTIDOFFSET)
    return fOAEventObj;

  ND::THandle<ND::TReconBase> Output;
  int NbVerticalClu = 0;
  int NbHorizontalClu = 0;
  int NbFittedVerticalClu = 0;
  int NbFittedHorizontalClu = 0;
  //////////////////////////////////////////////////////////////////////
  // Tracks contained in one drift volume
  if ( fTrackType != kCATHCROSSTRACK){
    ND::THandle<ND::TReconTrack> ttrack( new TReconTrack((*this)) ); 

    for (ND::THitSelection::const_iterator tmpClu = GetHits()->begin(); tmpClu != GetHits()->end(); tmpClu++) {
      ND::THandle<ND::TTPCHVCluster> oClu = *tmpClu; 
      if (oClu->IsVertical()){
        NbVerticalClu++;
        if (oClu->isOkForFit())
          NbFittedVerticalClu++;
      } else {
        NbHorizontalClu++;
        if (oClu->isOkForFit())
          NbFittedHorizontalClu++;
      }
    }

    // Create all the states, for the nodes and the track.
    if ( CheckStatus(ND::TReconBase::kChi2Fit) &&
       (! CheckStatus(ND::TReconBase::kLikelihoodFit))){
      // We have only the seeding results
      PrepareTReconTrack(fFrontSeedState, fBackSeedState, ttrack);
    } else if (CheckStatus(ND::TReconBase::kLikelihoodFit) ){
      // We have a likelihood fit (Not necessarily a reliable one !)
      PrepareTReconTrack(fFrontFitState, fBackFitState, ttrack);
    } else {
      // The track clusters are filled in PrepareTReconTrack when states are available.
      // SO do it here by hand when the seeding failed.
      ND::THandle<ND::THitSelection> newClusters = ttrack->GetHits();
      newClusters->clear();
      for (ND::THitSelection::const_iterator tmpClu = GetHits()->begin(); tmpClu != GetHits()->end(); tmpClu++) {
        ND::THandle<ND::TTPCHVCluster> oClu = *tmpClu; 
        ND::THandle<ND::TComboHit> nClu = oClu->ConvertToOAEvent(); 
        newClusters->push_back(nClu);
      }
    }
    // TODO: else, make one state at the first cluster with dummy direction and curvature
    // At least we'll have the position and time for the T0

    // With or without PID information ?
    Output = ttrack;
    if (fPID != ND::TReconPID::kNotSet){
      ND::THandle<ND::TReconPID> tpid(new ND::TReconPID());
      // Should we use the constructor taking TReconTrack here instead ???
      bool TODO_ConvertTrackToPID = ND::converter().TReconTrack_to_TReconPID(ttrack, tpid);
      tpid->SetParticleId(fPID);
      tpid->SetPIDWeight(fPIDweight);
      // TODO: Check that the conversion worked ... Exception if it didn't ?
      Output = tpid;
    }

    // Clear the TTPCPaths constituents which were copied in the constructor
    ND::THandle<ND::TReconObjectContainer> tmpConstituents = ttrack->GetConstituents();
    if (tmpConstituents)
      ttrack->GetConstituents()->clear();

    ND::THandle<ND::TReconObjectContainer> Constituents = GetConstituents();
    if (Constituents){
      for (ND::TReconObjectContainer::iterator ptc = Constituents->begin(); ptc != Constituents->end(); ptc++){
        ND::THandle<ND::TTPCPath> ConstPath = (*ptc);
        if (!ConstPath)
          continue;
        ND::THandle<ND::TReconBase> ConstOAEvt = ConstPath->GetOAEventObj();
        if (!ConstOAEvt)
          continue;
        // Add to ttrack and not output because tpid will have ttrack as constituent already.
        ttrack->AddConstituent(ConstOAEvt);
      }
    }
    // This path doesn't have constituents so it wasn't merged, 
    // but it has a fOAEventObj so we are currently preparing the StandardOutput
    else if (fOAEventObj){
      ttrack->AddConstituent(fOAEventObj);
    }
    // No constituent and no previous oaEvent conversion
    // so keep this for book keeping downstream after globalRecon
    else {
      fOAEventObj = Output;
    }
  }
  //////////////////////////////////////////////////////////////////////
  // Cathode crossers
  else {
    ND::THandle<ND::TReconObjectContainer> Constituents = GetConstituents();
    if (!Constituents){
      std::cerr<<"TTPCPath ERROR: no constituents in a cathode crosser !!!"<<std::endl;
      throw;
    }
    ND::THandle<ND::TTPCPath> ConstPath[2];
    ND::THandle<ND::TReconBase> ConstOAEvt[2];

    int NbFits = 0;
    int cstIdx = 0;
    for (ND::TReconObjectContainer::iterator ptc = Constituents->begin(); ptc != Constituents->end(); ptc++, cstIdx++){
      ConstPath[cstIdx] = *ptc;
      // Look at the segment and decide if they are "good enough", i.e. fitted and long enough
      ConstOAEvt[cstIdx] = ConstPath[cstIdx]->ConvertToOAEvent();
      ND::THandle<ND::TReconPID> tmpPID = ConstOAEvt[cstIdx];
      if(tmpPID && tmpPID->GetNodes().size() > 1){
        NbFits++;
      }
      for (ND::THitSelection::const_iterator tmpClu = ConstPath[cstIdx]->GetHits()->begin(); tmpClu != ConstPath[cstIdx]->GetHits()->end(); tmpClu++) {
        ND::THandle<ND::TTPCHVCluster> oClu = *tmpClu; 
        if (oClu->IsVertical()){
          NbVerticalClu++;
          if (oClu->isOkForFit())
            NbFittedVerticalClu++;
        } else {
          NbHorizontalClu++;
          if (oClu->isOkForFit())
            NbFittedHorizontalClu++;
        }
      }

    }
    if (NbFits > 1){
      // If both segments are good enough => KalmanFilter
      Output = TTPCUtils::MergeAndFitObjectsWithRecPack(ConstOAEvt[0], ConstOAEvt[1]);

    }

    if (!Output){
      // Only one segment as a good enough fit or the Kalman Filter fit failed
      // Create a new TReconPID of the good segment from scratch and simply add the clusters of the bad segment.
      // 1) Find the good (fitted) segment of the best of the two
      // First convert the constituents again to oaEvent containers because
      // the best of them will be used as output with the added clusters
      ND::THandle<ND::TReconBase> tmpA = ConstOAEvt[0];
      ND::THandle<ND::TReconBase> tmpB = ConstOAEvt[1];

      ND::THandle<ND::TReconBase> best;
      ND::THandle<ND::TReconBase> worst;
      ND::tman().GetBestObject( tmpA, tmpB, best, worst);
      ND::THandle<ND::TReconPID> bestPID = best;
      if( bestPID){
        Output = ND::THandle<ND::TReconPID> (new TReconPID(*bestPID));
      } else {
        ND::THandle<ND::TReconTrack> bestTrack = best;
        Output = ND::THandle<ND::TReconTrack> (new TReconTrack(*bestTrack));
      }

      // 2) Rearrange clusters properly
      // First copy the clusters of best which have the right orientation
      ND::THandle<ND::THitSelection> bestClusters  = best->GetHits();
      ND::THandle<ND::THitSelection> worstClusters = worst->GetHits();
      unsigned int bestEnd, worstEnd;
      TTPCUtils::FindClosestEnds(best, worst, bestEnd, worstEnd);

      ND::THitSelection* newClusters;
      if (bestEnd == 1 ){
        newClusters = new ND::THitSelection(*bestClusters);
        if (worstEnd == 0 ){
          for (ND::THitSelection::const_iterator hit = worst->GetHits()->begin(); hit != worst->GetHits()->end(); ++hit)
            newClusters->push_back(*hit);
        } else {
          for (ND::THitSelection::const_reverse_iterator hit = worst->GetHits()->rbegin(); hit != worst->GetHits()->rend(); ++hit)
            newClusters->push_back(*hit);
        }
      } else {
        if (worstEnd == 0 ){
          newClusters = new ND::THitSelection(*worstClusters);
        } else {
          newClusters = new ND::THitSelection();
          for (ND::THitSelection::const_reverse_iterator hit = worst->GetHits()->rbegin(); hit != worst->GetHits()->rend(); ++hit)
            newClusters->push_back(*hit);
        }
        for (ND::THitSelection::const_iterator hit = best->GetHits()->begin(); hit != best->GetHits()->end(); ++hit)
          newClusters->push_back(*hit);
      }
      bestClusters->clear();

      Output->AddHits(newClusters);


    }

    // Set proper constituents ConstOAEvt
    ND::THandle<ND::TReconObjectContainer> tmpConstituents = Output->GetConstituents();
    if (tmpConstituents)
      Output->GetConstituents()->clear();
    Output->AddConstituent(ConstOAEvt[0]);
    Output->AddConstituent(ConstOAEvt[1]);
  }
  //////////////////////////////////////////////////////////////////////

  Output->ClearStatus(ND::TReconBase::kRan);

  // Save all the TRealDatum and TIntegerDatum info needed downstream
  for (std::map<std::string,double>::iterator it=fInRealDatum.begin(); it!=fInRealDatum.end(); ++it){
    ND::TRealDatum* tmpRDatum = new ND::TRealDatum(it->first.c_str(), it->second);
    Output->AddDatum(tmpRDatum);
  }

  ND::TIntegerDatum* pathId = new ND::TIntegerDatum("PathId",fId);
  Output->AddDatum(pathId);

  ND::TIntegerDatum* patternId = new ND::TIntegerDatum("PatternId",fPatternId);
  Output->AddDatum(patternId);

  if (fJunctionId.size()){
    std::vector<unsigned int>::iterator JuIdIt = fJunctionId.begin();
    ND::TIntegerDatum* junctionIds = new ND::TIntegerDatum("JunctionIds", (*JuIdIt));
    JuIdIt++;
    for ( ; JuIdIt < fJunctionId.end(); JuIdIt++)
      junctionIds->push_back(*JuIdIt);
    Output->AddDatum(junctionIds);
  }
 
  ND::TIntegerDatum* TrackType = new ND::TIntegerDatum("TrackType", fTrackType);
  Output->AddDatum(TrackType);

  // Matching to the clusters of other paths,
  // so do that only if there is a successful fit.
  if (CheckStatus(ND::TReconBase::kLikelihoodFit) ){
    ND::TRealDatum* FitLklhd = new ND::TRealDatum("FitLogLikelihood", fFitResults.LogLikelihood.Total);
    Output->AddDatum(FitLklhd);
    ND::TRealDatum* FitLklhdX = new ND::TRealDatum("FitLogLikelihoodX", fFitResults.LogLikelihood.X);
    Output->AddDatum(FitLklhdX);
    ND::TRealDatum* FitLklhdHV = new ND::TRealDatum("FitLogLikelihoodHV", fFitResults.LogLikelihood.HV);
    Output->AddDatum(FitLklhdHV);

    if(fPathsMatched.size()){
      std::vector<PathMatchInfo>::iterator PaMaIt = fPathsMatched.begin();
      ND::TIntegerDatum* PathIdMatch = new ND::TIntegerDatum("PathIdMatch", (*PaMaIt).PathId);
      ND::TRealDatum* PathMatchLklhd   = new ND::TRealDatum("PathMatchLklhd", (*PaMaIt).MatchLikelihood.Total);
      ND::TRealDatum* PathMatchLklhdX  = new ND::TRealDatum("PathMatchLklhdX", (*PaMaIt).MatchLikelihood.X);
      ND::TRealDatum* PathMatchLklhdHV = new ND::TRealDatum("PathMatchLklhdHV", (*PaMaIt).MatchLikelihood.HV);
      PaMaIt++;
      for ( ; PaMaIt < fPathsMatched.end(); PaMaIt++){
        PathIdMatch->push_back((*PaMaIt).PathId);
        PathMatchLklhd->push_back((*PaMaIt).MatchLikelihood.Total);
        PathMatchLklhdX->push_back((*PaMaIt).MatchLikelihood.X);
        PathMatchLklhdHV->push_back((*PaMaIt).MatchLikelihood.HV);
      }
      Output->AddDatum(PathIdMatch);
      Output->AddDatum(PathMatchLklhd);
      Output->AddDatum(PathMatchLklhdX);
      Output->AddDatum(PathMatchLklhdHV);
    }
    if(fPatternsMatched.size()){
      std::vector<PatternMatchInfo>::iterator PaMaIt = fPatternsMatched.begin();
      ND::TIntegerDatum* PatternIdMatch     = new ND::TIntegerDatum("PatternIdMatch", (*PaMaIt).PatternId);
      ND::TIntegerDatum* PatternPathIdMatch = new ND::TIntegerDatum("PatternPathIdMatch", (*PaMaIt).PathId);
      ND::TRealDatum* PatternMatchLklhd     = new ND::TRealDatum("PatternMatchLklhd", (*PaMaIt).MatchLikelihood.Total);
      ND::TRealDatum* PatternMatchLklhdX    = new ND::TRealDatum("PatternMatchLklhdX", (*PaMaIt).MatchLikelihood.X);
      ND::TRealDatum* PatternMatchLklhdHV   = new ND::TRealDatum("PatternMatchLklhdHV", (*PaMaIt).MatchLikelihood.HV);
      PaMaIt++;
      for ( ; PaMaIt < fPatternsMatched.end(); PaMaIt++){
        PatternIdMatch->push_back((*PaMaIt).PatternId);
        PatternPathIdMatch->push_back((*PaMaIt).PathId);
        PatternMatchLklhd->push_back((*PaMaIt).MatchLikelihood.Total);
        PatternMatchLklhdX->push_back((*PaMaIt).MatchLikelihood.X);
        PatternMatchLklhdHV->push_back((*PaMaIt).MatchLikelihood.HV);
      }
      Output->AddDatum(PatternIdMatch);
      Output->AddDatum(PatternPathIdMatch);
      Output->AddDatum(PatternMatchLklhd);
      Output->AddDatum(PatternMatchLklhdX);
      Output->AddDatum(PatternMatchLklhdHV);
    }
  }


  ND::TRealDatum* length = new ND::TRealDatum("Length",fLength);
  Output->AddDatum(length);

  ND::TIntegerDatum* NbVerticalCluDatum = new ND::TIntegerDatum("NbVerticalClusters",NbVerticalClu);
  Output->AddDatum(NbVerticalCluDatum);

  ND::TIntegerDatum* NbHorizontalCluDatum = new ND::TIntegerDatum("NbHorizontalClusters",NbHorizontalClu);
  Output->AddDatum(NbHorizontalCluDatum);

  ND::TIntegerDatum* NbFittedVerticalCluDatum = new ND::TIntegerDatum("NbFittedVerticalClusters",NbFittedVerticalClu);
  Output->AddDatum(NbFittedVerticalCluDatum);

  ND::TIntegerDatum* NbFittedHorizontalCluDatum = new ND::TIntegerDatum("NbFittedHorizontalClusters",NbFittedHorizontalClu);
  Output->AddDatum(NbFittedHorizontalCluDatum);

  fT0.FillTRealData(Output);

  return Output;
}

