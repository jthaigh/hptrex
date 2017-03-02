#include "TTRExPath.hxx"

void trex::TTRExPath::AddConnectedJunction(trex::TTRExJunction* junct){

  if(find(fConnectedJunctions.begin(),fConnectedJunctions.end(), junct) != fConnectedJunctions.end()){
    std::cout << "This Junction is already connected" << std::endl;
    return;
  }
  if(fFrontIsConnected&&fBackIsConnected){
    std::cout<<"Attempted to add more than two junctions to path!!"<<std::endl;
    exit(0);
  }
  else if(fFrontIsConnected){
    fConnectedJunctions.push_back(junct);
    fConnectedJunctionsId.push_back(junct->GetId());
    fBackIsConnected=true;
  }
  else if(fBackIsConnected){
    trex::TTRExJunction* tmpJunct=fConnectedJunctions.front();
    fConnectedJunctions.clear();
    fConnectedJunctionsId.clear();
    fConnectedJunctions.push_back(junct);
    fConnectedJunctionsId.push_back(junct->GetId());
    fConnectedJunctions.push_back(tmpJunct);
    fConnectedJunctionsId.push_back(tmpJunct->GetId());
    fFrontIsConnected=true;
  }
  else{
    fConnectedJunctions.push_back(junct);
    fConnectedJunctionsId.push_back(junct->GetId());
    TVector3 frontPos=fClusters.front()->GetPosition();
    TVector3 backPos=fClusters.back()->GetPosition();
    TVector3 junctPos=junct->GetPosition();
    if((frontPos-junctPos).Mag()<(backPos-junctPos).Mag()){
      fFrontIsConnected=true;
    }
    else{
      fBackIsConnected=true;
    }
  }
}


void trex::TTRExPath::SetEndClustersToNodes(){  
  std::vector<trex::TTRExHVCluster*>::iterator tmpClu = GetClusters().begin();
  trex::TTRExHVCluster* Clu = *tmpClu;
  Clu->SetEndNode();
  std::vector<trex::TTRExHVCluster*>::reverse_iterator tmpClu2 = GetClusters().rbegin();
  Clu = *tmpClu2;
  Clu->SetEndNode();
}

//This method is overloaded  
void trex::TTRExPath::SaveFitState(std::vector<double> inState){
  TTPCPathFitResults tmpFitRes;
  tmpFitRes.FitState = inState;
  tmpFitRes.IsFitReliable = true;
  SaveFitState(tmpFitRes);
}


void trex::TTRExPath::SaveFitState(TTPCPathFitResults &result){
  
  std::vector<double> localState = result.FitState;
  
  bool effIsVertical=true;

  if(fabs(localState[4])>fabs(localState[5])){
    effIsVertical=false;
  }

  std::vector<trex::TTRExHVCluster*>::iterator tmpTRExClu = this->GetClusters().begin();
  
  bool frontStateDone = false;

  trex::TTPCHelixPropagator& hp=trex::helixPropagator();
  hp.InitHelixPosDirQoP(localState, effIsVertical);
  
  for ( ; tmpTRExClu != this->GetClusters().end(); tmpTRExClu++){

    trex::TTRExHVCluster* Clu = *tmpTRExClu;

    bool ok = hp.PropagateToHVCluster(*Clu);
    
    if (!ok){
      continue;
    }
    
    hp.GetHelixPosDirQoP(localState);
    
    if(!frontStateDone){
      frontStateDone = true;
      fFrontFitState = localState; 
    }
    
    fBackFitState = localState;
    

  }
  
  if (!frontStateDone){
    fHasFitState = false;
    return;
  }
  
  fFitState = result;
  fHasFitState = true;
}


//Alternative way of testing for fit state...might not need this
/*bool trex::TTRExPath::HasFitState(){
  return (this->HasRunFit() && this->HasLikelihoodFit());

}
*/

bool trex::TTRExPath::HasReliableFitState(){
  return (this->HasFitState() && fFitState.IsFitReliable);
}

std::vector<double> trex::TTRExPath::GetFitState(){
  if(this->HasFitState()){
    return fFrontFitState;
  }
  return std::vector<double>(0);
}

std::vector<double> trex::TTRExPath::GetFrontFitState(){
  if(this->HasFitState()){
    return fFrontFitState;
  }
  return std::vector<double>(0);
}

std::vector<double> trex::TTRExPath::GetBackFitState(){
  if(this->HasFitState()){
    return fBackFitState;
  }
  return std::vector<double>(0);
}

std::vector<double> trex::TTRExPath::GetFrontSeedState(){
  if(this->HasSeedState()){
    return fFrontSeedState;
  }
  return std::vector<double>(0);
}

std::vector<double> trex::TTRExPath::GetBackSeedState(){
  if(this->HasSeedState()){
    return fBackSeedState;
  }
  return std::vector<double>(0);
}

void trex::TTRExPath::SaveSeedStates(std::vector<double>& frontSeedState, std::vector<double>& backSeedState){
  SetHasChi2Fit(true);
  fFrontSeedState=frontSeedState;
  fBackSeedState=backSeedState;
}


int trex::TTRExPath::GetConnectedEnd(unsigned int JunctionId){
  if(fConnectedJunctionsId.size() ==1){
    if(fFrontIsConnected)
      return -1;
    else
      return 1;
  }
  if((*(fConnectedJunctionsId.begin()) == JunctionId))
    return -1;
  else if((*(fConnectedJunctionsId.rbegin()) ==JunctionId))
    return 1;
  else
    return 0;
}


void trex::TTRExPath::SaveMatchedPath(unsigned int id, trex::TTPCLogLikelihood likelihood){

  PathMatchInfo newMatch;
  newMatch.PathId = id;
  newMatch.MatchLikelihood = likelihood;
  fPathsMatched.push_back(newMatch);
}

void trex::TTRExPath::SaveMatchedPattern(unsigned int mPatternId,unsigned int mPathId, TTPCLogLikelihood matchLklhd){
  PatternMatchInfo newMatch;
  newMatch.PatternId = mPatternId;
  newMatch.PathId = mPathId;
  newMatch.MatchLikelihood = matchLklhd;
  fPatternsMatched.push_back(newMatch);
}

unsigned int trex::TTRExPath::GetNMatchedPattern(){
  
  return fPatternsMatched.size();
}

unsigned int trex::TTRExPath::GetNMatchedPath(){
  
  return fPathsMatched.size();
}


unsigned int trex::TTRExPath::GetMatchPatternId(unsigned int i){

  return fPatternsMatched[i].PatternId;
} 

unsigned int trex::TTRExPath::GetPatternMatchPathId(unsigned int i){
  
  return fPatternsMatched[i].PathId;
}


unsigned int trex::TTRExPath::NbEndsFreeToMatch(){
  unsigned int count = 0;
  for ( int i = 0; i < 2; i++)
    if (fEndFreeToMatch[i])
      count++;
  return count;
}

unsigned int trex::TTRExPath::GetMatchPathIdIndex(unsigned int pathId){
  int idx = 0;
  for( std::vector<PathMatchInfo>::iterator pmi = fPathsMatched.begin(); pmi != fPathsMatched.end(); pmi++){
    if ( (*pmi).PathId == pathId)
      return idx;
    idx++;
  }
  return -1;
}


double trex::TTRExPath::GetPathMatchLikelihood(unsigned int n){
  return fPathsMatched[n].MatchLikelihood.Total;
}



bool trex::TTRExPath::IsEndFreeToMatch(unsigned int end){
  // TODO: Proper exception needed                           
  if (end > 1)
    throw;
  return fEndFreeToMatch[end];
}

double trex::TTRExPath::GetPatternMatchLikelihood(unsigned int n){
  return fPatternsMatched[n].MatchLikelihood.Total;
}


double trex::TTRExPath::GetLogLikelihood(){
  return fFitState.LogLikelihood.Total;
}

void trex::TTRExPath::SetEndNotFreeToMatch(int end){
  // TODO: Proper exception needed                                                     
  if (end > 1)
    throw;
  fEndFreeToMatch[end] = false;
}








