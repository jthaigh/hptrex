#include "TTPCLikelihoodMatch.hxx"

#include "TTPCJunction.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCDebug.hxx"
#include "TTPCRecPackUtils.hxx"
#include "TTPCUtils.hxx"

#include <TGeomInfo.hxx>


//*****************************************************************************
trex::TTPCLikelihoodMatch::TTPCLikelihoodMatch( ){
//*****************************************************************************

  fLklhdCalc = new TTPCLikFitPath(true);
}

//*****************************************************************************
trex::TTPCLikelihoodMatch::~TTPCLikelihoodMatch( ){
//*****************************************************************************
  delete fLklhdCalc;
}


//*****************************************************************************
void trex::TTPCLikelihoodMatch::Process( std::vector<trex::TTRExPattern>& allPatterns){
//*****************************************************************************

  for (auto patit allPatterns.begin(); patit != allPatterns.end(); patit++) {
    trex::TTRExPattern&=*patit;
    if ( Pattern.GetConstituents()->size() != 1)
      MatchAcrossJunctions(Pattern);

    for (auto constit = Pattern.GetConstituents()->begin(); constit != Pattern.GetConstituents()->end(); constit++) {
      path->SetEndClustersToNodes();
    }
  }

  MatchBrokenPaths(allPatterns);
}


//*****************************************************************************
void trex::TTPCLikelihoodMatch::MatchAcrossJunctions(trex::TTRExPattern& Pattern){
//*****************************************************************************

  // Loop through the paths
  std::vector< std::vector<trex::TTPCHitPad> >& Junctions = Pattern.GetJunctions();
  std::vector<trex::TTRExPath> Paths = Pattern.GetPaths();
  for (std::vector< auto jct = Junctions.begin(); jct != Junctions.end(); jct++) {
      std::vector<trex::TTPCHitPad>& junction = *jct;
    
      //MDH TODO: Write new code to get paths connected to this junction
      // Extract the connected paths
      std::vector< trex::TTRExPath* > ConnectedPaths;
    for (std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin(); pth != Paths.end(); pth++) {
      ND::THandle<ND::TTPCPath> path = *pth;
      if(junction->IsPathConnected(path->GetId()))
        ConnectedPaths.push_back(path);
    }

    std::vector< trex::TTRExPath* >::iterator coPthA;
    std::vector< trex::TTRExPath* >::iterator coPthB;
    for (coPthA = ConnectedPaths.begin(); coPthA != ConnectedPaths.end(); coPthA++) {
      trex::TTRExPath& pathA = **coPthA;
      coPthB = coPthA + 1;
      for (; coPthB != ConnectedPaths.end(); coPthB++) {
	trex::TTRExPath& pathB = *coPthB;
        // Propagate A to B
        if (pathA.HasFitState()){
          MatchPathsAtJunction(pathA, pathB, junction.GetId());
        }
        // Propagate B to A
        if (pathB.HasFitState()){
          MatchPathsAtJunction(pathB, pathA, junction.GetId());
        }
      } 
    }
  } 
}



//*****************************************************************************
  void trex::TTPCLikelihoodMatch::MatchPathsAtJunction(trex::TTRExPath& Path1, trex::TTRExPath& Path2, int JunctionId){
//*****************************************************************************
  TTPCLogLikelihood Likelihood;
  Likelihood.Total = 0.0;
  Likelihood.X = 0.0;
  Likelihood.HV = 0.0;

  std::vector<double> propagState;

  //MDH TODO: Fix this to get junction connections properly
  // Front state is connected to this junction
  if ( Path1.GetConnectedEnd(JunctionId) == -1){
    propagState = Path1->GetFrontFitState();
    // We want to propagate towards the "outside" of the track.
    ND::tman().ReverseStateSenseAndCharge(propagState);
  // Back state is connected to this junction
  } else if ( Path1->GetConnectedEnd(JunctionId) == 1){
    propagState = Path1->GetBackFitState();
  } else {
    // PROBLEM. Don't do anything.
    return;
  }

  trex::TTRExHVCluster* targetClu;
  if ( Path2.GetConnectedEnd(JunctionId) == -1) {
    targetClu = *(Path2.GetClusterHits()->begin());
  } else if ( Path2->GetConnectedEnd(JunctionId) == 1) {
    targetClu = *(Path2->GetHits()->rbegin());
  } else {
    // PROBLEM. Don't do anything.
    return;
  }

  int targetSense = TTPCUtils::SenseFromTwoClusters(Path2, targetClu);

  if( CheckClusterMatch(propagState, targetClu, targetSense, true)) {
    ND::THandle<ND::THitSelection> orderedClusters;
    if ( Path2->GetConnectedEnd(JunctionId) == 1) {
      orderedClusters = ND::THandle<ND::THitSelection>(new ND::THitSelection());
      for (ND::THitSelection::const_reverse_iterator Clu = Path2->GetHits()->rbegin(); Clu != Path2->GetHits()->rend(); Clu++) {
        orderedClusters->push_back(*Clu);
      }
    } else {
      orderedClusters = Path2->GetHits();
    }
    Likelihood = GetMatchLikelihood(propagState, orderedClusters, false);
    if (!Likelihood.Total){
      Likelihood.Total = 2e14;
      Likelihood.X = 1e14;
      Likelihood.HV = 1e14;
    }
  } else {
    Likelihood.Total = 2e13;
    Likelihood.X = 1e13;
    Likelihood.HV = 1e13;
  }
  if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
    std::cout<<"      --> Likelihood match = "<<Likelihood.Total<<std::endl;
  Path1->SaveMatchedPath(Path2->GetId(), Likelihood);

}

//*****************************************************************************
void ND::TTPCLikelihoodMatch::MatchBrokenPaths(std::vector< ND::THandle<ND::TTPCPattern> > inPatterns){
//*****************************************************************************
  if ( ND::tpcDebug().LikelihoodMatch(DB_INFO))
    std::cout<<" ----- MatchBrokenPaths "<<std::endl;
  for (std::vector< ND::THandle<ND::TTPCPattern> >::iterator patitA = inPatterns.begin(); patitA != inPatterns.end(); patitA++) {
    ND::THandle<ND::TTPCPattern> PatternA = *patitA;
    if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
      std::cout<<"  -*- Pattern "<<PatternA->GetId()<<std::endl;
    // Loop over the paths in this pattern
    for (ND::TReconObjectContainer::iterator constitA = PatternA->GetConstituents()->begin(); constitA != PatternA->GetConstituents()->end(); constitA++) {
      ND::THandle<ND::TTPCPath> PathA = *constitA;
      if ( !PathA)
        continue;
      if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
        std::cout<<"    - Path "<<PathA->GetId()<<std::endl;
      if ( !PathA->HasFitState())
        continue;

      bool ForceSense = false;
      State StateA;
      // Find which ends of the pathA (if any) are connected to a junction
      // 3 possibilities:
      //   - 2 junctions => Don't try to match to another pattern
      //   - No junction => Propagate in both sense
      //   - 1 junction => Take the state from the free end and force sense
      if ( PathA->IsFrontConnected() && PathA->IsBackConnected() ){
        continue;
      } else if ( (!PathA->IsFrontConnected()) && (!PathA->IsBackConnected()) ){
        ForceSense = false;
        StateA = PathA->GetFrontFitState();
        if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
          std::cout<<"    > Use front state"<<std::endl;
      } else {
        ForceSense = true;
        if ( (!PathA->IsFrontConnected()) && PathA->IsBackConnected() ){
          StateA = PathA->GetFrontFitState();
          // We want to propagate towards the "outside" of the track.
          ND::tman().ReverseStateSenseAndCharge(StateA);
          if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
            std::cout<<"    > Use reversed front state"<<std::endl;
        } else {
          StateA = PathA->GetBackFitState();
          if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
            std::cout<<"    > Use back state"<<std::endl;
        }
      }

      // Loop over the other patterns -> patternB
      for (std::vector< ND::THandle<ND::TTPCPattern> >::iterator patitB = inPatterns.begin(); patitB != inPatterns.end(); patitB++) {
        ND::THandle<ND::TTPCPattern> PatternB = *patitB;
        if (PatternA == PatternB)
          continue;
        if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
          std::cout<<"    Matching to pattern "<<PatternB->GetId()<<std::endl;
        // Loop over the paths in patternB
        for (ND::TReconObjectContainer::iterator constitB = PatternB->GetConstituents()->begin(); constitB != PatternB->GetConstituents()->end(); constitB++) {
          ND::THandle<ND::TTPCPath> PathB = *constitB;
          if ( !PathB)
            continue;
          if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
            std::cout<<"      Matching to path "<<PathB->GetId()<<std::endl;
          // Find which ends of the pathB (if any) are connected to a junction
          // 3 possibilities:
          //   - 2 junctions => Don't try to match
          //   - No junction => Try to match to both ends
          //   - 1 junction => Use only one end
          bool matchedLastClu = false;
          State firstState = StateA;
          State lastState = StateA;
          if ( PathB->IsFrontConnected() && PathB->IsBackConnected() ){
            continue;
          } else if ( (!PathB->IsFrontConnected()) && (!PathB->IsBackConnected()) ){
            // 1) Try matching to the first B cluster
            ND::THandle<ND::TTPCHVCluster> firstClu = *(PathB->GetHits()->begin());
            double firstMatchDist = 999999.;
            int targetSense = TTPCUtils::SenseFromTwoClusters(PathB, firstClu);
            bool firstOk = CheckClusterMatch(firstState, firstClu, targetSense, ForceSense, firstMatchDist);
            // 2) Try matching to the last B cluster
            ND::THandle<ND::TTPCHVCluster> lastClu = *(PathB->GetHits()->rbegin());
            double lastMatchDist = 999999.;
            targetSense = TTPCUtils::SenseFromTwoClusters(PathB, lastClu);
            bool lastOk = CheckClusterMatch(lastState, lastClu, targetSense, ForceSense, lastMatchDist);
            // Which one is best ?
            if ( (!firstOk) && (!lastOk)) {
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        Failed basic match to end clusters"<<std::endl;
              continue;
            } else if ( firstOk && (!lastOk)) {
              // ===> Use first cluster !
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        State->cluster match to first cluster of path "<<PathB->GetId()<<" is successful"<<std::endl;
            } else if ( lastOk && (!firstOk)) {
              // ===> Use last cluster !
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        State->cluster match to last cluster of path "<<PathB->GetId()<<" is successful"<<std::endl;
              matchedLastClu = true;
            } else {
              matchedLastClu = firstMatchDist > lastMatchDist;
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE)){
                  std::cout<<"        State->cluster matches to first and last clusters are successful."<<std::endl;
                if ( matchedLastClu)
                  std::cout<<"        Based on distance, start with last cluster of path "<<PathB->GetId()<<std::endl;
                else
                  std::cout<<"        Based on distance, start with first cluster of path "<<PathB->GetId()<<std::endl;
              }
            }
          } else {
            ForceSense = true;
            if ( (!PathB->IsFrontConnected()) && PathB->IsBackConnected() ){
              // Try matching to the first B cluster
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        Only front end is available. Check match with first cluster of path "<<PathB->GetId()<<std::endl;
              ND::THandle<ND::TTPCHVCluster> firstClu = *(PathB->GetHits()->begin());
              int targetSense = TTPCUtils::SenseFromTwoClusters(PathB, firstClu);
              if( !CheckClusterMatch(firstState, firstClu, targetSense, ForceSense)){
                if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                  std::cout<<"        => Failure."<<std::endl;
                continue;
              }
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        => Success."<<std::endl;
            } else {
              // Try matching to the last B cluster
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        Only back end is available. Check match with first cluster of path "<<PathB->GetId()<<std::endl;
              ND::THandle<ND::TTPCHVCluster> lastClu = *(PathB->GetHits()->rbegin());
              int targetSense = TTPCUtils::SenseFromTwoClusters(PathB, lastClu);
              if( !CheckClusterMatch(lastState, lastClu, targetSense, ForceSense)){
                if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                  std::cout<<"        => Failure."<<std::endl;
                continue;
              }
              matchedLastClu = true;
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
                std::cout<<"        => Success."<<std::endl;
            }
          }
          // Reorder the clusters of PathB if necessary (if matched to last cluster)
          ND::THandle<ND::THitSelection> orderedClusters;
          State propagState;
          if ( matchedLastClu) {
            orderedClusters = ND::THandle<ND::THitSelection>(new ND::THitSelection());
            for (ND::THitSelection::const_reverse_iterator Clu = PathB->GetHits()->rbegin(); Clu != PathB->GetHits()->rend(); Clu++) {
              orderedClusters->push_back(*Clu);
            }
            propagState = lastState;
          } else {
            orderedClusters = PathB->GetHits();
            propagState = firstState;
          }
          // TODO: Should I worry about the direction/sense at the matched cluster ???
          // Calculate the likelihood 
          TTPCLogLikelihood Likelihood;
          Likelihood.Total = 0.0;
          Likelihood.X = 0.0;
          Likelihood.HV = 0.0;
          Likelihood = GetMatchLikelihood(propagState, orderedClusters, false);
          if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
            std::cout<<"      --> Likelihood match = "<<Likelihood.Total<<std::endl;
          if (Likelihood.Total){
            // std::cout<<" =======>>>> "<<Likelihood.Total<<std::endl;
            PathA->SaveMatchedPattern(PatternB->GetId(), PathB->GetId(), Likelihood);
          }
        }
      }
    }
  }
}
 

//*****************************************************************************
bool ND::TTPCLikelihoodMatch::CheckClusterMatch(State &propState, ND::THandle<ND::TTPCHVCluster> Target, int targetSense, bool ForcePropagSense){
//*****************************************************************************
  double matchDistance = 999999.;
  return CheckClusterMatch(propState, Target, targetSense, ForcePropagSense, matchDistance);
}

//*****************************************************************************
bool ND::TTPCLikelihoodMatch::CheckClusterMatch(State &propState, ND::THandle<ND::TTPCHVCluster> Target, int targetSense, bool ForcePropagSense, double &matchDistance){
//*****************************************************************************
  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  if (ForcePropagSense)
    ND::rpman("TREx").model_svc().model().intersector().set_length_sign(1);
  else
    ND::rpman("TREx").model_svc().model().intersector().set_length_sign(0);

  bool ok = TTPCRecPackUtils::FullPropagateToHVCluster(Target, propState);

  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);

  // The propagation failed !
  if (!ok){
    return false;
  }

  // SAFETY CUT. We shouldn't try to use this match if the
  // state is super far from the cluster's position !
  const TVector3 CluPos = Target->GetCalibPosition();
  TVector3 PropPos(propState.vector()[0], propState.vector()[1], propState.vector()[2]);
  matchDistance = (CluPos - PropPos).Mag();
  if ( matchDistance > 100.){
    return false;
  }

  // Check that the propagate state sense matches the crude sense given by the matched cluster and the next one in the track.
  // Very helpful to prevent matches of colinear tracks.
  if (targetSense != 0){
    if(propState.name(RP::representation) != RP::pos_dir_curv)
      RP::rep().convert(propState, RP::pos_dir_curv);
    if (  Target->IsVertical()  && (propState.vector()[5]*targetSense > 0))
      return true;
    else if ((!Target->IsVertical()) && (propState.vector()[4]*targetSense > 0))
      return true;
  }

  return false;
} 

//*****************************************************************************
TTPCLogLikelihood ND::TTPCLikelihoodMatch::GetMatchLikelihood(State trkState, ND::THandle<ND::THitSelection> cluToMatch, double trkLength){
//*****************************************************************************
  ND::THandle<ND::THitSelection> selectedClu(new THitSelection());
  for (ND::THitSelection::const_iterator tmpClu = cluToMatch->begin(); tmpClu != cluToMatch->end(); tmpClu++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
    if( !Cluster->isOkForFit() ) { continue;}  // Check that the plane is actually enabled.
    selectedClu->push_back(Cluster);
  }

  if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
    std::cout<<"        Calculate likelihood match to "<<selectedClu->size()<<" clusters"<<std::endl;
  if ( !fLklhdCalc->SetupLogLklhdCalculator(trkState, selectedClu, trkLength)){
    if ( ND::tpcDebug().LikelihoodMatch(DB_VERBOSE))
      std::cout<<"        SetupLogLklhdCalculator FAILED !"<<std::endl;
    TTPCLogLikelihood tmpLogLklhd;
    tmpLogLklhd.Total = 2e12;
    tmpLogLklhd.X = 1e12;
    tmpLogLklhd.HV = 1e12;
    return tmpLogLklhd;
  }

  // Calc likelihood
  TTPCLogLikelihood LogLikelihood = fLklhdCalc->LogLklhdCalculator();

  // Make sure to clear up the memory
  if( selectedClu->size() > 0 ) selectedClu->erase(selectedClu->begin(),selectedClu->end());

  return LogLikelihood;
}


