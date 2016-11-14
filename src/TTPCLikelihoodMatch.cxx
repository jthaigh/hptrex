#include "TTPCLikelihoodMatch.hxx"

#include "TTPCHitPad.hxx"
#include "TTPCUtils.hxx"

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
  Path1->SaveMatchedPath(Path2->GetId(), Likelihood);

}

//*****************************************************************************
void trex::TTPCLikelihoodMatch::MatchBrokenPaths(std::vector< ND::TTRExPattern >& inPatterns){
//*****************************************************************************

  for (auto patitA = inPatterns.begin(); patitA != inPatterns.end(); patitA++) {
    trex::TTRExPattern& PatternA = *patitA;
    // Loop over the paths in this pattern
    for (auto constitA = PatternA.GetPaths().begin(); constitA != PatternA.GetPaths().end(); constitA++) {
      ND::TTRExPath& PathA = *constitA;
      if ( !PathA.HasFitState())
        continue;

      bool ForceSense = false;
      std::vector<double> StateA;
      // Find which ends of the pathA (if any) are connected to a junction
      // 3 possibilities:
      //   - 2 junctions => Don't try to match to another pattern
      //   - No junction => Propagate in both sense
      //   - 1 junction => Take the state from the free end and force sense
      if ( PathA.IsFrontConnected() && PathA.IsBackConnected() ){
        continue;
      } else if ( (!PathA.IsFrontConnected()) && (!PathA.IsBackConnected()) ){
        ForceSense = false;
        StateA = PathA.GetFrontFitState();
      } 
      else {
        ForceSense = true;
        if ( (!PathA.IsFrontConnected()) && PathA.IsBackConnected() ){
          StateA = PathA.GetFrontFitState();
          // We want to propagate towards the "outside" of the track.
          
	  //MDH TODO: Implement reversing state and charge ourselves
	  //ND::tman().ReverseStateSenseAndCharge(StateA);
        } else {
          StateA = PathA.GetBackFitState();
        }
      }

      // Loop over the other patterns -> patternB
      for (auto patitB = inPatterns.begin(); patitB != inPatterns.end(); patitB++) {
        trex::TTRExPattern& PatternB = *patitB;
        if (&PatternA == &PatternB)
          continue;
        // Loop over the paths in patternB
        for (auto constitB = PatternB.GetPaths().begin(); constitB != PatternB.GetPaths().end(); constitB++) {
          trex::TTRExPath& PathB = *constitB;
          // Find which ends of the pathB (if any) are connected to a junction
          // 3 possibilities:
          //   - 2 junctions => Don't try to match
          //   - No junction => Try to match to both ends
          //   - 1 junction => Use only one end
          bool matchedLastClu = false;
          State firstState = StateA;
          State lastState = StateA;
          if ( PathB.IsFrontConnected() && PathB.IsBackConnected() ){
            continue;
          } else if ( (!PathB.IsFrontConnected()) && (!PathB.IsBackConnected()) ){
            // 1) Try matching to the first B cluster
            trex::TTPCHVCluster& firstClu = *(PathB.GetHits().begin());
            double firstMatchDist = 999999.;

	    //MDH TODO: Implement this...
	    //int targetSense = TTPCUtils::SenseFromTwoClusters(PathB, firstClu);
            bool firstOk = CheckClusterMatch(firstState, firstClu, targetSense, ForceSense, firstMatchDist);
            // 2) Try matching to the last B cluster
            trex::TTRExHVCluster> lastClu = *(PathB.GetHits().rbegin());
            double lastMatchDist = 999999.;
            //targetSense = TTPCUtils::SenseFromTwoClusters(PathB, lastClu);
            bool lastOk = CheckClusterMatch(lastState, lastClu, targetSense, ForceSense, lastMatchDist);
            // Which one is best ?
            if ( (!firstOk) && (!lastOk)) {
              continue;
            } else if ( firstOk && (!lastOk)) {
              // ===> Use first cluster !
              if ( ND::tpcDebug().LikelihoodMatch(DB_VVERBOSE))
            } else if ( lastOk && (!firstOk)) {
              // ===> Use last cluster !
              matchedLastClu = true;
            } else {
              matchedLastClu = firstMatchDist > lastMatchDist;
              }
            }
          } else {
            ForceSense = true;
            if ( (!PathB.IsFrontConnected()) && PathB.IsBackConnected() ){
              // Try matching to the first B cluster
              trex::TTRExHVCluster& firstClu = *(PathB.GetHits().begin());
              //MDH TODO: As above
	      //int targetSense = TTPCUtils::SenseFromTwoClusters(PathB, firstClu);
              if( !CheckClusterMatch(firstState, firstClu, targetSense, ForceSense)){
                continue;
              }
            } else {
              // Try matching to the last B cluster
              trex::TTRExHVCluster& lastClu = *(PathB.GetHits().rbegin());
              //MDH TODO: As above
	      //int targetSense = TTPCUtils::SenseFromTwoClusters(PathB, lastClu);
              if( !CheckClusterMatch(lastState, lastClu, targetSense, ForceSense)){
                continue;
              }
              matchedLastClu = true;
            }
          }
          // Reorder the clusters of PathB if necessary (if matched to last cluster)
          std::vector<trex::TTRExHVCluster*> orderedClusters;
	  std::vector<double> propagState;
          if ( matchedLastClu) {
            for (auto Clu = PathB.GetHits().rbegin(); Clu != PathB.GetHits().rend(); Clu++) {
              orderedClusters.push_back(&Clu);
            }
            propagState = lastState;
          } else {
	    for (auto Clu = PathB.GetHits().begin(); Clu != PathB.GetHits().end(); Clu++) {
              orderedClusters.push_back(&Clu);            
	      propagState = firstState;
          }
          // TODO: Should I worry about the direction/sense at the matched cluster ???
          // Calculate the likelihood 
          TTPCLogLikelihood Likelihood;
          Likelihood.Total = 0.0;
          Likelihood.X = 0.0;
          Likelihood.HV = 0.0;
          Likelihood = GetMatchLikelihood(propagState, orderedClusters, false);
          if (Likelihood.Total){
            // std::cout<<" =======>>>> "<<Likelihood.Total<<std::endl;
            PathA.SaveMatchedPattern(PatternB.GetId(), PathB.GetId(), Likelihood);
          }
        }
      }
    }
  }
}
 

//*****************************************************************************
bool trex::TTPCLikelihoodMatch::CheckClusterMatch(std::vector<double> &propState, trex::TTRExHVCluster& Target, int targetSense, bool ForcePropagSense){
//*****************************************************************************
  double matchDistance = 999999.;
  return CheckClusterMatch(propState, Target, targetSense, ForcePropagSense, matchDistance);
}

//*****************************************************************************
bool trex::TTPCLikelihoodMatch::CheckClusterMatch(std::vector<double>& &propState, trex::TTRExHVCluster& Target, int targetSense, bool ForcePropagSense, double &matchDistance){
//*****************************************************************************

  trex::TTPCHelixPropagator& hp=trex::helixPropagator();
  hp.InitHelixPosDirQoP(propState,Target.IsVertical());
  double pthLength = 0.0;

  bool ok = hp.PropagateToHVCluster(Target, propState);

  // The propagation failed !
  if (!ok){
    return false;
  }

  // SAFETY CUT. We shouldn't try to use this match if the
  // state is super far from the cluster's position !
  const TVector3 CluPos = Target.GetCalibPosition();
  TVector3 PropPos(propState[0], propState[1], propState[2]);
  matchDistance = (CluPos - PropPos).Mag();
  if ( matchDistance > 100.){
    return false;
  }

  // Check that the propagate state sense matches the crude sense given by the matched cluster and the next one in the track.
  // Very helpful to prevent matches of colinear tracks.
  if (targetSense != 0){
    if (  Target.IsVertical()  && (propState[5]*targetSense > 0))
      return true;
    else if ((!Target.IsVertical()) && (propState[4]*targetSense > 0))
      return true;
  }

  return false;
}

//*****************************************************************************
TTPCLogLikelihood trex::TTPCLikelihoodMatch::GetMatchLikelihood(std::vector<double> trkState, std::vector<TTRExHVCluster>& cluToMatch, double trkLength){
//*****************************************************************************
  std::vector<trex::TTRExHVCluster> selectedClu;
  for (auto tmpClu = cluToMatch.begin(); tmpClu != cluToMatch.end(); tmpClu++) {
    ND::TTPCHVCluster Cluster = *tmpClu;
    if( !Cluster.isOkForFit() ) { continue;}  // Check that the plane is actually enabled.
    selectedClu.push_back(Cluster);
  }

  if ( !fLklhdCalc->SetupLogLklhdCalculator(trkState, selectedClu, trkLength)){
    TTPCLogLikelihood tmpLogLklhd;
    tmpLogLklhd.Total = 2e12;
    tmpLogLklhd.X = 1e12;
    tmpLogLklhd.HV = 1e12;
    return tmpLogLklhd;
  }

  // Calc likelihood
  TTPCLogLikelihood LogLikelihood = fLklhdCalc->LogLklhdCalculator();

  return LogLikelihood;
}


