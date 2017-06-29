#include "TTPCLikelihoodMerge.hxx"
#include "TTPCUtils.hxx"
#include <map>

//*****************************************************************************
trex::TTPCLikelihoodMerge::TTPCLikelihoodMerge(){
//*****************************************************************************

  fRunThroughGoingMerging = true;
  fRunBrokenTracksMerging = true;
  fThruGoDeltaLogLklhdCut = true;
  fBrkTrkDeltaLogLklhdCut = true;

  fThruGoDeltaLogLklhdCut = 200;
  fBrkTrkDeltaLogLklhdCut = 500;
}


//*****************************************************************************
void trex::TTPCLikelihoodMerge::Process(std::vector<trex::TTRExPattern>& inputPatterns, std::vector<trex::TTRExPattern>& mergedPatterns){
//*****************************************************************************

  for (auto patit = inputPatterns.begin(); patit != inputPatterns.end(); patit++) {
    trex::TTRExPattern& Pattern = *patit;
    Pattern.SetUsable(true);
  }

  // 1) Loop over the patterns
  for (auto basePatit = inputPatterns.begin(); basePatit != inputPatterns.end(); basePatit++) {
    trex::TTRExPattern& BasePattern = *basePatit;
    // Did we already use this pattern in another set of matchings ?
    if ( !BasePattern.IsUsable())
      continue;
    // if not, start a chain with this pattern.
    fPatternChain[0] = &BasePattern;
    BasePattern.SetUsable(false);
    fNbPatternChained = 1;
    MatchBrokenPaths(inputPatterns);
    
    // Now check all the patterns in the chain to see if they have matches through the junctions
    for (int pidx = 0; pidx < fNbPatternChained; pidx++){
      MatchThroughJunctions(*(fPatternChain[pidx]));
    }
    
    bool mergingNeeded = false;
    for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
      if (fMTracker[mt].NeedsMerging()){
	mergingNeeded = true;
	break;
      }
    }
    if (! mergingNeeded){
      // We have nothing to merge in this pattern, just copy it.
      mergedPatterns.push_back(BasePattern);
      // Clean up after each chain of matching.
      CleanUp();
      continue;
    }
    // All matches are found so now merge !
    MergeAll(mergedPatterns);
    
    // Clean up after each chain of matching.
    CleanUp();
  }
}


//*****************************************************************************
void trex::TTPCLikelihoodMerge::CleanUp(){
//*****************************************************************************
  fMTracker.clear();
  fNbPatternChained = 0;
  for( int i = 0; i < MAXPATTERNCHAIN; i++){
    fPatternChain[i] = 0;
  }
}


//*****************************************************************************
double trex::TTPCLikelihoodMerge::PathToPatternMatch( trex::TTRExPath& PathA, trex::TTRExPattern& PatternB, trex::TTRExPath* bestPathB){
  double bestDLL = 1.e13;
  unsigned int bestEndA = 2;
  unsigned int bestEndB = 2;
  for ( int m = 0; m < PathA.GetNMatchedPattern(); m++){
    if (int(PatternB.GetId()) == PathA.GetMatchPatternId(m)){
      // Get the corresponding matched PathB likelihood fit
      trex::TTRExPath* PathB=0;
      for (auto mconst = PatternB.GetPaths().begin(); mconst != PatternB.GetPaths().end(); mconst++) {
        trex::TTRExPath& tmpPath = *mconst;
        if ( int(tmpPath.GetId()) == PathA.GetPatternMatchPathId(m)){
          PathB = &tmpPath;
          // If there are no free ends on this path, stop here.
          if ( !PathB->NbEndsFreeToMatch()){
            break;
          }
          break;
        }
      }
      if ( !PathB){
        return bestDLL;
      }

      // Check that the closest ends of the paths are actually free for matching
      unsigned int UseEndA, UseEndB;

      TTPCUtils::FindClosestEnds(PathA, *PathB, UseEndA, UseEndB);
      if ( (!PathA.IsEndFreeToMatch(UseEndA)) || (!PathB->IsEndFreeToMatch(UseEndB)) )
        continue;

      // Calculate deltaLogLikelihood and apply cut
      double DeltaLogLikelihood = PathA.GetPatternMatchLikelihood(m) - PathB->GetLogLikelihood();
      if (DeltaLogLikelihood < fBrkTrkDeltaLogLklhdCut){
        if (DeltaLogLikelihood < bestDLL){
          bestDLL = DeltaLogLikelihood;
          bestPathB = PathB;
          bestEndA = UseEndA;
          bestEndB = UseEndB;
          break;
        }
      }
    }
  }

  // If we have a match, make sure that we won't reuse these track ends.
  if (bestPathB){
    PathA.SetEndNotFreeToMatch(bestEndA);
    bestPathB->SetEndNotFreeToMatch(bestEndB);
  }

  return bestDLL;
}



//*****************************************************************************
void trex::TTPCLikelihoodMerge::MatchBrokenPaths(std::vector< trex::TTRExPattern>& patterns){
//*****************************************************************************
  if ( !fRunBrokenTracksMerging)
    return;

  // We need to make sure that paths are not used more than twice 

  // We will loop on the patterns in the chain
  // until we can't find anymore matches to other patterns.
  for (int pidx = 0; pidx < fNbPatternChained && fNbPatternChained < MAXPATTERNCHAIN; pidx++){
    trex::TTRExPattern* PatternA = fPatternChain[pidx];
  
    for (auto patitB = patterns.begin(); patitB != patterns.end() && fNbPatternChained < MAXPATTERNCHAIN; patitB++) {
      trex::TTRExPattern* PatternB = &*patitB;
      if ( !PatternB->IsUsable())
        continue;

      double BestMatchDLL = 1.e13;
      trex::TTRExPath* BestMatchPaths[2]={0,0};

      // 1.a) Try to match A->B
      for (auto constit = PatternA->GetPaths().begin(); constit != PatternA->GetPaths().end(); constit++) {
        trex::TTRExPath& PathA = *constit;

	trex::TTRExPath* tmpPath=0;
        double tmpDLL = PathToPatternMatch(PathA, *PatternB, tmpPath);
        if (tmpPath && tmpDLL < BestMatchDLL){
          BestMatchDLL = tmpDLL;
          BestMatchPaths[0] = &PathA;
          BestMatchPaths[1] = tmpPath;
        } 

      }
      
      // 1.a) Try to match B->A only if A->B didn't yield any match
      for (auto constit = PatternB->GetPaths().begin(); constit != PatternB->GetPaths().end(); constit++) {
        trex::TTRExPath& PathB = *constit;

        trex::TTRExPath* tmpPath=0;
        double tmpDLL = PathToPatternMatch(PathB, *PatternA, tmpPath);
        if (tmpPath && tmpDLL < BestMatchDLL){
          BestMatchDLL = tmpDLL;
          BestMatchPaths[0] = tmpPath;
          BestMatchPaths[1] = &PathB;
        } 

      }
      // Did we find a match between these patterns ?
      if (BestMatchPaths[0] && BestMatchPaths[1]){
        fMTracker.push_back( MatchingTracker(*(BestMatchPaths[0]), *(BestMatchPaths[1])) );
        fPatternChain[fNbPatternChained] = PatternB;
        PatternB->SetUsable(false);
        fNbPatternChained++;
      }
    }
  }
}


//*****************************************************************************
void trex::TTPCLikelihoodMerge::MatchThroughJunctions(trex::TTRExPattern& pattern){
//*****************************************************************************
  if ( !fRunThroughGoingMerging)
    return;

  if (pattern.GetPaths().size() == 1)
    return;

  //bool atLeastOneMatch = false;
  ////// Look for a match across each junction
  for (auto constit = pattern.GetJunctions().begin(); constit != pattern.GetJunctions().end(); constit++) {
    trex::TTRExJunction& junction = (*constit);
    std::vector< trex::TTRExPath* > Paths;
    int NbPath = junction.GetConnectedPaths().size();
    for (auto pathit = junction.GetConnectedPaths().begin(); pathit != junction.GetConnectedPaths().end(); pathit++) {
      std::cout<<"NConnected="<<junction.GetConnectedPaths().size()<<std::endl;
      trex::TTRExPath* Path = *pathit;
      Paths.push_back(Path);
      }
    // 1) Create dynamic size LogLikelihood matrix NxN with N = number of paths connected to junction
    double** LklhdMatrix = new double*[NbPath];
    for ( int i = 0;  i < NbPath; ++i){
      // Build rows
      LklhdMatrix[i] = new double[NbPath];
      // and initialize them
      for ( int j = 0;  j < NbPath; ++j)
        LklhdMatrix[i][j] = 0xABCDEF;
    }
    /*
    The matrix is built to have:
              clusters
             1  2  3  4
           1
    states 2
           3
           4
    Therefore L13, the likelihood for the propagation
    of the state 1 onto the clusters 3 will be row 1, column 3.
    */

    double trace = 0.0;
    for ( int row = 0; row < NbPath; row++){
      trex::TTRExPath* pathA = Paths[row];
      if ( ! pathA->HasFitState()){
        // No reliable fit state, fill big values to be sure
        // the corresponding matching is discarded
        for ( int col = 0; col < NbPath; col++)
          LklhdMatrix[row][col] = 0xABCDEF;
        LklhdMatrix[row][row] = 0.0;
        continue;
      }
      LklhdMatrix[row][row] = pathA->GetLogLikelihood();
      trace += LklhdMatrix[row][row];
      for ( int col = 0; col < NbPath; col++){
        if (row == col)
          continue;
        trex::TTRExPath* pathB = Paths[col];
        int matchIdx = pathA->GetMatchPathIdIndex(pathB->GetId());
        // If there is a match likelihood, save it.
        if (matchIdx > -1)
          LklhdMatrix[row][col] = pathA->GetPathMatchLikelihood(matchIdx);
      }
    }

    // 2) Loop over Match hypotheses and select best
    double bestLogLklhd = 99999999999.9;
    int bestCol = -1;
    int bestRow = -1;
    for ( int col = 0; col < NbPath; col++){
      for ( int row = 0; row < NbPath; row++){
        if(LklhdMatrix[row][col] == 0xABCDEF)
          continue;
        if(col == row)
          continue;
        double thisLklhd = trace - LklhdMatrix[col][col] + LklhdMatrix[row][col];
        if (thisLklhd < bestLogLklhd){
          bestLogLklhd = thisLklhd;
          bestCol = col;
          bestRow = row;
        }
      }
    }
      
    // bool MatchFound = true;
    bool MatchFound = false;
    // 3) Apply cut on delta LogLikelihood if we have a best match
    if ( bestCol > -1 && bestRow > -1){
      double DeltaLogLklhd = bestLogLklhd - trace;
      if (DeltaLogLklhd < fThruGoDeltaLogLklhdCut){
        MatchFound = true;
      }
    }

    // 4) Create new Match object add to list of matches
    // even if there isn't a merge across the junction so
    // we can keep track of all the junctions.
    if (MatchFound){
      fMTracker.push_back( MatchingTracker(*(Paths[bestRow]), *(Paths[bestCol]), junction) );
      //atLeastOneMatch = true;
    } else {
      fMTracker.push_back( MatchingTracker(junction) );
    }

    // Delete properly the matrix
    for( int i = 0; i < NbPath; ++i)
      delete LklhdMatrix[i];
    delete LklhdMatrix;
    
  }

}


//*****************************************************************************
void trex::TTPCLikelihoodMerge::MergeAll(std::vector<trex::TTRExPattern>& outputVector){
//*****************************************************************************

  /////// New pattern !!!!!!
  outputVector.emplace_back();
  trex::TTRExPattern* NewPattern = &(outputVector.back());
  //TTPCT0 T0(ND::tpcCalibration().GetDefaultT0());
  
  std::vector< trex::TTRExPath* > MergedPaths;
  std::map< trex::TTRExPath*,unsigned int > pathIndexMap;

  // 1) Create the new merged paths by creating a chain of paths that need to be merged.
  int Chain = 0;
  for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
    if( ! fMTracker[mt].NeedsMerging() )
      continue;
    // 1.a) First list the fMTracker entries that will be merged

    // Defines the two paths at the end of the chain of paths
    // which can be of any length.
    trex::TTRExPath* NewPathEnds[2];
    for ( unsigned int i = 0; i < 2; i++){
      NewPathEnds[i] = fMTracker[mt].GetRawPath(i);
    }
    fMTracker[mt].SetMergingChain(Chain);

    // Find the other matches involving the two paths just merged
    for ( int submt = 0; submt < int(fMTracker.size()); submt++){
      if( ! fMTracker[submt].NeedsMerging() )
        continue;
      for ( unsigned int i = 0; i < 2; i++){
        if (fMTracker[submt].HasThisPath(*(NewPathEnds[i]))){
          // Find which of the two paths must be added to find the new end of the chain.
          for (unsigned int j = 0; j < 2; j++){
            trex::TTRExPath* PathB = fMTracker[submt].GetRawPath(j);
            if (PathB == NewPathEnds[i])
              continue;
            NewPathEnds[i] = PathB;
            break;
          }
          fMTracker[submt].SetMergingChain(Chain);
          // Must loop over all the fMTracker to avoid skipping one too early in long chains.
          submt = -1;
          break;
        }
      }
    }

    // 1.b) Create a NewPath each time you add a path from another fMTracker entry. Only the last one matters.
    trex::TTRExPath* NewPath=0;
    for ( unsigned int submt = 0; submt < fMTracker.size(); submt++){
      if ( ! fMTracker[submt].IsMergingChainOk(Chain))
        continue;
      // First use both paths in the fMTracker entry
      if (!NewPath){

	//MDH TODO: This call creates a new path! Need to manage persistence...
	std::cout<<fMTracker[submt].GetRawPath(0)<<" "<<fMTracker[submt].GetRawPath(1)<<std::endl;
        NewPath = TTPCUtils::MergePaths(*(fMTracker[submt].GetRawPath(0)), *(fMTracker[submt].GetRawPath(1)));

        for ( unsigned int i = 0; i < 2; i++){
          MergedPaths.push_back(fMTracker[submt].GetRawPath(i));
        }
        continue;
      }
      
      // Here we need to find which of the two paths must be merged
      for (unsigned int j = 0; j < 2; j++){
        trex::TTRExPath* PathB = fMTracker[submt].GetRawPath(j);
        std::vector< trex::TTRExPath* >::iterator it = find (MergedPaths.begin(), MergedPaths.end(), PathB);
        if (it != MergedPaths.end())
          continue;
        int tmpId = NewPath->GetId();
        trex::TTRExPath* prevPath = NewPath;

	//MDH TODO: See above
        NewPath = TTPCUtils::MergePaths(*NewPath, *PathB);

        // Now this new path is at one end of NewPath
        MergedPaths.push_back(PathB);
      }

    }
    for ( unsigned int submt = 0; submt < fMTracker.size(); submt++){
      if ( ! fMTracker[submt].IsMergingChainOk(Chain))
        continue;
      fMTracker[submt].SetMergedPath(*NewPath);
    }

    Chain++;
  }

  // 2) Recreate/copy the unmerged paths
  // Use map to keep match between old and new paths
  std::map< trex::TTRExPath*,trex::TTRExPath* > PathCopy;
  NewPattern->GetPaths().reserve(1000);
  for (int pidx = 0; pidx < fNbPatternChained; pidx++){
    for (auto constit = fPatternChain[pidx]->GetPaths().begin(); constit != fPatternChain[pidx]->GetPaths().end(); constit++) {
      trex::TTRExPath* tmpPath = &*constit;

      // Only copy the unmerged ones
      auto it = find (MergedPaths.begin(), MergedPaths.end(), tmpPath);
      if (it == MergedPaths.end()){
	NewPattern->GetPaths().emplace_back(*tmpPath);
	trex::TTRExPath* newPath=&(NewPattern->GetPaths().back());
        //newPath->AddConstituent(tmpPath);
	//MDH TODO: Do we need to set ID properly?
	newPath->SetId(0);
	//        newPath->SetId(ND::tpcCalibration().GetPathId());
        //newPath->ClearJunctionIds();
        PathCopy[tmpPath] = newPath;
	pathIndexMap[newPath]=NewPattern->GetPaths().size()-1;
      }
    }
  }

  // 3) Recreate all the junctions and copy the unmerged paths and the new merged paths correctly.
  int NbNewJunctions = 0;
  for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
    // Skip matching between two segments of a broken which doesn't have a junction.
    trex::TTRExJunction* tmpJunc= fMTracker[mt].GetJunction();
    if( !tmpJunc){
      continue;
    }

    NbNewJunctions++;
    NewPattern->GetJunctions().emplace_back();
    NewPattern->GetMap().emplace_back();

    trex::TTRExJunction* NewJunction=&(NewPattern->GetJunctions().back());
    trex::TTRExJunction* OldJunction = fMTracker[mt].GetJunction();

    NewJunction->SetId(0);
    //MDH TODO: See above
    //    NewJunction->SetId(ND::tpcCalibration().GetJunctionId());
    // 3.a) Pass the hits
    NewJunction->SetHits(OldJunction->GetHits());
    // 3.b) Copy the paths not merged so that we don't modify the version saved in the GasInteractionOutput
    for (auto pathtit = OldJunction->GetConnectedPaths().begin(); pathtit != OldJunction->GetConnectedPaths().end(); pathtit++) {
      trex::TTRExPath* tmpPath = (*pathtit);
        
      std::vector< trex::TTRExPath* >::iterator it = find (MergedPaths.begin(), MergedPaths.end(), tmpPath);
      if (it == MergedPaths.end()){
        NewJunction->AddConnectedPath(PathCopy[tmpPath]);
	NewPattern->GetMap().back().push_back(pathIndexMap[PathCopy[tmpPath]]);	
      }
    }

    // 3.b) Is this junction connected to a new merged path ?
    std::vector< trex::TTRExPath*> AlreadyAdded;
    for ( unsigned int submt = 0; submt < fMTracker.size(); submt++){
      for (auto constit = OldJunction->GetConnectedPaths().begin(); constit != OldJunction->GetConnectedPaths().end(); constit++) {
	trex::TTRExPath* tmpPath = *constit;
        if (fMTracker[submt].HasThisPath(**constit)){
          auto it = find (AlreadyAdded.begin(), AlreadyAdded.end(), fMTracker[submt].GetMergedPath());
          if (it == AlreadyAdded.end()){
	    if(PathCopy.count(fMTracker[submt].GetMergedPath())==0){
	      NewPattern->GetPaths().emplace_back(*(fMTracker[submt].GetMergedPath()));
	      PathCopy[fMTracker[submt].GetMergedPath()]=&(NewPattern->GetPaths().back());
	      pathIndexMap[&(NewPattern->GetPaths().back())]=NewPattern->GetPaths().size()-1;
	    }
	    NewJunction->AddConnectedPath(PathCopy[fMTracker[submt].GetMergedPath()]);
	    AlreadyAdded.push_back(fMTracker[submt].GetMergedPath());
	  }
	}
      }
    }
  }

  // If there are no new junctions, we need to create another pattern for a clean track
  if (!NbNewJunctions) {
    TTRExPath pathCopy=*(fMTracker[0].GetMergedPath());
    outputVector.pop_back();
    outputVector.emplace_back();
    NewPattern = &(outputVector.back());
    NewPattern->GetPaths().push_back(pathCopy);
  }
  
  //MDH TODO: Do we need to set the pattern ID to something nontrivial?
  //  NewPattern->SetId(ND::tpcCalibration().GetPatternId());
  NewPattern->SetId(0);
  //  NewPattern->InitialSetup();
  //NewPattern->SetT0(0);

}


