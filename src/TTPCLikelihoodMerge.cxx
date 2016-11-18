#include "TTPCLikelihoodMerge.hxx"
#include "TTPCUtils.hxx"

//*****************************************************************************
trex::TTPCLikelihoodMerge::TTPCLikelihoodMerge(){
//*****************************************************************************

  //MDH TODO: Implement this properly - which of these do we need?

  fRunThroughGoingMerging = true;
  fRunBrokenTracksMerging = true;
  fThruGoDeltaLogLklhdCut = true;
  fBrkTrkDeltaLogLklhdCut = true;


  /*
  fRunThroughGoingMerging = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunThroughGoingTracksMerging");
  fRunBrokenTracksMerging = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunBrokenTracksMerging");
  if (!fRunThroughGoingMerging)
    std::cout<<"TRexRecon WARNING: Merging of through going tracks disabled"<<std::endl;
  if (!fRunBrokenTracksMerging)
    std::cout<<"TRexRecon WARNING: Merging of broken tracks disabled"<<std::endl;

  fThruGoDeltaLogLklhdCut = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.LikelihoodMerge.ThroughGoing.DeltaLogLklhdCut");
  fBrkTrkDeltaLogLklhdCut = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.LikelihoodMerge.BrokenTracks.DeltaLogLklhdCut");
  */
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
      trex::TTRExPath* PathB;
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

      //MDH TODO: Implement this
      //TTPCUtils::FindClosestEnds(PathA, PathB, UseEndA, UseEndB);
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
      trex::TTRExPath* BestMatchPaths[2];

      // 1.a) Try to match A->B
      for (auto constit = PatternA->GetPaths().begin(); constit != PatternA->GetPaths().end(); constit++) {
        trex::TTRExPath& PathA = *constit;

	trex::TTRExPath* tmpPath;
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

        trex::TTRExPath* tmpPath;
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


  //MDH TODO: Implement this properly using path-to-junction map in pattern
  bool atLeastOneMatch = false;
  ////// Look for a match across each junction
  for (auto constit = pattern.GetJunctions().begin(); constit != pattern.GetJunctions().end(); constit++) {
    trex::TTRExJunction& junction = (*constit);
    std::vector< trex::TTRExPath* > Paths;
    int NbPath = 0;//junction.GetConstituents()->size();
    /*for (ND::TReconObjectContainer::iterator pathit = junction->GetConstituents()->begin(); pathit != junction->GetConstituents()->end(); pathit++) {
      ND::THandle<ND::TTPCPath> Path = *pathit;
      Paths.push_back(Path);
      }*/
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
      atLeastOneMatch = true;
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

//MDH TODO: Fixing this to make sure we get the right objects out and keep the old ones
//that are needed is a big job. Need to make sure our new objects have all the necessary
//bells and whistles...
/*
  std::vector< ND::TTPCPath* > MergedPaths;


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
        if (fMTracker[submt].HasThisPath(NewPathEnds[i])){
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

	//MDH TODO: Implement this
        NewPath = TTPCUtils::MergePaths(fMTracker[submt].GetRawPath(0), fMTracker[submt].GetRawPath(1));

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

	//MDH TODO: As above
        NewPath = TTPCUtils::MergePaths(NewPath, PathB);

        // Now this new path is at one end of NewPath
        MergedPaths.push_back(PathB);
      }

    }
    for ( unsigned int submt = 0; submt < fMTracker.size(); submt++){
      if ( ! fMTracker[submt].IsMergingChainOk(Chain))
        continue;
      fMTracker[submt].SetMergedPath(NewPath);
    }

    Chain++;
  }

  // 2) Recreate/copy the unmerged paths
  // Use map to keep match between old and new paths
  std::map< ND::THandle<ND::TTPCPath>,ND::THandle<ND::TTPCPath> > PathCopy;
  for (int pidx = 0; pidx < fNbPatternChained; pidx++){
    for (ND::TReconObjectContainer::iterator constit = fPatternChain[pidx]->GetConstituents()->begin(); constit != fPatternChain[pidx]->GetConstituents()->end(); constit++) {
      ND::THandle<ND::TTPCPath> tmpPath = (*constit);
      if (!tmpPath)
        continue;

      // Only copy the unmerged ones
      std::vector< ND::THandle<ND::TTPCPath> >::iterator it = find (MergedPaths.begin(), MergedPaths.end(), tmpPath);
      if (it == MergedPaths.end()){
        ND::THandle<ND::TTPCPath> newPath( new TTPCPath(*tmpPath));
        newPath->AddConstituent(tmpPath);
        newPath->SetId(ND::tpcCalibration().GetPathId());
        newPath->ClearJunctionIds();
        PathCopy[tmpPath] = newPath;
        if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
          std::cout<<"   + Old path Id "<<tmpPath->GetId()<<" copied into path Id "<<newPath->GetId()<<std::endl;

      }
    }
  }

  // 3) Recreate all the junctions and copy the unmerged paths and the new merged paths correctly.
  int NbNewJunctions = 0;
  for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
    // Skip matching between two segments of a broken which doesn't have a junction.
    ND::THandle<ND::TTPCJunction> tmpJunc= fMTracker[mt].GetJunction();
    if( !tmpJunc){
      continue;
    }
    ND::THandle<ND::TTPCJunction> OldJunction = fMTracker[mt].GetJunction();
    ND::THandle<ND::TTPCJunction> NewJunction( new ND::TTPCJunction(OldJunction->GetPosition()));
    NewJunction->SetId(ND::tpcCalibration().GetJunctionId());
    if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
      std::cout<<"   + NewJunction Id "<<NewJunction->GetId()<<" replacing junction Id "<<OldJunction->GetId()<<std::endl;
    // 3.a) Pass the hits
    ND::THandle<ND::THitSelection> oldHits = OldJunction->GetHits();
    ND::THitSelection *newHits = new ND::THitSelection();
    for (ND::THitSelection::const_iterator tmpHit = oldHits->begin(); tmpHit != oldHits->end(); tmpHit++)
      newHits->push_back(*tmpHit);
    NewJunction->AddHits(newHits);
    // 3.b) Copy the paths not merged so that we don't modify the version saved in the GasInteractionOutput
    for (ND::TReconObjectContainer::iterator pathtit = OldJunction->GetConstituents()->begin(); pathtit != OldJunction->GetConstituents()->end(); pathtit++) {
      ND::THandle<ND::TTPCPath> tmpPath = (*pathtit);
        
      std::vector< ND::THandle<ND::TTPCPath> >::iterator it = find (MergedPaths.begin(), MergedPaths.end(), tmpPath);
      if (it == MergedPaths.end()){
        NewJunction->AddConstituent(PathCopy[tmpPath]);
      }
    }

    // 3.b) Is this junction connected to a new merged path ?
    std::vector< ND::THandle<ND::TTPCPath> > AlreadyAdded;
    for ( unsigned int submt = 0; submt < fMTracker.size(); submt++){
      for (ND::TReconObjectContainer::iterator constit = OldJunction->GetConstituents()->begin(); constit != OldJunction->GetConstituents()->end(); constit++) {
        ND::THandle<ND::TTPCPath> tmpPath = *constit;
        if (fMTracker[submt].HasThisPath(*constit)){
          std::vector< ND::THandle<ND::TTPCPath> >::iterator it = find (AlreadyAdded.begin(), AlreadyAdded.end(), fMTracker[submt].GetMergedPath());
          if (it == AlreadyAdded.end()){
            NewJunction->AddConstituent(fMTracker[submt].GetMergedPath());
            AlreadyAdded.push_back(fMTracker[submt].GetMergedPath());
            if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
              std::cout<<"    + Add merged path Id "<<fMTracker[submt].GetMergedPath()->GetId()<<std::endl;
          }
        }
      }
    }

    NbNewJunctions++;
    NewPattern->AddJunction(NewJunction);
  }

  /////// New pattern !!!!!!
  outputVector.emplace_back();
  ND::THandle<ND::TTPCPattern> NewPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern() );
  TTPCT0 T0(ND::tpcCalibration().GetDefaultT0());

  // If there are new junctions, we need to create another pattern for a clean track
  if (!NbNewJunctions) {
    NewPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern(fMTracker[0].GetMergedPath()) );
  }
  
  NewPattern->SetId(ND::tpcCalibration().GetPatternId());
  NewPattern->InitialSetup();
  NewPattern->SetT0(T0);

*/
}


