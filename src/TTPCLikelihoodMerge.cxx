#include "TTPCLikelihoodMerge.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCUtils.hxx"
#include "TTPCDebug.hxx"

#include <TOARuntimeParameters.hxx>
#include <TGeomInfo.hxx>


//*****************************************************************************
ND::TTPCLikelihoodMerge::TTPCLikelihoodMerge(){
//*****************************************************************************
  fRunThroughGoingMerging = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunThroughGoingTracksMerging");
  fRunBrokenTracksMerging = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.RunBrokenTracksMerging");
  if (!fRunThroughGoingMerging)
    std::cout<<"TRexRecon WARNING: Merging of through going tracks disabled"<<std::endl;
  if (!fRunBrokenTracksMerging)
    std::cout<<"TRexRecon WARNING: Merging of broken tracks disabled"<<std::endl;

  fThruGoDeltaLogLklhdCut = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.LikelihoodMerge.ThroughGoing.DeltaLogLklhdCut");
  fBrkTrkDeltaLogLklhdCut = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.LikelihoodMerge.BrokenTracks.DeltaLogLklhdCut");

}


//*****************************************************************************
void ND::TTPCLikelihoodMerge::Process(ND::TReconObjectContainer *inputPatterns, ND::TReconObjectContainer *mergedPatterns){
//*****************************************************************************
  std::vector< ND::THandle<ND::TTPCPattern> > PatternPerDriftVol[6];
  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
    std::cout<<" ===== TTPCLikelihoodMerge"<<std::endl;
  for (ND::TReconObjectContainer::iterator patit = inputPatterns->begin(); patit != inputPatterns->end(); patit++) {
    ND::THandle<ND::TTPCPattern> Pattern = *patit;
    Pattern->SetUsable(true);
    for (ND::TReconObjectContainer::iterator constit = Pattern->GetConstituents()->begin(); constit != Pattern->GetConstituents()->end(); constit++) {
      ND::THandle<ND::TTPCPath> path = *constit;
      if ( !path)
        continue;
      ND::THandle<ND::TTPCHVCluster> cluster = *(path->GetHits()->begin());
      ND::THandle<ND::TTPCHitPad> hitPad = *(cluster->GetHits().begin());
      int tpc, half, mm, pad;
      ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();
      tpcGeom.GetGeometryInfo(hitPad->GetGeomId(), tpc, half, mm, pad);
      PatternPerDriftVol[(tpc*2)+half].push_back(Pattern);
      break;
    }
  }

  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE)){
    for( int i = 0; i < 6; i++){
      std::cout<<" -- Found "<<PatternPerDriftVol[i].size()<<" patterns in drift volume "<<i<<std::endl;
    }
  }

  for( int i = 0; i < 6; i++){
    // 1) Loop over the patterns
    for (std::vector< ND::THandle<ND::TTPCPattern> >::iterator basePatit = PatternPerDriftVol[i].begin(); basePatit != PatternPerDriftVol[i].end(); basePatit++) {
      ND::THandle<ND::TTPCPattern> BasePattern = *basePatit;
      // Did we already use this pattern in another set of matchings ?
      if ( !BasePattern->IsUsable())
        continue;
      // if not, start a chain with this pattern.
      if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
        std::cout<<" -----------------------------------------------------------------------"<<std::endl;
      fPatternChain[0] = BasePattern;
      BasePattern->SetUsable(false);
      fNbPatternChained = 1;
      MatchBrokenPaths(PatternPerDriftVol[i]);

      // Now check all the patterns in the chain to see if they have matches through the junctions
      for (int pidx = 0; pidx < fNbPatternChained; pidx++){
        MatchThroughJunctions(fPatternChain[pidx]);
      }

      if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE)){
        for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
          if( fMTracker[mt].NeedsMerging() ){
            if( fMTracker[mt].GetJunction()){
              std::cout<<"   => Merge path "<<fMTracker[mt].GetRawPath(0)->GetId()<<" with path "<<fMTracker[mt].GetRawPath(1)->GetId()<<" across junction "<<fMTracker[mt].GetJunction()->GetId()<<std::endl;
            } else {
              std::cout<<"   => Merge path "<<fMTracker[mt].GetRawPath(0)->GetId()<<" with path "<<fMTracker[mt].GetRawPath(1)->GetId()<<" from two different patterns"<<std::endl;
            }
          } else {
            std::cout<<"   => No merging needed across junction "<<fMTracker[mt].GetJunction()->GetId()<<std::endl;
          }
        }
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
        mergedPatterns->push_back(BasePattern);
        // Clean up after each chain of matching.
        CleanUp();
        continue;
      }
      // All matches are found so now merge !
      mergedPatterns->push_back(MergeAll());

      // Clean up after each chain of matching.
      CleanUp();
    }
  }
}


//*****************************************************************************
void ND::TTPCLikelihoodMerge::CleanUp(){
//*****************************************************************************
  fMTracker.clear();
  fNbPatternChained = 0;
  for( int i = 0; i < MAXPATTERNCHAIN; i++){
    fPatternChain[i] = ND::THandle<ND::TTPCPattern>();
  }
}


//*****************************************************************************
double ND::TTPCLikelihoodMerge::PathToPatternMatch( ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPattern> PatternB, ND::THandle<ND::TTPCPath> &bestPathB){
  double bestDLL = 1.e13;
  unsigned int bestEndA = 2;
  unsigned int bestEndB = 2;
  for ( int m = 0; m < PathA->GetNMatchedPattern(); m++){
    if (int(PatternB->GetId()) == PathA->GetMatchPatternId(m)){
      // Get the corresponding matched PathB likelihood fit
      ND::THandle<ND::TTPCPath> PathB;
      for (ND::TReconObjectContainer::iterator mconst = PatternB->GetConstituents()->begin(); mconst != PatternB->GetConstituents()->end(); mconst++) {
        ND::THandle<ND::TTPCPath> tmpPath = *mconst;
        if ( !tmpPath)
          continue;
        if ( int(tmpPath->GetId()) == PathA->GetPatternMatchPathId(m)){
          PathB = tmpPath;
          // If there are no free ends on this path, stop here.
          if ( !PathB->NbEndsFreeToMatch()){
            if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE) )
              std::cout<<"     SKIP the match from path "<<PathA->GetId()<<" to path "<<PathB->GetId()<<" in pattern "<<PatternB->GetId()<<" because path "<<PathB->GetId()<<" has connections at both ends already."<<std::endl;
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
      TTPCUtils::FindClosestEnds(PathA, PathB, UseEndA, UseEndB);
      if ( (!PathA->IsEndFreeToMatch(UseEndA)) || (!PathB->IsEndFreeToMatch(UseEndB)) )
        continue;

      // Calculate deltaLogLikelihood and apply cut
      double DeltaLogLikelihood = PathA->GetPatternMatchLikelihood(m) - PathB->GetLogLikelihood();
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
    PathA->SetEndNotFreeToMatch(bestEndA);
    bestPathB->SetEndNotFreeToMatch(bestEndB);
  }

  return bestDLL;
}



//*****************************************************************************
void ND::TTPCLikelihoodMerge::MatchBrokenPaths(std::vector< ND::THandle<ND::TTPCPattern> > patterns){
//*****************************************************************************
  if ( !fRunBrokenTracksMerging)
    return;

  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
    std::cout<<" ------ Match broken paths, "<<patterns.size()<<" patterns in input."<<std::endl;

  // We need to make sure that paths are not used more than twice 

  // We will loop on the patterns in the chain
  // until we can't find anymore matches to other patterns.
  for (int pidx = 0; pidx < fNbPatternChained && fNbPatternChained < MAXPATTERNCHAIN; pidx++){
    ND::THandle<ND::TTPCPattern> PatternA = fPatternChain[pidx];
    if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE) ){
      if ( pidx == 0)
        std::cout<<"  >>> Start pattern matching chain with pattern "<<PatternA->GetId()<<std::endl;
      else
        std::cout<<"  >>> Next pattern in pattern matching chain is pattern "<<PatternA->GetId()<<std::endl;
    }
  
    for (std::vector< ND::THandle<ND::TTPCPattern> >::iterator patitB = patterns.begin(); patitB != patterns.end() && fNbPatternChained < MAXPATTERNCHAIN; patitB++) {
      ND::THandle<ND::TTPCPattern> PatternB = *patitB;
      if ( !PatternB->IsUsable())
        continue;

      double BestMatchDLL = 1.e13;
      ND::THandle<ND::TTPCPath> BestMatchPaths[2];

      // 1.a) Try to match A->B
      for (ND::TReconObjectContainer::iterator constit = PatternA->GetConstituents()->begin(); constit != PatternA->GetConstituents()->end(); constit++) {
        ND::THandle<ND::TTPCPath> PathA = *constit;
        if ( !PathA)
          continue;

        ND::THandle<ND::TTPCPath> tmpPath;
        double tmpDLL = PathToPatternMatch(PathA, PatternB, tmpPath);
        if (tmpPath && tmpDLL < BestMatchDLL){
          if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
            std::cout<<"  -> pattern "<<PatternA->GetId()<<", path "<<PathA->GetId()<<" matched to pattern "<<PatternB->GetId()<<", path "<<tmpPath->GetId()<<" with DeltaLogLikelihood of "<<tmpDLL<<std::endl;
          BestMatchDLL = tmpDLL;
          BestMatchPaths[0] = PathA;
          BestMatchPaths[1] = tmpPath;
        } else {
          if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
            std::cout<<"  -> no match from pattern "<<PatternA->GetId()<<" path "<<PathA->GetId()<<" to pattern "<<PatternB->GetId()<<std::endl;
        }

      }
      
      // 1.a) Try to match B->A only if A->B didn't yield any match
      for (ND::TReconObjectContainer::iterator constit = PatternB->GetConstituents()->begin(); constit != PatternB->GetConstituents()->end(); constit++) {
        ND::THandle<ND::TTPCPath> PathB = *constit;
        if ( !PathB)
          continue;

        ND::THandle<ND::TTPCPath> tmpPath;
        double tmpDLL = PathToPatternMatch(PathB, PatternA, tmpPath);
        if (tmpPath && tmpDLL < BestMatchDLL){
          if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
            std::cout<<"  -> pattern "<<PatternB->GetId()<<", path "<<PathB->GetId()<<" matched to pattern "<<PatternA->GetId()<<", path "<<tmpPath->GetId()<<" with DeltaLogLikelihood of "<<tmpDLL<<std::endl;
          BestMatchDLL = tmpDLL;
          BestMatchPaths[0] = tmpPath;
          BestMatchPaths[1] = PathB;
        } else {
          if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
            std::cout<<"  -> no match from pattern "<<PatternB->GetId()<<" path "<<PathB->GetId()<<" to pattern "<<PatternA->GetId()<<std::endl;
        }

      }
      // Did we find a match between these patterns ?
      if (BestMatchPaths[0] && BestMatchPaths[1]){
        fMTracker.push_back( MatchingTracker(BestMatchPaths[0], BestMatchPaths[1]) );
        fPatternChain[fNbPatternChained] = PatternB;
        PatternB->SetUsable(false);
        fNbPatternChained++;
        if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
          std::cout<<"  >>> Add pattern "<<PatternB->GetId()<<" to the matching chain."<<std::endl;
      }
    }
  }

  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
    std::cout<<"  >>> Chain contains "<<fNbPatternChained<<" patterns."<<std::endl;
  if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE)){
    std::cout<<"      Pattern Ids: ";
    for (int pidx = 0; pidx < fNbPatternChained; pidx++){
      std::cout<<fPatternChain[pidx]->GetId()<<" ";
    }
    std::cout<<std::endl;
  }

}


//*****************************************************************************
void ND::TTPCLikelihoodMerge::MatchThroughJunctions(ND::THandle<ND::TTPCPattern> pattern){
//*****************************************************************************
  if ( !fRunThroughGoingMerging)
    return;

  if (pattern->GetConstituents()->size() == 1)
    return;

  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
    std::cout<<" ------ Match through junctions for pattern"<<pattern->GetId()<<std::endl;

  bool atLeastOneMatch = false;
  ////// Look for a match across each junction
  for (ND::TReconObjectContainer::iterator constit = pattern->GetConstituents()->begin(); constit != pattern->GetConstituents()->end(); constit++) {
    ND::THandle<ND::TTPCJunction> junction = (*constit);
    if (!junction)
      continue;
    int NbPath = junction->GetConstituents()->size();
    std::vector< ND::THandle<ND::TTPCPath> > Paths;
    for (ND::TReconObjectContainer::iterator pathit = junction->GetConstituents()->begin(); pathit != junction->GetConstituents()->end(); pathit++) {
      ND::THandle<ND::TTPCPath> Path = *pathit;
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
      ND::THandle<ND::TTPCPath> pathA = Paths[row];
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
        ND::THandle<ND::TTPCPath> pathB = Paths[col];
        int matchIdx = pathA->GetMatchPathIdIndex(pathB->GetId());
        // If there is a match likelihood, save it.
        if (matchIdx > -1)
          LklhdMatrix[row][col] = pathA->GetPathMatchLikelihood(matchIdx);
      }
    }

    if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE)){
      std::cout<<"   LogLikelihood matrix for junction "<<junction->GetId()<<std::endl;
      for ( int row = 0; row < NbPath; row++){
        for ( int col = 0; col < NbPath; col++){
          std::cout<<"     "<<LklhdMatrix[row][col];
        }
        std::cout<<std::endl;
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
        if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE)){
          std::cout<<" Through going hypothesis LogLikelihood = "<<trace<<" - "<<LklhdMatrix[col][col]<<" + "<<LklhdMatrix[row][col]<<" = "<<thisLklhd<<std::endl;
        }
        if (thisLklhd < bestLogLklhd){
          if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE)){
            std::cout<<" --> Current best match: row "<<row<<"  col "<<col<<std::endl;
          }
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
      if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE)){
        std::cout<<"   => Best DeltaLogLikelihood = "<<bestLogLklhd<<" - "<<trace<<" = "<<DeltaLogLklhd<<std::endl;
      }
      if (DeltaLogLklhd < fThruGoDeltaLogLklhdCut){
        MatchFound = true;
        if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
          std::cout<<"   => DeltaLogLikelihood SELECTED !"<<std::endl;
      }
    }

    // 4) Create new Match object add to list of matches
    // even if there isn't a merge across the junction so
    // we can keep track of all the junctions.
    if (MatchFound){
      fMTracker.push_back( MatchingTracker(Paths[bestRow], Paths[bestCol], junction) );
      atLeastOneMatch = true;
    } else {
      fMTracker.push_back( MatchingTracker(junction) );
    }

    // Delete properly the matrix
    for( int i = 0; i < NbPath; ++i)
      delete LklhdMatrix[i];
    delete LklhdMatrix;
    
  }

  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
    std::cout<<"   MatchingTracker size == "<<fMTracker.size()<<" after through-junction matching of pattern "<<pattern->GetId()<<std::endl;

}


//*****************************************************************************
ND::THandle<ND::TTPCPattern> ND::TTPCLikelihoodMerge::MergeAll(){
//*****************************************************************************

  if ( ND::tpcDebug().LikelihoodMerge(DB_VERBOSE))
    std::cout<<" ------ Merge the paths and patterns !"<<std::endl;

  std::vector< ND::THandle<ND::TTPCPath> > MergedPaths;

  /////// New pattern !!!!!!
  ND::THandle<ND::TTPCPattern> NewPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern() );
  TTPCT0 T0(ND::tpcCalibration().GetDefaultT0());

  // 1) Create the new merged paths by creating a chain of paths that need to be merged.
  int Chain = 0;
  for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
    if( ! fMTracker[mt].NeedsMerging() )
      continue;
    // 1.a) First list the fMTracker entries that will be merged

    // Defines the two paths at the end of the chain of paths
    // which can be of any length.
    ND::THandle<ND::TTPCPath> NewPathEnds[2];
    for ( unsigned int i = 0; i < 2; i++){
      NewPathEnds[i] = fMTracker[mt].GetRawPath(i);
    }
    fMTracker[mt].SetMergingChain(Chain);
    if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
      std::cout<<" Start merging chain with paths "<<NewPathEnds[0]->GetId()<<" and "<<NewPathEnds[1]->GetId()<<std::endl;


    // We need to copy the T0 information to the new pattern.
    // Just take the T0 from the first Path.
    T0 = fMTracker[mt].GetRawPath(0)->GetTTPCT0();

    // Find the other matches involving the two paths just merged
    for ( int submt = 0; submt < int(fMTracker.size()); submt++){
      if( ! fMTracker[submt].NeedsMerging() )
        continue;
      for ( unsigned int i = 0; i < 2; i++){
        if (fMTracker[submt].HasThisPath(NewPathEnds[i])){
          // Find which of the two paths must be added to find the new end of the chain.
          for (unsigned int j = 0; j < 2; j++){
            ND::THandle<ND::TTPCPath> PathB = fMTracker[submt].GetRawPath(j);
            if (PathB == NewPathEnds[i])
              continue;
            NewPathEnds[i] = PathB;
            if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
              std::cout<<"       merging chain with paths "<<NewPathEnds[0]->GetId()<<" and "<<NewPathEnds[1]->GetId()<<" at the ends"<<std::endl;
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
    ND::THandle<ND::TTPCPath> NewPath;
    for ( unsigned int submt = 0; submt < fMTracker.size(); submt++){
      if ( ! fMTracker[submt].IsMergingChainOk(Chain))
        continue;
      // First use both paths in the fMTracker entry
      if (!NewPath){
        NewPath = TTPCUtils::MergePaths(fMTracker[submt].GetRawPath(0), fMTracker[submt].GetRawPath(1));
        for ( unsigned int i = 0; i < 2; i++){
          MergedPaths.push_back(fMTracker[submt].GetRawPath(i));
          NewPath->AddConstituent(fMTracker[submt].GetRawPath(i));
        }
        if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
          std::cout<<"   + NewPath Id "<<NewPath->GetId()<<" created with path Id "<<fMTracker[submt].GetRawPath(0)->GetId()<<" and path Id "<<fMTracker[submt].GetRawPath(1)->GetId()<<std::endl;

        continue;
      }
      
      // Here we need to find which of the two paths must be merged
      for (unsigned int j = 0; j < 2; j++){
        ND::THandle<ND::TTPCPath> PathB = fMTracker[submt].GetRawPath(j);
        std::vector< ND::THandle<ND::TTPCPath> >::iterator it = find (MergedPaths.begin(), MergedPaths.end(), PathB);
        if (it != MergedPaths.end())
          continue;
        int tmpId = NewPath->GetId();
        ND::THandle<ND::TTPCPath> prevPath = NewPath;
        NewPath = TTPCUtils::MergePaths(NewPath, PathB);
        // Now this new path is at one end of NewPath
        if ( ND::tpcDebug().LikelihoodMerge(DB_VVERBOSE))
          std::cout<<"   + NewPath Id "<<NewPath->GetId()<<" created with path Id "<<PathB->GetId()<<" and path Id "<<tmpId<<std::endl;
        MergedPaths.push_back(PathB);
        // We recreate every time the NewPath so we need to carry over the constituents by hand.
        for (ND::TReconObjectContainer::iterator ptc = prevPath->GetConstituents()->begin(); ptc != prevPath->GetConstituents()->end(); ptc++){
          NewPath->AddConstituent(*ptc);
        }
        NewPath->AddConstituent(PathB);
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

  // If there are new junctions, we need to create another pattern for a clean track
  if (!NbNewJunctions) {
    NewPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern(fMTracker[0].GetMergedPath()) );
  }
  
  NewPattern->SetId(ND::tpcCalibration().GetPatternId());
  NewPattern->InitialSetup();
  NewPattern->SetT0(T0);

  return NewPattern;
}


