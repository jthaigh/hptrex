#include "TTPCLinearMerge.hxx"
#include "TTPCUtils.hxx"
#include <map>
#include "TGraph.h"
#include "TFitResult.h"

#include "TFile.h"

//*****************************************************************************
trex::TTPCLinearMerge::TTPCLinearMerge(){
//*****************************************************************************

  fRunThroughGoingMerging = true;
  //fRunBrokenTracksMerging = true;
  fThruGoDeltaLogLklhdCut = true;
  //fBrkTrkDeltaLogLklhdCut = true;

  fThruGoDeltaLogLklhdCut = 200;
  //fBrkTrkDeltaLogLklhdCut = 500;
}


//*****************************************************************************
void trex::TTPCLinearMerge::Process(std::vector<trex::TTRExPattern>& inputPatterns, std::vector<trex::TTRExPattern>& mergedPatterns){
//*****************************************************************************

  for (auto patit = inputPatterns.begin(); patit != inputPatterns.end(); patit++) {
    trex::TTRExPattern& Pattern = *patit;
    Pattern.SetUsable(true);
  }

  // 1) Loop over the patterns
  for (auto basePatit = inputPatterns.begin(); basePatit != inputPatterns.end(); basePatit++) {
    trex::TTRExPattern& BasePattern = *basePatit;
    
    MatchThroughJunctions(BasePattern);
    
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
    MergeAll(BasePattern,mergedPatterns);
    
    // Clean up after each chain of matching.
    CleanUp();
  }
}


//*****************************************************************************
void trex::TTPCLinearMerge::CleanUp(){
//*****************************************************************************
  fMTracker.clear();
 }

double trex::TTPCLinearMerge::CalcChi2PerDOF( trex::TTRExPath& PathA, trex::TTRExPath& PathB){
  TGraph grMatch(1);
  unsigned int iPt=0;

  for(auto clIter=PathA.GetClusters().begin();clIter!=PathA.GetClusters().end();++clIter){
    for(auto hitIter=(*clIter)->GetClusterHits().begin();hitIter!=(*clIter)->GetClusterHits().end();++hitIter){
      TVector3 hitPos=(*hitIter)->GetPosition();
      grMatch.SetPoint(iPt++,hitPos.Z(),hitPos.Y());
    }
  }

  for(auto clIter=PathB.GetClusters().begin();clIter!=PathB.GetClusters().end();++clIter){
    for(auto hitIter=(*clIter)->GetClusterHits().begin();hitIter!=(*clIter)->GetClusterHits().end();++hitIter){
      TVector3 hitPos=(*hitIter)->GetPosition();
      grMatch.SetPoint(iPt++,hitPos.Z(),hitPos.Y());
    }
  }

  TFitResultPtr rpMatch = grMatch.Fit("pol3","S");

  return rpMatch->Chi2()/rpMatch->Ndf();

}


//*****************************************************************************
void trex::TTPCLinearMerge::MatchThroughJunctions(trex::TTRExPattern& pattern){
//*****************************************************************************
  if ( !fRunThroughGoingMerging)
    return;

  if (pattern.GetPaths().size() == 1)
    return;

  ////// Look for a match across each junction
  for (auto constit = pattern.GetJunctions().begin(); constit != pattern.GetJunctions().end(); constit++) {
    bool atLeastOneMatch = false;
    trex::TTRExJunction& junction = (*constit);
    TVector3 posJunct=junction.GetPosition();

    std::vector< trex::TTRExPath* > Paths;
    for (auto pathIter = junction.GetConnectedPaths().begin(); pathIter != junction.GetConnectedPaths().end(); pathIter++) {
      Paths.push_back(*pathIter);
    }
    
    while(Paths.size()>1){
      auto bestPathA=junction.GetConnectedPaths().begin();
      auto bestPathB=junction.GetConnectedPaths().begin();
      double bestChi2=1.E13;
      bool hasMatch=false;

      for (auto pathA = Paths.begin(); pathA != Paths.end(); pathA++) {
	
	TVector3 posAStart=(*pathA)->GetClusters().front()->GetPosition();
	TVector3 posAEnd=(*pathA)->GetClusters().back()->GetPosition();

	if((posAStart-posJunct).Mag()>(posAEnd-posJunct).Mag()){
	  posAEnd=posAStart;
	}

	for (auto pathB = pathA+1; pathB != Paths.end(); pathB++) {
	  TVector3 posBStart=(*pathB)->GetClusters().front()->GetPosition();
	  TVector3 posBEnd=(*pathB)->GetClusters().back()->GetPosition();
	  
	  if((posBStart-posJunct).Mag()>(posBEnd-posJunct).Mag()){
	    posBEnd=posBStart;
	  }
	  
	  //If both track ends are above or below the junction then
	  //sense is wrong - skip
	  if((posAEnd-posJunct).Z()*(posBEnd-posJunct).Z()>0){
	    continue;
	  }
	  
	  double thisChi2=CalcChi2PerDOF(**pathA,**pathB);
	  if(thisChi2<bestChi2){
	    bestChi2=thisChi2;
	    bestPathA=pathA;
	    bestPathB=pathB;
	    hasMatch=true;
	  }
	}
      }
      
      if((!hasMatch)||bestChi2>10.){
	break;
      }
      fMTracker.push_back( MatchingTracker(**bestPathA, **bestPathB, junction) );
      Paths.erase(bestPathB);
      Paths.erase(bestPathA);
      atLeastOneMatch=true;
    }
    if(!atLeastOneMatch){
      fMTracker.push_back( MatchingTracker(junction) );
    }
  }
}     

//*****************************************************************************
void trex::TTPCLinearMerge::MergeAll(trex::TTRExPattern& inputPattern, std::vector<trex::TTRExPattern>& outputVector){
//*****************************************************************************

  /////// New pattern !!!!!!
  outputVector.emplace_back();
  trex::TTRExPattern* NewPattern = &(outputVector.back());

  std::vector< trex::TTRExPath* > MergedPaths;

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
  for (auto constit = inputPattern.GetPaths().begin(); constit != inputPattern.GetPaths().end(); constit++) {
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
    }
  }

  // 3) Recreate all the junctions and copy the unmerged paths and the new merged paths correctly.
  int NbNewJunctions = 0;
  std::map<TTRExPath*,TTRExPath*> tmpToPatternPathsMap;
  for ( unsigned int mt = 0; mt < fMTracker.size(); mt++){
    // Skip matching between two segments of a broken which doesn't have a junction.
    trex::TTRExJunction* tmpJunc= fMTracker[mt].GetJunction();
    if( !tmpJunc){
      continue;
    }

    NbNewJunctions++;
    NewPattern->GetJunctions().emplace_back();
    NewPattern->GetPaths().reserve(1000);
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
	    if(!tmpToPatternPathsMap.count(fMTracker[submt].GetMergedPath())){
	      NewPattern->GetPaths().emplace_back(*(fMTracker[submt].GetMergedPath()));
	      NewJunction->AddConnectedPath(&(NewPattern->GetPaths().back()));
	      tmpToPatternPathsMap[fMTracker[submt].GetMergedPath()]=&(NewPattern->GetPaths().back());
	    }
	    else{
	      NewJunction->AddConnectedPath(tmpToPatternPathsMap[fMTracker[submt].GetMergedPath()]);
	    }	    
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


/*
//*****************************************************************************
double trex::TTPCLikelihoodMerge::PathToPatternMatch( trex::TTRExPath& PathA, trex::TTRExPattern& PatternB, trex::TTRExPath* bestPathB){
  double bestDLL = 1.e13;
  unsigned int bestEndA = 2;
  unsigned int bestEndB = 2;

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
*/
