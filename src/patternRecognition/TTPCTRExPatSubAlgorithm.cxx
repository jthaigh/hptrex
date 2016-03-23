// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"

ND::TTPCTRExPatSubAlgorithm::TTPCTRExPatSubAlgorithm(ND::TTPCLayout* layout) {
  // initial values for containing hits, tpc number and event and sub-event numbers
  fHasHits = false;
  fHasValidPaths = false;
  fPrimary = false;
  fTPC = 0;

  // set up objects for layout, hit group manager, feature finder and path finder
  fLayout = layout;

  fVolGroupMan = new ND::TTPCVolGroupMan(fLayout);
  fAStar = new ND::TTPCAStar(fLayout);  
  //  fPattern = new ND::TTPCPattern;
}
ND::TTPCTRExPatSubAlgorithm::~TTPCTRExPatSubAlgorithm(){
  // delete groups of hits and path finders
  delete fVolGroupMan;
  delete fAStar;

  //MDH
  //Nobody had better use the pattern after we destroy the subalgorithm...
  //This is an output object that needs to be thoroughly re-engineered anyway
  //delete fPattern;
}

void ND::TTPCTRExPatSubAlgorithm::SetUpHits(std::map<long, ND::TTPCUnitVolume*>& map, ND::TTPCAStar* aStarCopy){
  // add new pattern recognition cells if the charge cut is met
  fHitMap=map;
  
  fHasHits = fHitMap.size() > 0;
  // add map of hits to group manager object
  fVolGroupMan->AddPrimaryHits(fHitMap);
  // add map of hits to path finder
  fAStar->AddHits(fVolGroupMan, fHitMap);
}

//MDH
//Appears not to be used anywhere and is problematic since it
//creates/deletes TTPCUnitVolumes
/*
void ND::TTPCTRExPatSubAlgorithm::AbsorbCell(long id, ND::TTPCUnitVolume* cell){
  // if cell doesn't already exist in map, create it; if it does, increment its charge by cell's charge
  std::map<long, ND::TTPCUnitVolume*>::iterator el = fHitMap.find(id);
  if(el == fHitMap.end()){
    fHitMap[id] = cell;
  }
  else{
    fHitMap[id]->AddCharge(cell->GetQ());
    fHitMap[id]->AddHits(cell->GetHits());
    delete cell;
  };
}

void ND::TTPCTRExPatSubAlgorithm::AppendHits(ND::TTPCVolGroup& hits){
  // loop over all cells in hit group and absorb each one into this
  for(std::map<long, ND::TTPCUnitVolume*>::iterator it = hits->begin(); it != hits->end(); ++it){
    AbsorbCell(it->first, it->second);
  };
}
*/

void ND::TTPCTRExPatSubAlgorithm::ProduceContainers(){
  // ignore if hits don't already exist
  if(!fHasHits) return;

  // start out having no valid paths by default
  fHasValidPaths = false;

  // find list of hit groups on edge of paths
  std::vector< ND::TTPCVolGroup > edgeGroups=fVolGroupMan->GetEdgeGroups();

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  fAStar->ClearRedundancies(fVolGroupMan, edgeGroups);

  // connect pairs of edge groups to find missing hits
  std::vector< ND::TTPCOrderedVolGroup > edgePaths;
  fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, edgePaths, false, false);

  // find unmatched hits
  std::map<long, ND::TTPCUnitVolume*> unmatchedVolumes;
  unmatchedVolumes = fHitMap;
  
  for(std::vector< ND::TTPCOrderedVolGroup>::iterator edgePathIt = edgePaths.begin(); edgePathIt != edgePaths.end(); ++edgePathIt){
    ND::TTPCOrderedVolGroup& edgePath = *edgePathIt;

    fVolGroupMan->BuildGroupFriends(edgePath, ND::TTPCConnection::extraHits);
    ND::TTPCVolGroup pathHits=edgePath.GetExtendedHits();

    for(std::map<long, ND::TTPCUnitVolume*>::iterator foundVolIt = pathHits.begin(); foundVolIt != pathHits.end(); ++foundVolIt){
      long id = foundVolIt->first;
      std::map<long, ND::TTPCUnitVolume*>::iterator foundIdLoc = unmatchedVolumes.find(id);
      if(foundIdLoc != unmatchedVolumes.end()) unmatchedVolumes.erase(foundIdLoc);
    };
  };

  // split and group unmatched hits
  ND::TTPCVolGroupMan secondPassVolume(fLayout);
  secondPassVolume.AddPrimaryHits(unmatchedVolumes);
  std::vector< ND::TTPCVolGroup > extraSubVolumes;
  secondPassVolume.GetConnectedHits(extraSubVolumes,ND::TTPCConnection::path, ND::TTPCHitGroupings::all, true);

  for(std::vector< ND::TTPCVolGroup >::iterator extraSubVolumeIt = extraSubVolumes.begin(); extraSubVolumeIt != extraSubVolumes.end(); ++extraSubVolumeIt){
    ND::TTPCVolGroup& extraSubVolume = *extraSubVolumeIt;

    ND::TTPCVolGroupMan subVolumeMan(fLayout);

    //MDH
    //Hmm, bit clumsy. Add a getter for the TTPCVolGroup hit map?
    std::map<long, ND::TTPCUnitVolume*> subVolumeHits;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator volIt = extraSubVolume.begin(); volIt != extraSubVolume.end(); ++volIt){
      subVolumeHits[volIt->first] = volIt->second;
    }

    subVolumeMan.AddPrimaryHits(subVolumeHits);
    std::vector<ND::TTPCVolGroup> extraEdgeGroups=subVolumeMan.GetEdgeGroups();
    // add unmatched edges to main group
    for(std::vector< ND::TTPCVolGroup >::iterator extraEdgeGroupIt = extraEdgeGroups.begin(); extraEdgeGroupIt != extraEdgeGroups.end(); ++extraEdgeGroupIt){
      edgeGroups.emplace_back(std::move(*extraEdgeGroupIt));
    }
  }

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  fAStar->ClearRedundancies(fVolGroupMan, edgeGroups);


  // connect pairs of edge groups

  // define container for true paths before corner detection
  std::vector< ND::TTPCOrderedVolGroup > truePaths;

  // define vertices for later use
  std::vector< ND::TTPCVolGroup > vertices;

  // if there are at least three edge groups, look for vertices to connect to them
  if(edgeGroups.size() > 2){

    edgePaths.clear();

    // connect pairs of edge groups
    fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, edgePaths, false, true);
    // look for vertices from edge group pairs
    fVolGroupMan->GetFoci(edgePaths, vertices ,.1);
    // make sure number of vertices is two less than number of track ends
    fVolGroupMan->CleanUpVertices(edgeGroups, vertices);
    // make sure vertices contain all the right hits
    fVolGroupMan->BulkGroups(vertices);

    if(vertices.size() > 0){
      // connect vertices to path ends to get paths
      fAStar->ConnectVertexGroupsOrdered(fVolGroupMan, vertices, edgeGroups,truePaths);
      // clear redundant paths
      fAStar->ClearVertexConnectionRedundancies(fVolGroupMan, truePaths, vertices);
    }
    // otherwise, return empty handed
    else{
      return;
    };
  }
  // if there are two, connect the two edge groups
  else if (edgeGroups.size() > 1){
    fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, truePaths);
  }
  // otherwise, return empty handed
  else{
    return;
  };

  // ensure that paths don't share hits with each other
  fVolGroupMan->BuildAllFriends(truePaths);

  // extend and cluster paths
  for(std::vector< ND::TTPCOrderedVolGroup >::iterator truePathIt = truePaths.begin(); truePathIt != truePaths.end(); ++truePathIt){
    fVolGroupMan->ClusterGroupFriends(*truePathIt, false, true);
  };
  // clear empties
  std::vector< ND::TTPCOrderedVolGroup >::iterator pathDel = truePaths.begin();
  while(pathDel != truePaths.end()){
    if(!pathDel->size()) truePaths.erase(pathDel);
    else pathDel ++;
  }

  // look for kinks in paths
  fVolGroupMan->BreakPathsAboutKinks(truePaths, true);
  

  for(std::vector< ND::TTPCOrderedVolGroup >::iterator brokenPathIt = truePaths.begin(); brokenPathIt != truePaths.end(); ++brokenPathIt){
    fVolGroupMan->ClusterGroupFriends(*brokenPathIt, true, true);
  };

  // temporarily merge x-paths into junctions to make the next two steps easier
  fVolGroupMan->SeparateXPathHits(truePaths);

  // clean up any and all instances of hits beiing placed in the wrong object
  fVolGroupMan->SeparateEmptyClusters(truePaths);
  fVolGroupMan->SeparateClusterHits(truePaths);
  fVolGroupMan->SeparateAnomHits(truePaths);
  fVolGroupMan->SeparateJunctionHits(truePaths);

  // enforce path ordering
  fVolGroupMan->EnforceOrdering(truePaths);

  // add all paths that make sense and store if none are found
  fVolGroupMan->SanityFilter(truePaths);

  // associate any remaining unused hits to junctions where possible
  ND::TTPCVolGroup unusedHits(fLayout);
  fVolGroupMan->GetUnusedHits(truePaths,unusedHits);
  fVolGroupMan->AssociateUnusedWithJunctions(unusedHits, truePaths);

  // reset vertex status based on how many paths a junction is connected to
  fVolGroupMan->ResetVertexStatuses(truePaths);

  fTracks=std::move(truePaths);
  fHasValidPaths = fTracks.size() > 0;
}

//ND::TTPCPattern* ND::TTPCTRExPatSubAlgorithm::GetPattern(){
//  if(!fPattern){
//    return 0;
//  };
//  return fPattern;
//}

//MDH
//Get rid of this for now since it needs to be rewritten for new output objects
/*
void ND::TTPCTRExPatSubAlgorithm::ProducePattern(ND::THitSelection* used){
  if(fPattern){
    return;
  };

  // final paths and junctions
  std::vector< ND::THandle<ND::TTPCPath> > paths = std::vector< ND::THandle<ND::TTPCPath> >();
  std::vector< ND::THandle<ND::TTPCJunction> > junctions = std::vector< ND::THandle<ND::TTPCJunction> >();

  // temporary vectors for holding candidate paths and junctions, and parallel vector of groups associated with junctions
  std::vector< ND::THandle<ND::TTPCPath> > tempPaths = std::vector< ND::THandle<ND::TTPCPath> >();
  std::vector< ND::THandle<ND::TTPCJunction> > tempJunctions = std::vector< ND::THandle<ND::TTPCJunction> >();
  std::vector< ND::THandle<ND::TTPCVolGroup> > tempJunctionGroups = std::vector< ND::THandle<ND::TTPCVolGroup> >();

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    ND::THandle<ND::TTPCOrderedVolGroup> track = *trackIt;

    // ignore if empty
    ND::THandle<ND::THitSelection> clusters = track->GetClusters();
    if(!clusters->size()) continue;

    // create new path
    ND::THandle<ND::TTPCPath> path( new ND::TTPCPath() );
    tempPaths.push_back(path);

    // otherwise add all hits to path
    ND::THitSelection* copyClusters = new ND::THitSelection(*clusters);
    path->AddHits(copyClusters);

    // also save if delta criteria met
    if(track->GetDeltaCriteriaMet()){
      path->SetPID(ND::TReconPID::kEM, 1.);
    };
    // and save if the path is in the x direction
    if(track->GetIsXPath()){
      path->SetIsXPath(true);
    };

    // save groups for junctions
    std::vector< ND::THandle<ND::TTPCVolGroup> > junctionGroups = std::vector< ND::THandle<ND::TTPCVolGroup> >();
    if(track->HasBackHits()) junctionGroups.push_back(track->GetBackHits());
    if(track->HasFrontHits()) junctionGroups.push_back(track->GetFrontHits());
    // determine if junction candidate has already been found and create a group for it if it hasn't
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionGroupIt = junctionGroups.begin(); junctionGroupIt != junctionGroups.end(); ++junctionGroupIt){
      ND::THandle<ND::TTPCVolGroup> junctionGroup = *junctionGroupIt;
      bool found = false;

      int iMax = tempJunctions.size();
      for(int i=0; i<iMax; i++){
        ND::THandle<ND::TTPCJunction> tempJunction = tempJunctions[i];
        ND::THandle<ND::TTPCVolGroup> tempJunctionGroup = tempJunctionGroups[i];

        // if the junction has the same id as the temporary, add this path to it
        if(junctionGroup->GetID() == tempJunctionGroup->GetID()){
          // if already found, check this path isn't already added and add this path
          bool pathInJunction = false;
          for(ND::TReconObjectContainer::iterator junctionConstIt = tempJunction->GetConstituents()->begin(); junctionConstIt != tempJunction->GetConstituents()->end(); ++junctionConstIt){
            if(path == *junctionConstIt){
              pathInJunction = true;
              break;
            };
          };
          if(!pathInJunction) tempJunction->AddConstituent(path);

          found = true;
        };
      };
      // if not alread found, create a new one
      if(!found){
        ND::THitSelection* junctionHitSelection = junctionGroup->GetHitSelection();
        TVector3 tempVect3 = junctionGroup->GetAverageVol()->GetPos();
        double tempTime = junctionGroup->GetAverageVol()->GetTime();
        // Get a guess X position based on a default T0 of 0.0 just to have a useful X position downstream.
        double tempX = TTPCUtils::GetXWithT0(tempTime, 0.0, (*(junctionHitSelection->begin())));
        tempVect3.SetX(tempX);
        TLorentzVector tempPosition(tempVect3, tempTime);

        ND::THandle<ND::TTPCJunction>  junct( new ND::TTPCJunction(tempPosition) );
        junct->AddHits(junctionHitSelection);
        junct->AddConstituent(path);

        // make sure these two vectors remain parallel
        tempJunctions.push_back(junct);
        tempJunctionGroups.push_back(junctionGroup);
      };
    };
  };

  // add paths that meet criteria
  for(std::vector< ND::THandle<ND::TTPCPath> >::iterator path = tempPaths.begin(); path != tempPaths.end(); ++path){
    (*path)->SetId(ND::tpcCalibration().GetPathId());
    paths.push_back(*path);
  };

  // add junctions that meet criteria
  for(std::vector< ND::THandle<ND::TTPCJunction> >::iterator junction = tempJunctions.begin(); junction != tempJunctions.end(); ++junction){
    if((*junction)->GetNPaths() > 1){
      (*junction)->SetId(ND::tpcCalibration().GetJunctionId());
      junctions.push_back(*junction);
    }
  };

  // add junctions only.  paths will be extracted from the junctions
  if(junctions.size() == 0){
    if(paths.size() == 1){
      ND::THandle<ND::TTPCPath> thepath = *(paths.begin());
      fPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern(thepath) );
      // One clean single track. Special treatment.
    }
    else{
      // something went wrong - probably an unreconstructable event.  save a null pointer so the main algorithm knows not to add it to list of patterns.
      if(ND::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << "WARNING: no pattern reconstructed " << std::endl;
      fPattern = ND::THandle<ND::TTPCPattern>();
    };
  } else {
    fPattern = ND::THandle<ND::TTPCPattern>( new ND::TTPCPattern() );
    for(std::vector< ND::THandle<ND::TTPCJunction> >::iterator junction = junctions.begin(); junction != junctions.end(); ++junction)
      fPattern->AddJunction(*junction);
  };
  // fill used and unused hits
  FillUsedHits(used);

  // Call the setup method to define things like the detector bit.
  if(fPattern){
    fPattern->SetId(ND::tpcCalibration().GetPatternId());
    fPattern->InitialSetup();
  }
  }*/

std::vector<ND::TTPCHitPad*> ND::TTPCTRExPatSubAlgorithm::GetHits(){
  return fVolGroupMan->GetHits();
}
std::vector<ND::TTPCHitPad*> ND::TTPCTRExPatSubAlgorithm::GetHits(ND::TTPCOrderedVolGroup& path){
  return std::move(path.GetClusters());
}

void ND::TTPCTRExPatSubAlgorithm::GetRegions(std::vector< ND::TTPCVolGroup >& regions){

  if(!fHasHits){
    regions.clear();
    return;
  }
  // get list of hit groups for sub-events
  fVolGroupMan->GetConnectedHits(regions,ND::TTPCConnection::path);
}

//MDH
//Not currently used. May be used later if we resurrect ProducePattern.
/*
void ND::TTPCTRExPatSubAlgorithm::FillUsedHits(std::vector<ND::TTPCHitPad*>& used){
  // if pattern hasn't been correctly built, hits can't possibly be used
  if(!fPattern) return;

  // matched (initially empty) volumes
  std::map<long, ND::TTPCUnitVolume*> matchedVolumes = std::map<long, ND::TTPCUnitVolume*>();

  // find those contained in any paths or junctions
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    ND::THandle<ND::TTPCOrderedVolGroup> track = *trackIt;
    // path hits
    for(std::vector<ND::TTPCPathVolume*>::iterator pathIt = track->begin(); pathIt != track->end(); ++pathIt){
      ND::TTPCPathVolume* path = *pathIt;

      std::vector<ND::TTPCUnitVolume*> pathClusterVols = path->GetExtendedCell(0);
      for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathClusterVols.begin(); volIt != pathClusterVols.end(); ++volIt){
        std::map<long, ND::TTPCUnitVolume*>::iterator foundVol = fHitMap.find((*volIt)->GetID());
        matchedVolumes[foundVol->first] = foundVol->second;
      };
    };
  };

  // junction hits
  std::vector< ND::THandle<ND::TTPCVolGroup> > pathVertices = fVolGroupMan->GetJunctionsFromPaths(fTracks);
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator pathVertexIt = pathVertices.begin(); pathVertexIt != pathVertices.end(); ++pathVertexIt){
    ND::THandle<ND::TTPCVolGroup> pathVertex = *pathVertexIt;

    for(std::map<long, ND::TTPCUnitVolume*>::iterator volIt = pathVertex->begin(); volIt != pathVertex->end(); ++volIt){
      std::map<long, ND::TTPCUnitVolume*>::iterator foundVol = fHitMap.find(volIt->first);
      matchedVolumes[foundVol->first] = foundVol->second;
    };
  };

  // add used hits
  for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = matchedVolumes.begin(); volEl != matchedVolumes.end(); ++volEl){
    ND::TTPCUnitVolume* vol = volEl->second;
    if(!vol) continue;

    std::vector< ND::THandle<ND::TTPCHitPad> > hits = vol->GetHits();
    for(std::vector< ND::THandle<ND::TTPCHitPad> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){
      ND::THandle<ND::TTPCHitPad> hitPad = *hitIt;
      if(hitPad) used->AddHit(hitPad);
    };
  };
}
*/
void ND::TTPCTRExPatSubAlgorithm::GetTrackExtendedHits(std::vector<ND::TTPCVolGroup>& extHits){

  for(std::vector< ND::TTPCOrderedVolGroup >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    ND::TTPCOrderedVolGroup& track = *trackIt;
    extHits.push_back(track.GetExtendedHits());
  }
  
}
