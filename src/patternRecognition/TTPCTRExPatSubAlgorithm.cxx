// eddy
#include "TTPCTRExPatSubAlgorithm.hxx"

trex::TTPCTRExPatSubAlgorithm::TTPCTRExPatSubAlgorithm(trex::TTPCLayout* layout) {
  // initial values for containing hits, tpc number and event and sub-event numbers
  fHasHits = false;
  fHasValidPaths = false;
  fPrimary = false;
  fTPC = 0;

  // set up objects for layout, hit group manager, feature finder and path finder
  fLayout = layout;

  fVolGroupMan = new trex::TTPCVolGroupMan(fLayout);
  fAStar = new trex::TTPCAStar(fLayout);  
  //  fPattern = new trex::TTPCPattern;
}
trex::TTPCTRExPatSubAlgorithm::~TTPCTRExPatSubAlgorithm(){
  // delete groups of hits and path finders
  if(fVolGroupMan) delete fVolGroupMan;
  if(fAStar) delete fAStar;

  //MDH
  //Nobody had better use the pattern after we destroy the subalgorithm...
  //This is an output object that needs to be thoroughly re-engineered anyway
  //delete fPattern;
}

trex::TTPCTRExPatSubAlgorithm::TTPCTRExPatSubAlgorithm(trex::TTPCTRExPatSubAlgorithm&& in){
  fHitMap=std::move(in.fHitMap);
  fVolGroupMan=in.fVolGroupMan;
  fAStar=in.fAStar;
  fTracks=std::move(in.fTracks);
  fHasHits=in.fHasHits;
  fHasValidPaths=in.fHasValidPaths;
  fTPC=in.fTPC;
  fLayout=in.fLayout;
  fPrimary=in.fPrimary;

  in.fAStar=0;
  in.fVolGroupMan=0;
}

void trex::TTPCTRExPatSubAlgorithm::SetUpHits(std::map<long, trex::TTPCUnitVolume*>& map, trex::TTPCAStar* aStarCopy){

  fHitMap=map;
  
  fHasHits = fHitMap.size() > 0;
  // add map of hits to group manager object
  fVolGroupMan->AddPrimaryHits(fHitMap);
  // add map of hits to path finder
  fAStar->AddHits(fVolGroupMan, fHitMap);
}

void trex::TTPCTRExPatSubAlgorithm::ProduceContainers(){
  // ignore if hits don't already exist
  if(!fHasHits) return;

  // start out having no valid paths by default
  fHasValidPaths = false;

  // find list of hit groups on edge of paths
  std::vector< trex::TTPCVolGroup > edgeGroups=fVolGroupMan->GetEdgeGroups();

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  fAStar->ClearRedundancies(fVolGroupMan, edgeGroups);

  // connect pairs of edge groups to find missing hits
  std::vector< trex::TTPCOrderedVolGroup > edgePaths;
  fAStar->ConnectGroupsOrdered(fVolGroupMan, edgeGroups, edgePaths, false, false);

  // find unmatched hits
  std::map<long, trex::TTPCUnitVolume*> unmatchedVolumes;
  unmatchedVolumes = fHitMap;
  
  for(std::vector< trex::TTPCOrderedVolGroup>::iterator edgePathIt = edgePaths.begin(); edgePathIt != edgePaths.end(); ++edgePathIt){
    trex::TTPCOrderedVolGroup& edgePath = *edgePathIt;

    fVolGroupMan->BuildGroupFriends(edgePath, trex::TTPCConnection::extraHits);
    trex::TTPCVolGroup pathHits=edgePath.GetExtendedHits();

    for(std::map<long, trex::TTPCUnitVolume*>::iterator foundVolIt = pathHits.begin(); foundVolIt != pathHits.end(); ++foundVolIt){
      long id = foundVolIt->first;
      std::map<long, trex::TTPCUnitVolume*>::iterator foundIdLoc = unmatchedVolumes.find(id);
      if(foundIdLoc != unmatchedVolumes.end()) unmatchedVolumes.erase(foundIdLoc);
    };
  };

  // split and group unmatched hits
  trex::TTPCVolGroupMan secondPassVolume(fLayout);
  secondPassVolume.AddPrimaryHits(unmatchedVolumes);
  std::vector< trex::TTPCVolGroup > extraSubVolumes;
  secondPassVolume.GetConnectedHits(extraSubVolumes,trex::TTPCConnection::path, trex::TTPCHitGroupings::all, true);

  for(std::vector< trex::TTPCVolGroup >::iterator extraSubVolumeIt = extraSubVolumes.begin(); extraSubVolumeIt != extraSubVolumes.end(); ++extraSubVolumeIt){
    trex::TTPCVolGroup& extraSubVolume = *extraSubVolumeIt;

    trex::TTPCVolGroupMan subVolumeMan(fLayout);

    //MDH
    //Hmm, bit clumsy. Add a getter for the TTPCVolGroup hit map?
    std::map<long, trex::TTPCUnitVolume*> subVolumeHits;
    for(std::map<long, trex::TTPCUnitVolume*>::iterator volIt = extraSubVolume.begin(); volIt != extraSubVolume.end(); ++volIt){
      subVolumeHits[volIt->first] = volIt->second;
    }

    subVolumeMan.AddPrimaryHits(subVolumeHits);
    std::vector<trex::TTPCVolGroup> extraEdgeGroups=subVolumeMan.GetEdgeGroups();
    // add unmatched edges to main group
    for(std::vector< trex::TTPCVolGroup >::iterator extraEdgeGroupIt = extraEdgeGroups.begin(); extraEdgeGroupIt != extraEdgeGroups.end(); ++extraEdgeGroupIt){
      edgeGroups.emplace_back(std::move(*extraEdgeGroupIt));
    }
  }

  // clear redundant edge groups (i.e. ones which are actually just half way points along existing paths)
  fAStar->ClearRedundancies(fVolGroupMan, edgeGroups);


  // connect pairs of edge groups

  // define container for true paths before corner detection
  std::vector< trex::TTPCOrderedVolGroup > truePaths;

  // define vertices for later use
  std::vector< trex::TTPCVolGroup > vertices;

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
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator truePathIt = truePaths.begin(); truePathIt != truePaths.end(); ++truePathIt){
    fVolGroupMan->ClusterGroupFriends(*truePathIt, false, true);
  };
  // clear empties
  std::vector< trex::TTPCOrderedVolGroup >::iterator pathDel = truePaths.begin();
  while(pathDel != truePaths.end()){
    if(!pathDel->size()) truePaths.erase(pathDel);
    else pathDel ++;
  }

  // look for kinks in paths
  fVolGroupMan->BreakPathsAboutKinks(truePaths);
  

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator brokenPathIt = truePaths.begin(); brokenPathIt != truePaths.end(); ++brokenPathIt){
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
  trex::TTPCVolGroup unusedHits(fLayout);
  fVolGroupMan->GetUnusedHits(truePaths,unusedHits);
  fVolGroupMan->AssociateUnusedWithJunctions(unusedHits, truePaths);

  // reset vertex status based on how many paths a junction is connected to
  fVolGroupMan->ResetVertexStatuses(truePaths);

  fTracks=std::move(truePaths);
  fHasValidPaths = fTracks.size() > 0;
}

//trex::TTPCPattern* trex::TTPCTRExPatSubAlgorithm::GetPattern(){
//  if(!fPattern){
//    return 0;
//  };
//  return fPattern;
//}

//MDH
//Get rid of this for now since it needs to be rewritten for new output objects
/*
void trex::TTPCTRExPatSubAlgorithm::ProducePattern(trex::THitSelection* used){
  if(fPattern){
    return;
  };

  // final paths and junctions
  std::vector< trex::THandle<trex::TTPCPath> > paths = std::vector< trex::THandle<trex::TTPCPath> >();
  std::vector< trex::THandle<trex::TTPCJunction> > junctions = std::vector< trex::THandle<trex::TTPCJunction> >();

  // temporary vectors for holding candidate paths and junctions, and parallel vector of groups associated with junctions
  std::vector< trex::THandle<trex::TTPCPath> > tempPaths = std::vector< trex::THandle<trex::TTPCPath> >();
  std::vector< trex::THandle<trex::TTPCJunction> > tempJunctions = std::vector< trex::THandle<trex::TTPCJunction> >();
  std::vector< trex::THandle<trex::TTPCVolGroup> > tempJunctionGroups = std::vector< trex::THandle<trex::TTPCVolGroup> >();

  for(std::vector< trex::THandle<trex::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    trex::THandle<trex::TTPCOrderedVolGroup> track = *trackIt;

    // ignore if empty
    trex::THandle<trex::THitSelection> clusters = track->GetClusters();
    if(!clusters->size()) continue;

    // create new path
    trex::THandle<trex::TTPCPath> path( new trex::TTPCPath() );
    tempPaths.push_back(path);

    // otherwise add all hits to path
    trex::THitSelection* copyClusters = new trex::THitSelection(*clusters);
    path->AddHits(copyClusters);

    // also save if delta criteria met
    if(track->GetDeltaCriteriaMet()){
      path->SetPID(trex::TReconPID::kEM, 1.);
    };
    // and save if the path is in the x direction
    if(track->GetIsXPath()){
      path->SetIsXPath(true);
    };

    // save groups for junctions
    std::vector< trex::THandle<trex::TTPCVolGroup> > junctionGroups = std::vector< trex::THandle<trex::TTPCVolGroup> >();
    if(track->HasBackHits()) junctionGroups.push_back(track->GetBackHits());
    if(track->HasFrontHits()) junctionGroups.push_back(track->GetFrontHits());
    // determine if junction candidate has already been found and create a group for it if it hasn't
    for(std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator junctionGroupIt = junctionGroups.begin(); junctionGroupIt != junctionGroups.end(); ++junctionGroupIt){
      trex::THandle<trex::TTPCVolGroup> junctionGroup = *junctionGroupIt;
      bool found = false;

      int iMax = tempJunctions.size();
      for(int i=0; i<iMax; i++){
        trex::THandle<trex::TTPCJunction> tempJunction = tempJunctions[i];
        trex::THandle<trex::TTPCVolGroup> tempJunctionGroup = tempJunctionGroups[i];

        // if the junction has the same id as the temporary, add this path to it
        if(junctionGroup->GetID() == tempJunctionGroup->GetID()){
          // if already found, check this path isn't already added and add this path
          bool pathInJunction = false;
          for(trex::TReconObjectContainer::iterator junctionConstIt = tempJunction->GetConstituents()->begin(); junctionConstIt != tempJunction->GetConstituents()->end(); ++junctionConstIt){
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
        trex::THitSelection* junctionHitSelection = junctionGroup->GetHitSelection();
        TVector3 tempVect3 = junctionGroup->GetAverageVol()->GetPos();
        double tempTime = junctionGroup->GetAverageVol()->GetTime();
        // Get a guess X position based on a default T0 of 0.0 just to have a useful X position downstream.
        double tempX = TTPCUtils::GetXWithT0(tempTime, 0.0, (*(junctionHitSelection->begin())));
        tempVect3.SetX(tempX);
        TLorentzVector tempPosition(tempVect3, tempTime);

        trex::THandle<trex::TTPCJunction>  junct( new trex::TTPCJunction(tempPosition) );
        junct->AddHits(junctionHitSelection);
        junct->AddConstituent(path);

        // make sure these two vectors remain parallel
        tempJunctions.push_back(junct);
        tempJunctionGroups.push_back(junctionGroup);
      };
    };
  };

  // add paths that meet criteria
  for(std::vector< trex::THandle<trex::TTPCPath> >::iterator path = tempPaths.begin(); path != tempPaths.end(); ++path){
    (*path)->SetId(trex::tpcCalibration().GetPathId());
    paths.push_back(*path);
  };

  // add junctions that meet criteria
  for(std::vector< trex::THandle<trex::TTPCJunction> >::iterator junction = tempJunctions.begin(); junction != tempJunctions.end(); ++junction){
    if((*junction)->GetNPaths() > 1){
      (*junction)->SetId(trex::tpcCalibration().GetJunctionId());
      junctions.push_back(*junction);
    }
  };

  // add junctions only.  paths will be extracted from the junctions
  if(junctions.size() == 0){
    if(paths.size() == 1){
      trex::THandle<trex::TTPCPath> thepath = *(paths.begin());
      fPattern = trex::THandle<trex::TTPCPattern>( new trex::TTPCPattern(thepath) );
      // One clean single track. Special treatment.
    }
    else{
      // something went wrong - probably an unreconstructable event.  save a null pointer so the main algorithm knows not to add it to list of patterns.
      if(trex::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << "WARNING: no pattern reconstructed " << std::endl;
      fPattern = trex::THandle<trex::TTPCPattern>();
    };
  } else {
    fPattern = trex::THandle<trex::TTPCPattern>( new trex::TTPCPattern() );
    for(std::vector< trex::THandle<trex::TTPCJunction> >::iterator junction = junctions.begin(); junction != junctions.end(); ++junction)
      fPattern->AddJunction(*junction);
  };
  // fill used and unused hits
  FillUsedHits(used);

  // Call the setup method to define things like the detector bit.
  if(fPattern){
    fPattern->SetId(trex::tpcCalibration().GetPatternId());
    fPattern->InitialSetup();
  }
  }*/

std::vector<trex::TTPCHitPad*> trex::TTPCTRExPatSubAlgorithm::GetHits(){
  return fVolGroupMan->GetHits();
}
std::vector<trex::TTPCHitPad*> trex::TTPCTRExPatSubAlgorithm::GetHits(trex::TTPCOrderedVolGroup& path){
  return std::move(path.GetClusters());
}

void trex::TTPCTRExPatSubAlgorithm::GetRegions(std::vector< trex::TTPCVolGroup >& regions){

  if(!fHasHits){
    regions.clear();
    return;
  }
  // get list of hit groups for sub-events
  fVolGroupMan->GetConnectedHits(regions,trex::TTPCConnection::path);
}

//MDH
//Not currently used. May be used later if we resurrect ProducePattern.
/*
void trex::TTPCTRExPatSubAlgorithm::FillUsedHits(std::vector<trex::TTPCHitPad*>& used){
  // if pattern hasn't been correctly built, hits can't possibly be used
  if(!fPattern) return;

  // matched (initially empty) volumes
  std::map<long, trex::TTPCUnitVolume*> matchedVolumes = std::map<long, trex::TTPCUnitVolume*>();

  // find those contained in any paths or junctions
  for(std::vector< trex::THandle<trex::TTPCOrderedVolGroup> >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    trex::THandle<trex::TTPCOrderedVolGroup> track = *trackIt;
    // path hits
    for(std::vector<trex::TTPCPathVolume*>::iterator pathIt = track->begin(); pathIt != track->end(); ++pathIt){
      trex::TTPCPathVolume* path = *pathIt;

      std::vector<trex::TTPCUnitVolume*> pathClusterVols = path->GetExtendedCell(0);
      for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathClusterVols.begin(); volIt != pathClusterVols.end(); ++volIt){
        std::map<long, trex::TTPCUnitVolume*>::iterator foundVol = fHitMap.find((*volIt)->GetID());
        matchedVolumes[foundVol->first] = foundVol->second;
      };
    };
  };

  // junction hits
  std::vector< trex::THandle<trex::TTPCVolGroup> > pathVertices = fVolGroupMan->GetJunctionsFromPaths(fTracks);
  for(std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator pathVertexIt = pathVertices.begin(); pathVertexIt != pathVertices.end(); ++pathVertexIt){
    trex::THandle<trex::TTPCVolGroup> pathVertex = *pathVertexIt;

    for(std::map<long, trex::TTPCUnitVolume*>::iterator volIt = pathVertex->begin(); volIt != pathVertex->end(); ++volIt){
      std::map<long, trex::TTPCUnitVolume*>::iterator foundVol = fHitMap.find(volIt->first);
      matchedVolumes[foundVol->first] = foundVol->second;
    };
  };

  // add used hits
  for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = matchedVolumes.begin(); volEl != matchedVolumes.end(); ++volEl){
    trex::TTPCUnitVolume* vol = volEl->second;
    if(!vol) continue;

    std::vector< trex::THandle<trex::TTPCHitPad> > hits = vol->GetHits();
    for(std::vector< trex::THandle<trex::TTPCHitPad> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){
      trex::THandle<trex::TTPCHitPad> hitPad = *hitIt;
      if(hitPad) used->AddHit(hitPad);
    };
  };
}
*/
void trex::TTPCTRExPatSubAlgorithm::GetTrackExtendedHits(std::vector<trex::TTPCVolGroup>& extHits){

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator trackIt = fTracks.begin(); trackIt != fTracks.end(); ++trackIt){
    trex::TTPCOrderedVolGroup& track = *trackIt;
    extHits.push_back(track.GetExtendedHits());
  }
  
}
