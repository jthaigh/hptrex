// eddy
#include "TTPCVolGroupMan.hxx"

ND::TTPCVolGroupMan::TTPCVolGroupMan(ND::TTPCLayout* layout){
  fHitMap = std::map<long, ND::TTPCUnitVolume*>();
  fLayout = layout;

  // set up primary groups of hits
  fPrimaryHits = ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));
}

void ND::TTPCVolGroupMan::AddPrimaryHits(std::map<long, ND::TTPCUnitVolume*> hitMap){
  fHitMap = std::map<long, ND::TTPCUnitVolume*>(hitMap);
  fPrimaryHits->AddHitMap(hitMap);
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GroupDeltaHits(ND::THandle<ND::TTPCVolGroup> suspectHits){
  // output
  std::vector< ND::THandle<ND::TTPCVolGroup> > groupedDeltas;
  // dummy input and full list of hits constructed
  ND::THandle<ND::TTPCVolGroup> dummyInput (new ND::TTPCVolGroup(fLayout));
  dummyInput->AddHits(suspectHits->begin(), suspectHits->end());
  ND::THandle<ND::TTPCVolGroup> dummyAllHits (new ND::TTPCVolGroup(fLayout));
  dummyAllHits->AddHits(fHitMap.begin(), fHitMap.end());

  // loop until all dummy hits are gone
  while(dummyInput->size()){
    ND::THandle<ND::TTPCVolGroup> groupedDelta = GrabADeltaHitGroup(dummyInput, dummyAllHits);
    groupedDeltas.push_back(groupedDelta);
  };

  return groupedDeltas;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GrabADeltaHitGroup(ND::THandle<ND::TTPCVolGroup> inputGroup, ND::THandle<ND::TTPCVolGroup> allHits){
  ND::THandle<ND::TTPCVolGroup> spreadGroup (new ND::TTPCVolGroup(fLayout));
  // initial distance
  double spreadRate = fLayout->GetDeltaSpreadRate();

  // cathode side
  int spreadGroupSense = (inputGroup->begin()->second->GetX() > 0) ? 1 : -1;

  // add first hit to begin with
  spreadGroup->AddHit(inputGroup->begin()->second, false);
  inputGroup->RemoveHit(inputGroup->begin()->first);

  // loop untill no more hits are added
  bool hitsAdded = true;
  double range = 0; // always save range and mid y and z for loop
  double midY = 0;
  double midZ = 0;
  while(hitsAdded){
    // expand in such a way that all hits are always encompassed
    hitsAdded = false;
    midY = (double)(spreadGroup->GetYMin() + spreadGroup->GetYMax()) / 2.;
    midZ = (double)(spreadGroup->GetZMin() + spreadGroup->GetZMax()) / 2.;

    double range2 = 0;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator addedIt = spreadGroup->begin(); addedIt != spreadGroup->end(); ++addedIt){
      double diffY = midY - (double)addedIt->second->GetY();
      double diffZ = midZ - (double)addedIt->second->GetZ();
      range2 = std::max(range2, (diffY*diffY + diffZ*diffZ));
    };
    range = std::sqrt(range2) + spreadRate;

    // add all hits inside the new range
    for(std::map<long, ND::TTPCUnitVolume*>::iterator addableIt = inputGroup->begin(); addableIt != inputGroup->end(); ++addableIt){
      // break if marked for deletion
      if(!addableIt->second) continue;

      // break at cathode boundry
      if( (addableIt->second->GetX() * spreadGroupSense) < 0 ) continue;

      // break if out of range
      double diffY = midY - (double)addableIt->second->GetY();
      double diffZ = midZ - (double)addableIt->second->GetZ();
      double newRange2 = diffY*diffY + diffZ*diffZ;

      if(newRange2 > range*range) continue; 

      // otherwise, add to the hit group and mark it for deletion from its parent
      spreadGroup->AddHit(addableIt->second, false);
      inputGroup->MarkHit(addableIt->first);
      hitsAdded = true;
    };
    // clear hits marked for deletion
    inputGroup->ClearMarked();
  };

  // now add all normal hits in a cylinder within range as well
  //double range = 0; // always save range and mid y and z for loop
  //double midY = 0;
  //double midZ = 0;
  double minX = (double)spreadGroup->GetXMin() - spreadRate;
  double maxX = (double)spreadGroup->GetXMax() + spreadRate;
  for(std::map<long, ND::TTPCUnitVolume*>::iterator addableIt = allHits->begin(); addableIt != allHits->end(); ++addableIt){
    // break if marked for deletion
    if(!addableIt->second) continue;

    // break at cathode boundry
    if( (addableIt->second->GetX() * spreadGroupSense) < 0 ) continue;

    // break if out of range
    double diffY = midY - (double)addableIt->second->GetY();
    double diffZ = midZ - (double)addableIt->second->GetZ();
    double newRange2 = diffY*diffY + diffZ*diffZ;

    if(newRange2 > range*range) continue; 

    // break if out of x range
    double addableX = (double)addableIt->second->GetX();
    if(addableX < minX || addableX > maxX) continue;

    // otherwise, add to the hit group and mark it for deletion from its parent
    spreadGroup->AddHit(addableIt->second, false);
    allHits->MarkHit(addableIt->first);
  };
  // clear hits marked for deletion
  allHits->ClearMarked();


  return spreadGroup;
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetAllEdgeGroups(){
  std::vector< ND::THandle<ND::TTPCVolGroup> > output;
  std::vector< ND::THandle<ND::TTPCVolGroup> > allGroups;
  std::vector< ND::THandle<ND::TTPCVolGroup> > deltaGroups = GetConnectedHits(ND::TTPCConnection::path, ND::TTPCHitGroupings::delta);
  std::vector< ND::THandle<ND::TTPCVolGroup> > nonDeltaGroups = GetConnectedHits(ND::TTPCConnection::path, ND::TTPCHitGroupings::nonDelta);

  // add the non-delta groups as one group
  ND::THandle<ND::TTPCVolGroup> nonDeltaGroup (new ND::TTPCVolGroup(fLayout));
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = nonDeltaGroups.begin(); grpIt != nonDeltaGroups.end(); ++grpIt) nonDeltaGroup->AddHits(*grpIt);
  allGroups.push_back(nonDeltaGroup);
  // add all of the delta groups individually
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = deltaGroups.begin(); grpIt != deltaGroups.end(); ++grpIt) allGroups.push_back(*grpIt);

  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator groupIt = allGroups.begin(); groupIt != allGroups.end(); ++groupIt){
    std::vector< ND::THandle<ND::TTPCVolGroup> > edgeGroups = GetEdgeGroups(*groupIt, true);
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator edgeIt = edgeGroups.begin(); edgeIt != edgeGroups.end(); ++edgeIt){
      ND::THandle<ND::TTPCVolGroup> edge = *edgeIt;
      output.push_back(edge);
    };
  };

  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpCheckIt = output.begin(); grpCheckIt != output.end(); ++grpCheckIt){
    ND::THandle<ND::TTPCVolGroup> grpCheck = *grpCheckIt;
  };

  return output;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetEdgeGroups(){
  return GetEdgeGroups(fPrimaryHits);
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetEdgeGroups(ND::THandle<ND::TTPCVolGroup> inGroup, bool fiddlyLeans){
  int layers = fLayout->GetEdgeLayers();
  bool indirect = fLayout->GetUseIndirectEdges();
  ND::TTPCConnection::Type type = ND::TTPCConnection::path;

  int xLeanOverride = 0;
  int yLeanOverride = 0;
  int zLeanOverride = 0;
  // special checks for groups that're too small to get an edge for
  if(inGroup->GetXSize() < 2*layers){
    if(inGroup->GetXMax() == fPrimaryHits->GetXMax()) xLeanOverride += 1;
    if(inGroup->GetXMin() == fPrimaryHits->GetXMin()) xLeanOverride -= 1;
  };
  if(inGroup->GetYSize() < 2*layers){
    if(inGroup->GetYMax() == fPrimaryHits->GetYMax()) yLeanOverride += 1;
    if(inGroup->GetYMin() == fPrimaryHits->GetYMin()) yLeanOverride -= 1;
  };
  if(inGroup->GetZSize() < 2*layers){
    if(inGroup->GetZMax() == fPrimaryHits->GetZMax()) zLeanOverride += 1;
    if(inGroup->GetZMin() == fPrimaryHits->GetZMin()) zLeanOverride -= 1;
  };

  // get groups of all hits on x, y and z edges
  ND::THandle<ND::TTPCVolGroup> edgeHitsXHi(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> edgeHitsXLo(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> edgeHitsYHi(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> edgeHitsYLo(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> edgeHitsZHi(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> edgeHitsZLo(new ND::TTPCVolGroup(fLayout));

  for(std::map<long, ND::TTPCUnitVolume*>::iterator vol = inGroup->begin(); vol != inGroup->end(); ++vol){
    // fetch x, y and z id and edge status for current cell
    int x = vol->second->GetX();
    int y = vol->second->GetY();
    int z = vol->second->GetZ();
    int edgeX = vol->second->GetEdgeX();
    int edgeY = vol->second->GetEdgeY();
    int edgeZ = vol->second->GetEdgeZ();

    // add hits to relevant group, in x only in the case of deltas 
    if((x > (inGroup->GetXMax() - layers) ) && (indirect || (edgeX==1 )) ) edgeHitsXHi->AddHit(vol->second);
    if((x < (inGroup->GetXMin() + layers) ) && (indirect || (edgeX==-1)) ) edgeHitsXLo->AddHit(vol->second);
    if(!fLayout->GetUseAltEdgeDetection() || !vol->second->GetDeltaTagged()){
      if((y > (inGroup->GetYMax() - layers) ) && (indirect || (edgeY==1 )) ) edgeHitsYHi->AddHit(vol->second);
      if((y < (inGroup->GetYMin() + layers) ) && (indirect || (edgeY==-1)) ) edgeHitsYLo->AddHit(vol->second);
      if((z > (inGroup->GetZMax() - layers) ) && (indirect || (edgeZ==1 )) ) edgeHitsZHi->AddHit(vol->second);
      if((z < (inGroup->GetZMin() + layers) ) && (indirect || (edgeZ==-1)) ) edgeHitsZLo->AddHit(vol->second);
    };
  };

  std::vector< ND::THandle<ND::TTPCVolGroup> > edgeGroups;
  // recursively build groups of hits on x, y and z edges
  if(fiddlyLeans){
    FillWithSplitGroups(edgeGroups, edgeHitsXHi, type, 1*(inGroup->GetXMax() == fPrimaryHits->GetXMax()) ,0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsXLo, type, -1*(inGroup->GetXMin() == fPrimaryHits->GetXMin()),0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYHi, type, 0,1*(inGroup->GetYMax() == fPrimaryHits->GetYMax()) ,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYLo, type, 0,-1*(inGroup->GetYMin() == fPrimaryHits->GetYMin()),0);
    FillWithSplitGroups(edgeGroups, edgeHitsZHi, type, 0,0,1*(inGroup->GetZMax() == fPrimaryHits->GetZMax()) );
    FillWithSplitGroups(edgeGroups, edgeHitsZLo, type, 0,0,-1*(inGroup->GetZMin() == fPrimaryHits->GetZMin()));
  }
  else{
    FillWithSplitGroups(edgeGroups, edgeHitsXHi, type, 1 ,0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsXLo, type, -1,0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYHi, type, 0,1 ,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYLo, type, 0,-1,0);
    FillWithSplitGroups(edgeGroups, edgeHitsZHi, type, 0,0,1 );
    FillWithSplitGroups(edgeGroups, edgeHitsZLo, type, 0,0,-1);
  };

  // loop over all found groups on x, y and z edges
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp1It = edgeGroups.begin(); grp1It != edgeGroups.end(); ++grp1It){
    ND::THandle<ND::TTPCVolGroup> grp1 = *grp1It;
    // if group overlaps with another, delete the biggest (break ties in favour of z, then y)
    // smallest group is most likely to be close to end, of e.g. low angle track
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp2It = grp1It+1; grp2It != edgeGroups.end(); ++grp2It){
      ND::THandle<ND::TTPCVolGroup> grp2 = *grp2It;
      if(grp2->empty()) continue;
      for(std::map<long, ND::TTPCUnitVolume*>::iterator grp2El = grp2->begin(); grp2El != grp2->end(); ++grp2El){
        if(grp1->Contains(grp2El->first)){
          if(grp2->size() < grp1->size()){
            grp1->SetXLean(grp2->GetXLean());
            grp1->SetYLean(grp2->GetYLean());
            grp1->SetZLean(grp2->GetZLean());
            grp1->clear();
          }
          else{
            grp2->SetXLean(grp1->GetXLean());
            grp2->SetYLean(grp1->GetYLean());
            grp2->SetZLean(grp1->GetZLean());
            grp2->clear();
          };
          break;
        };
      };
    };
  };
  // clean up empty groups
  std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpDel = edgeGroups.begin();
  while(grpDel != edgeGroups.end()){
    if((*grpDel)->empty()) edgeGroups.erase(grpDel);
    else grpDel ++;
  };
  // override leans
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = edgeGroups.begin(); grpIt != edgeGroups.end(); ++grpIt){
    ND::THandle<ND::TTPCVolGroup> grp = *grpIt;
    if(xLeanOverride) grp->SetXLean(xLeanOverride, true);
    if(yLeanOverride) grp->SetYLean(yLeanOverride, true);
    if(zLeanOverride) grp->SetZLean(zLeanOverride, true);
  };
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp1It = edgeGroups.begin(); grp1It != edgeGroups.end(); ++grp1It){
    ND::THandle<ND::TTPCVolGroup> grp1 = *grp1It;
  };
  return edgeGroups;
}
void ND::TTPCVolGroupMan::FillWithSplitGroups(std::vector< ND::THandle<ND::TTPCVolGroup> >& container, ND::THandle<ND::TTPCVolGroup> inputHits, ND::TTPCConnection::Type type, int maxFilterX, int maxFilterY, int maxFilterZ){
  std::vector< ND::THandle<ND::TTPCVolGroup> > splitGroups = GetSplitGroups(inputHits, type);
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp = splitGroups.begin(); grp != splitGroups.end(); ++grp){
    if(maxFilterX) (*grp)->SetXLean(maxFilterX);
    if(maxFilterY) (*grp)->SetYLean(maxFilterY);
    if(maxFilterZ) (*grp)->SetZLean(maxFilterZ);

    container.push_back(*grp);
  };
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetSplitGroups(ND::THandle<ND::TTPCVolGroup> inputHits, ND::TTPCConnection::Type type){
  std::vector< ND::THandle<ND::TTPCVolGroup> > groups = std::vector< ND::THandle<ND::TTPCVolGroup> >();

  for(int i=0; i<9999; i++){
    if(inputHits->size() < 1) break;
    std::map<long, ND::TTPCUnitVolume*>::iterator el = inputHits->begin();

    ND::THandle<ND::TTPCVolGroup> newGroup(new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()));
    RecursiveFriendBuild(el->first, newGroup, inputHits, type);

    if(newGroup->size() > 0) groups.push_back(newGroup);
  };
  return groups;
}

ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetFocus(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // get focus with arbitrary threshold close to unity
  std::vector< ND::THandle<ND::TTPCVolGroup> > foci = GetFoci(paths, .99);
  // return first (should be only) element
  if(foci.size() > 0) return *foci.begin();
  else return ND::THandle<ND::TTPCVolGroup>(new TTPCVolGroup(fLayout));
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetFoci(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, float diffThreshold){
  // select pairs of groups to loop over
  std::vector< std::pair< ND::THandle<ND::TTPCOrderedVolGroup>, ND::THandle<ND::TTPCOrderedVolGroup> > > groupPairs;
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It){
    ND::THandle<ND::TTPCOrderedVolGroup> path1 = *path1It;
    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
      ND::THandle<ND::TTPCOrderedVolGroup> path2 = *path2It;

      // make sure paths are populated
      if(path1->size() < 1 || path2->size() < 1) continue;

      // make sure paths start in same place
      ND::TTPCUnitVolume* beginVol1 = (*(path1->begin()))->GetUnitVolume();
      ND::TTPCUnitVolume* beginVol2 = (*(path2->begin()))->GetUnitVolume();
      if(beginVol1 != beginVol2) continue;

      bool startAdded = false;
      // loop and make sure only the best pair is chosen for each start point
      for(std::vector< std::pair< ND::THandle<ND::TTPCOrderedVolGroup>, ND::THandle<ND::TTPCOrderedVolGroup> > >::iterator groupPairIt = groupPairs.begin(); groupPairIt != groupPairs.end(); ++groupPairIt){
        ND::THandle<ND::TTPCOrderedVolGroup> pairGrp1 = groupPairIt->first;
        ND::THandle<ND::TTPCOrderedVolGroup> pairGrp2 = groupPairIt->second;

        // do they start in the same place?
        if( (*(pairGrp1->begin()))->GetUnitVolume() != beginVol1) continue;
        startAdded = true;

        // end volumes
        ND::TTPCUnitVolume* endVol1 = (*(path1->rbegin()))->GetUnitVolume();
        ND::TTPCUnitVolume* endVol2 = (*(path2->rbegin()))->GetUnitVolume();
        ND::TTPCUnitVolume* endVol3 = (*(pairGrp1->rbegin()))->GetUnitVolume();
        ND::TTPCUnitVolume* endVol4 = (*(pairGrp2->rbegin()))->GetUnitVolume();

        // compare triangle areas to work out best one
        double area1 = GetTriangleArea(beginVol1, endVol1, endVol2);
        double area2 = GetTriangleArea(beginVol1, endVol3, endVol4);

        // replace old with new if new gives a bigger triangle
        if(area1 > area2){
          *groupPairIt = std::pair< ND::THandle<ND::TTPCOrderedVolGroup>, ND::THandle<ND::TTPCOrderedVolGroup> >(path1, path2);
        };
      };

      // push back if start not already found
      if(!startAdded) groupPairs.push_back( std::pair<ND::THandle<ND::TTPCOrderedVolGroup>, ND::THandle<ND::TTPCOrderedVolGroup> >(path1, path2) );
    };
  };

  // groups for holding hits at points of divergence
  ND::THandle<ND::TTPCVolGroup> grp(new ND::TTPCVolGroup(fLayout));
  // map for counting how many groups average at a given point
  std::map<long, int> avgPoints = std::map<long, int>();
  for(std::vector< std::pair< ND::THandle<ND::TTPCOrderedVolGroup>, ND::THandle<ND::TTPCOrderedVolGroup> > >::iterator groupPairIt = groupPairs.begin(); groupPairIt != groupPairs.end(); ++groupPairIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path1 = groupPairIt->first;
    ND::THandle<ND::TTPCOrderedVolGroup> path2 = groupPairIt->second;

    // get detatchment points
    ND::TTPCPathVolume* detatch1 = GetDetatchmentPoint(path1, path2, ND::TTPCConnection::vertexFind);
    ND::TTPCPathVolume* detatch2 = GetDetatchmentPoint(path2, path1, ND::TTPCConnection::vertexFind);

    // find the closest points between the two lines
    ND::TTPCUnitVolume* closestPoint = GetZero(path1, path2, detatch1, detatch2);

    // get closest point from projection points (defunct)
    //ND::TTPCPathVolume* project1 = GetProjectionPoint(path1, detatch1, fLayout->GetExtendedPathDist());
    //ND::TTPCPathVolume* project2 = GetProjectionPoint(path2, detatch2, fLayout->GetExtendedPathDist());
    //ND::TTPCUnitVolume* closestPoint = GetClosestPoint(project1, project2, detatch1, detatch2);

    // add those points and save number in each volume for future calculations
    if(closestPoint){
      grp->AddHit(closestPoint);
      if(!avgPoints.count(closestPoint->GetID())) avgPoints[closestPoint->GetID()] = 1;
      else avgPoints[closestPoint->GetID()]++;
    };
  };

  // break points of divergence into separate groups
  std::vector< ND::THandle<ND::TTPCVolGroup> > groups = GetConnectedHits(grp, ND::TTPCConnection::vertexMerge);

  int maxSize = -999999;
  // find group of biggest size
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = groups.begin(); grpIt != groups.end(); ++grpIt){
    maxSize = std::max(maxSize, int((*grpIt)->size()));
  };
  // set cutoff for group size at fraction of maxSize defined by threshold
  float minSize = float(maxSize) * diffThreshold;

  // build list of groups of size above the threshold
  std::vector< ND::THandle<ND::TTPCVolGroup> > groupsOut;
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grpIt = groups.begin(); grpIt != groups.end(); ++grpIt){
    ND::THandle<ND::TTPCVolGroup> grp = *grpIt;
    // find weighted average position of hit in group
    if(float(grp->size()) >= minSize){
      ND::THandle<ND::TTPCVolGroup> newGrp (new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()));
      int avgX = 0;
      int avgY = 0;
      int avgZ = 0;
      int weightSum = 0;
      for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = grp->begin(); volEl != grp->end(); ++volEl){
        int weight = avgPoints[volEl->first];

        avgX += weight * volEl->second->GetX();
        avgY += weight * volEl->second->GetY();
        avgZ += weight * volEl->second->GetZ();
        weightSum += weight;
      };
      avgX = int( (float)avgX/(float)weightSum + .5);
      avgY = int( (float)avgY/(float)weightSum + .5);
      avgZ = int( (float)avgZ/(float)weightSum + .5);

      // add average hit
      long id = grp->GetNearestHit(avgX, avgY, avgZ);
      newGrp->AddHit(grp->GetEl(id));

      // also add nearby hits
      ND::THandle<ND::TTPCVolGroup> nearHits = GetNearHits(fPrimaryHits, id, ND::TTPCConnection::vertexHit);
      newGrp->AddHits(nearHits);

      // close and add to main group
      groupsOut.push_back(newGrp);
    };
  };

  // return those groups
  return groupsOut;
}

ND::TTPCUnitVolume* ND::TTPCVolGroupMan::GetZero(ND::THandle<ND::TTPCOrderedVolGroup> path1, ND::THandle<ND::TTPCOrderedVolGroup> path2, ND::TTPCPathVolume* vol1, ND::TTPCPathVolume* vol2){
  ND::TTPCPathVolume* nearVol1 = GetPathZero(path1, path2, vol1);
  ND::TTPCPathVolume* nearVol2 = GetPathZero(path2, path1, vol2);

  if(nearVol1 && nearVol2){
    // average
    int ax = (int)((nearVol2->GetX() + nearVol1->GetX())*.5);
    int ay = (int)((nearVol2->GetY() + nearVol1->GetY())*.5);
    int az = (int)((nearVol2->GetZ() + nearVol1->GetZ())*.5);

    long id = fPrimaryHits->GetNearestHit(ax, ay, az);
    return fPrimaryHits->GetHit(id);
  }
  else if(nearVol1){
    return nearVol1->GetUnitVolume();
  }
  else if(nearVol2){
    return nearVol2->GetUnitVolume();
  }
  else if(vol1 && vol2){
    // average
    int ax = (int)((vol2->GetX() + vol1->GetX())*.5);
    int ay = (int)((vol2->GetY() + vol1->GetY())*.5);
    int az = (int)((vol2->GetZ() + vol1->GetZ())*.5);

    long id = fPrimaryHits->GetNearestHit(ax, ay, az);
    return fPrimaryHits->GetHit(id);
  }
  else if(vol1){
    return vol1->GetUnitVolume();
  }
  else if(vol2){
    return vol2->GetUnitVolume();
  }
  else{
    return 0;
  };
}
ND::TTPCPathVolume* ND::TTPCVolGroupMan::GetPathZero(ND::THandle<ND::TTPCOrderedVolGroup> path1, ND::THandle<ND::TTPCOrderedVolGroup> path2, ND::TTPCPathVolume* vol){
  if(!vol) vol = *(path1->rbegin());

  // start checking after first points fed in
  bool go=false;

  // initial distances and number of iterations
  double initDist = -1.;
  int nIterations = 0;

  ND::TTPCPathVolume* nearVol = 0;

  // check backwards along first path
  for(std::vector<ND::TTPCPathVolume*>::reverse_iterator pathVolRit = path1->rbegin(); pathVolRit != path1->rend(); ++pathVolRit){
    ND::TTPCPathVolume* pathVol = *pathVolRit;
    if(go){
      if(nIterations<1){
        initDist = GetMinDistance(path2, pathVol);
      }
      else{
        double curDist = GetMinDistance(path2, pathVol);
        // if more than half way to the zero, project
        if( (curDist/initDist) / .5){
          while(nIterations > 0 && pathVolRit != (path1->rend()-1)){
            nIterations--;
            pathVolRit++;
          };
          nearVol = *pathVolRit;
          break;
        };
      };
      nIterations++;
    }
    else go = (pathVol==vol);
  };

  return nearVol;
}
double ND::TTPCVolGroupMan::GetTriangleArea(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, ND::TTPCUnitVolume* vol3){
  TVector3 sideA = vol2->GetPos() - vol1->GetPos();
  TVector3 sideB = vol3->GetPos() - vol1->GetPos();

  double a2 = sideA.Mag2();
  double b2 = sideB.Mag2();
  double ab = std::sqrt(a2*b2);
  double sinC = sideA.Angle(sideB);
  return .5*ab*sinC;
}
double ND::TTPCVolGroupMan::GetMinDistance(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCPathVolume* point){
  return GetMinDistance(path, point->GetUnitVolume());
}
double ND::TTPCVolGroupMan::GetMinDistance(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCUnitVolume* vol){
  double minDist2 = 99999999.;
  for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
    ND::TTPCPathVolume* pathVol = *pathVolIt;

    double dist2 = (pathVol->GetPosXYZ() - vol->GetPosXYZ()).Mag2();
    minDist2 = std::min(minDist2, dist2);
  };
  return std::sqrt(minDist2);
}
ND::TTPCPathVolume* ND::TTPCVolGroupMan::GetDetatchmentPoint(ND::THandle<ND::TTPCOrderedVolGroup> path1, ND::THandle<ND::TTPCOrderedVolGroup> path2, ND::TTPCConnection::Type type){
  int sizeX;
  int sizeY;
  int sizeZ;
  fLayout->GetTypeDistances(sizeX, sizeY, sizeZ, ND::TTPCConnection::vertexFind);

  // check most likely point first to save time
  int prevPoint=0;
  int prevMax=path2->size();

  for(std::vector<ND::TTPCPathVolume*>::iterator point1It = path1->begin(); point1It != path1->end(); ++point1It){
    bool wasFound = false;

    ND::TTPCPathVolume* point1 = *point1It;
    ND::TTPCPathVolume* point2;

    int prevOld=prevPoint;
    int prevOff=0;

    // check same index first to save time, since this should usually be true if anything is
    while(prevOld+prevOff < prevMax || prevOld-prevOff >= 0){
      for (int sign=-1; sign<=1; sign+=2){
        int dummy = prevOld+sign*prevOff;
        if(dummy < 0) continue;
        if(dummy >= prevMax) continue;
        prevPoint = dummy;
        point2 = (*path2)[prevPoint];
        if(IsInRange(point1, point2, sizeX, sizeY, sizeZ)){
          wasFound = true;
          break;
        };
      };
      if(wasFound) break;

      prevOff++;
    };
    if(!wasFound) return point1;
  };
  return *(path1->end()-1);
}
ND::TTPCPathVolume* ND::TTPCVolGroupMan::GetProjectionPoint(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCPathVolume* pathVol, int dist){
  // get projection point
  bool foundDetatch = false;
  for(std::vector<ND::TTPCPathVolume*>::reverse_iterator pointRit = path->rbegin(); pointRit != path->rend(); ++pointRit){
    ND::TTPCPathVolume* point=*pointRit;
    if(!foundDetatch){
      if(point != pathVol) continue;
      foundDetatch = true;
    }
    else{
      int dx = point->GetX() - pathVol->GetX();
      int dy = point->GetY() - pathVol->GetY();
      int dz = point->GetZ() - pathVol->GetZ();
      if(dx*dx + dy*dy + dz*dz > dist*dist) return point;
    };
  };
  return *(path->rend()-1);
}
ND::TTPCUnitVolume* ND::TTPCVolGroupMan::GetClosestPoint(ND::TTPCPathVolume* begin1, ND::TTPCPathVolume* end1, ND::TTPCPathVolume* begin2, ND::TTPCPathVolume* end2){
  // closest position
  int ax;
  int ay;
  int az;
  // set up start points and directions
  TVector3 path1_0 = begin1->GetPos();
  TVector3 path2_0 = begin2->GetPos();
  TVector3 path1_1 = end1->GetPos();
  TVector3 path2_1 = end2->GetPos();
  TVector3 path1_d = path1_1 - path1_0;
  TVector3 path2_d = path2_1 - path2_0;
  TVector3 begin_d = path1_0 - path2_0;

  // vector between lines is minimum when (min_v . path1_d) = 0 and (min_v . path2_d) = 0
  double a = path1_d.Mag2();
  double b = path1_d.Dot(path2_d);
  double c = path2_d.Mag2();
  double d = path1_0.Dot(begin_d);
  double e = path2_0.Dot(begin_d);

  double disc = (a*c - b*b);
  if(std::abs(disc) < .2){
    // basically parallel - just use average of end points
    ax = (int)((end1->GetX() + end2->GetX())*.5);
    ay = (int)((end1->GetY() + end2->GetY())*.5);
    az = (int)((end1->GetZ() + end2->GetZ())*.5);
    return 0;
  }
  else{
    double s = b*e - c*d / disc;
    double t = a*e - b*d / disc;

    if(s < 0 || t < 0){
      // require s and t to be positive - otherwise just use average of start points
      ax = (int)((begin1->GetX() + begin2->GetX())*.5);
      ay = (int)((begin1->GetY() + begin2->GetY())*.5);
      az = (int)((begin1->GetZ() + begin2->GetZ())*.5);
      return 0;
    }
    else{
      // average of closest points along two lines
      ax = (int)( ( begin1->GetX() + s*(end1->GetX()-begin1->GetX()) + begin2->GetX() + t*(end2->GetX()-begin2->GetX()) )*.5);
      ay = (int)( ( begin1->GetY() + s*(end1->GetY()-begin1->GetY()) + begin2->GetY() + t*(end2->GetY()-begin2->GetY()) )*.5);
      az = (int)( ( begin1->GetZ() + s*(end1->GetZ()-begin1->GetZ()) + begin2->GetZ() + t*(end2->GetZ()-begin2->GetZ()) )*.5);
    };
  };

  long id = fPrimaryHits->GetNearestHit(ax, ay, az);
  return fPrimaryHits->GetHit(id);
}
void ND::TTPCVolGroupMan::CleanUpVertices(std::vector< ND::THandle<ND::TTPCVolGroup> > edgeGroups, std::vector< ND::THandle<ND::TTPCVolGroup> > vertices){
  while(vertices.size() > edgeGroups.size()-2){
    // get rid of spurious vertices, one at a time
    std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator markedIt;
    double minDistance = 99999999.;

    // find closest vertex to any edge group and mark it as spurious
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
      ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;

      if(vertex)
      for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator edgeIt = edgeGroups.begin(); edgeIt != edgeGroups.end(); ++edgeIt){
        ND::THandle<ND::TTPCVolGroup> edge = *edgeIt;
        TVector3 diff = edge->GetAveragePosXYZ() - vertex->GetAveragePosXYZ(); 
        double dist = diff.Mag2();

        if(dist < minDistance){
          minDistance = dist;
          markedIt = vertexIt;
        };
      };
    };

    vertices.erase(markedIt);
  };
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetFarHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, bool tryChargeCut){
  ND::THandle<ND::TTPCVolGroup> inGroup (new ND::TTPCVolGroup(fLayout));

  // add all far hits from all paths to group
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path = paths.begin(); path != paths.end(); ++path){
    ND::THandle<ND::TTPCVolGroup> farHitsGroup;
    farHitsGroup = GetFarHitsPreGroup(*path, tryChargeCut);
    if(!farHitsGroup->size()){
      farHitsGroup = GetFarHitsGroup(*path, tryChargeCut);
    };
    inGroup->AddHits(farHitsGroup);
  };
  std::vector< ND::THandle<ND::TTPCVolGroup> > groups = GetConnectedHits(inGroup, ND::TTPCConnection::vertexMerge);

  return groups;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetFarHitsPreGroup(ND::THandle<ND::TTPCOrderedVolGroup> path, bool tryChargeCut){
  double edgePreDist = fLayout->GetEdgePreDist();
  double edgePreDist2 = edgePreDist*edgePreDist;
  double edgePreAng = fLayout->GetEdgePreAng();

  ND::THandle<ND::TTPCVolGroup> outGroup (new ND::TTPCVolGroup(fLayout));

  ND::TTPCPathVolume* kinkPoint = 0;
  double kinkAngle = 360.;

  // look for local kinks around each point
  for(int i=0; i<path->size(); ++i){
    ND::TTPCPathVolume* pathVol = path->at(i);

    // look for back and forward check points
    ND::TTPCPathVolume* backCheck = 0;
    ND::TTPCPathVolume* frontCheck = 0;
    for(int j=i; j>=0; --j){
      ND::TTPCPathVolume* backCheckVol = path->at(j);
      TVector3 posDiff = backCheckVol->GetPosXYZ() - pathVol->GetPosXYZ();

      if(posDiff.Mag2() >= edgePreDist2){
        backCheck = backCheckVol;
        break;
      };
    };
    for(int j=i; j<path->size(); ++j){
      ND::TTPCPathVolume* frontCheckVol = path->at(j);
      TVector3 posDiff = frontCheckVol->GetPosXYZ() - pathVol->GetPosXYZ();

      if(posDiff.Mag2() >= edgePreDist2){
        frontCheck = frontCheckVol;
        break;
      };
    };

    // now do local angle check
    if(backCheck && frontCheck){
      TVector3 toBack = backCheck->GetPosXYZ() - pathVol->GetPosXYZ();
      TVector3 toFront = frontCheck->GetPosXYZ() - pathVol->GetPosXYZ();

      //double angle = toBack.Angle(toFront) * 180. / TMath::Pi();
      // checking for edges caused by delta spirals, so only need to check in yz
      double dotProd = toFront.Y()*toBack.Y() + toFront.Z()*toBack.Z();
      double magProd = std::sqrt( (toFront.Y()*toFront.Y() + toFront.Z()*toFront.Z()) * (toBack.Y()*toBack.Y() + toBack.Z()*toBack.Z()) );
      double angle = std::acos( dotProd / magProd ) * 180./TMath::Pi();

      if(angle < edgePreAng){
        if(angle < kinkAngle){
          kinkPoint = pathVol;
          kinkAngle = angle;
        };
      };
    };
  };

  // add kink point to group
  // TODO: maybe add things behind it as well?  could get messy
  if(kinkPoint){
    outGroup->AddHit(kinkPoint->GetUnitVolume());
  };

  return outGroup;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetFarHitsGroup(ND::THandle<ND::TTPCOrderedVolGroup> path, bool tryChargeCut){
  ND::TTPCPathVolume* startVol = *(path->begin());
  ND::TTPCPathVolume* endVol = *(path->rbegin());

  TVector3 startPos = startVol->GetPosXYZ();
  TVector3 endPos = endVol->GetPosXYZ();
  TVector3 parNorm = (startPos - endPos).Unit();
  float parPos = (startPos - endPos).Mag();

  // build group to look along when finding kinks
  ND::THandle<ND::TTPCVolGroup> inGroup = ND::THandle<ND::TTPCVolGroup> (new ND::TTPCVolGroup(fLayout));
  for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
    ND::TTPCPathVolume* pathVol = *pathVolIt;
    ND::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

    inGroup->AddHit(vol);
  };
  // build group to look along when finding maximum hit
  ND::THandle<ND::TTPCVolGroup> extendedGroup;
  if(tryChargeCut){
    extendedGroup = ND::THandle<ND::TTPCVolGroup> (new ND::TTPCVolGroup(fLayout));

    ND::THandle<ND::TTPCVolGroup> extendedGroupAll = path->GetExtendedHits();
    for(std::map<long, ND::TTPCUnitVolume*>::iterator el = extendedGroupAll->begin(); el != extendedGroupAll->end(); ++el){
      if(!el->second->GetPathology()){
        extendedGroup->AddHit(el->second);
      };
    };
  }
  else{
    extendedGroup = path->GetExtendedHits();
  };

  // find single path max dist
  float maxDist2 = -1.;

  for(std::map<long, ND::TTPCUnitVolume*>::iterator el = inGroup->begin(); el != inGroup->end(); ++el){
    TVector3 diff = startPos - el->second->GetPosXYZ();
    TVector3 parallel = parNorm*diff.Dot(parNorm);
    TVector3 dist = diff - parallel;

    float dist2 = dist.Mag2();
    if(dist2 > maxDist2){
      maxDist2 = dist2;
    };
  };

  if (maxDist2 < 0.) return ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));
  float threshDist2 = maxDist2 * fLayout->GetEdgeRange()*fLayout->GetEdgeRange();
  float maxDist = std::sqrt(maxDist2);

  // find extended group max dist
  float extendedMaxDist2 = -1.;
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el = extendedGroup->begin(); el != extendedGroup->end(); ++el){
    TVector3 diff = startPos - el->second->GetPosXYZ();
    TVector3 parallel = parNorm*diff.Dot(parNorm);
    TVector3 dist = diff - parallel;

    float dist2 = dist.Mag2();
    if(dist2 > extendedMaxDist2){
      extendedMaxDist2 = dist2;
    };
  };

  float extendedThreshDist2 = extendedMaxDist2 * fLayout->GetEdgeRange()*fLayout->GetEdgeRange();

  // count fraction of hits above threshold
  int nAbove = 0;
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el = inGroup->begin(); el != inGroup->end(); ++el){
    TVector3 diff = startPos - el->second->GetPosXYZ();
    TVector3 parallel = parNorm*diff.Dot(parNorm);
    TVector3 dist = diff - parallel;

    double distMag2 = dist.Mag2();
    if ((float)distMag2 > threshDist2){
      nAbove ++;
    };
  };

  float threshFrac = (float)nAbove / (float)inGroup->size();
  if(!threshFrac) return ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));

  // vector for parallel and perpendicular positions of all hits above extended threshold
  std::vector<float> parPositions;
  std::vector<float> perpPositions;
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el = extendedGroup->begin(); el != extendedGroup->end(); ++el){
    TVector3 diff = startPos - el->second->GetPosXYZ();
    TVector3 parallel = parNorm*diff.Dot(parNorm);
    TVector3 dist = diff - parallel;

    double distMag2 = dist.Mag2();
    if ((float)distMag2 > extendedThreshDist2){
      parPositions.push_back( parallel.Mag()/parPos );
      perpPositions.push_back( dist.Mag()/maxDist );
    };
  };

  // calculate statistical mean and spread in positions above threshold, weighted by perpendicular distance
  float threshMean;
  float threshSigma;

  float norm = 0;
  float threshSum = 0;
  for(unsigned int i=0; i<parPositions.size(); ++i){
    float parDist = parPositions.at(i);
    float perpDist = perpPositions.at(i);
    perpDist = 1.;

    threshSum += parDist*perpDist;
    norm += perpDist;
  };
  threshMean = threshSum / norm;

  float threshVarSum = 0;
  int nonZeroWeights = 0;
  for(unsigned int i=0; i<parPositions.size(); ++i){
    float parDist = parPositions.at(i);
    float perpDist = perpPositions.at(i);
    perpDist = 1.;

    float threshDiff = parDist - threshMean;

    threshVarSum += threshDiff*threshDiff * perpDist;
    nonZeroWeights += int(perpDist > 0.);
  };
  threshSigma = std::sqrt(  threshVarSum / ( (nonZeroWeights - 1) * norm / nonZeroWeights )  );

  bool edge = false;

  // check one - standard kinks
  // ensure max hit is above max dist
  if(maxDist > fLayout->GetEdgeMinDist()){
    // ensure fraction above threshold is above minimum
    if(threshFrac > fLayout->GetEdgeMinHits()){
      // ensure fraction of hits above threshold is above minimum
      if(threshFrac < fLayout->GetEdgeThreshold()){
        edge = true;
      };
    };
  };
  // check two - small kinks at the end of short tracks
  // ensure max hit is above low max dist
  if(maxDist > fLayout->GetEdgeMinDistLow()){
    // ensure fraction above threshold is above minimum
    if(threshFrac > fLayout->GetEdgeMinHits()){
      // ensure sigma is low enough for this not to just be a straight path
      if(threshSigma < fLayout->GetEdgeMaxSigma()){
        // ensure maximum hit is offset suitably
        if(threshMean < fLayout->GetEdgeOffset() || threshMean > (1.-fLayout->GetEdgeOffset())){
          edge = true;
        };
      };
    };
  };

  if(edge){
    ND::THandle<ND::TTPCVolGroup> outGroup (new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()) );

    ND::TTPCUnitVolume* furthestHit = 0;
    double furthestDist = -1;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator el = extendedGroup->begin(); el != extendedGroup->end(); ++el){
      TVector3 diff = startPos - el->second->GetPosXYZ();
      TVector3 dist = diff - (parNorm*diff.Dot(parNorm));

      double distMag2 = dist.Mag2();
      if (distMag2 > furthestDist){
        furthestHit = el->second;
        furthestDist = distMag2;
      };
    };

    // add furthest hit
    if(furthestHit){
      outGroup->AddHit(furthestHit);
    };

    return outGroup;
  };

  // default
  return ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetDiscontinuity(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // get focus with arbitrary threshold close to unity
  std::vector< ND::THandle<ND::TTPCVolGroup> > discontinuities = GetDiscontinuities(paths, .99, 2.);
  // return first (should be only) element
  if(discontinuities.size() > 0) return *discontinuities.begin();
  else return ND::THandle<ND::TTPCVolGroup>(new TTPCVolGroup(fLayout));
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetDiscontinuities(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, float diffThreshold, float threshold){
  ND::THandle<ND::TTPCVolGroup> inGroup (new ND::TTPCVolGroup(fLayout));

  // add all discontinuity hits from all paths to group
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path = paths.begin(); path != paths.end(); ++path){
    int discX = GetDiscontinuityID(*path, 1, threshold);
    int discY = GetDiscontinuityID(*path, 2, threshold);
    int discZ = GetDiscontinuityID(*path, 3, threshold);

    if (discX > 0) inGroup->AddHit((**path)[discX]->GetUnitVolume());
    if (discY > 0) inGroup->AddHit((**path)[discY]->GetUnitVolume());
    if (discZ > 0) inGroup->AddHit((**path)[discZ]->GetUnitVolume());
  };

  std::vector< ND::THandle<ND::TTPCVolGroup> > groups = GetConnectedHits(inGroup, ND::TTPCConnection::path);

  int maxSize = -999999;
  // find group of biggest size
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp = groups.begin(); grp != groups.end(); ++grp){
    maxSize = std::max(maxSize, int((*grp)->size()));
  };
  // set cutoff for group size at fraction of maxSize defined by threshold
  float minSize = float(maxSize) * diffThreshold;

  // build list of groups of size above the threshold
  std::vector< ND::THandle<ND::TTPCVolGroup> > groupsOut;
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp = groups.begin(); grp != groups.end(); ++grp){
    if(float((*grp)->size()) >= minSize){
      groupsOut.push_back(*grp);
    };
  };

  return groupsOut;
}
int ND::TTPCVolGroupMan::GetDiscontinuityID(ND::THandle<ND::TTPCOrderedVolGroup> path, int dimension, float threshold){
  // set up histograms for x, y and z differences
  float size = float(path->size()-1);

  TH1F* diffs = new TH1F("diffs", "Diffs", size, 1., size);

  // loop over all points in path
  int i=1;
  for(std::vector<ND::TTPCPathVolume*>::iterator it = path->begin(); it != path->end()-1; ++it){
    // get difference of positions
    TVector3 diff = (*(it+1))->GetPos() - (*(it))->GetPos();

    if(dimension == 1){
      diffs->SetBinContent(i, diff.X());
      diffs->SetBinError(i, fLayout->GetXCellSize()/2.);
    }
    else if(dimension == 2){
      diffs->SetBinContent(i, diff.Y());
      diffs->SetBinError(i, fLayout->GetPadPitchY()/2.);
    }
    else if(dimension == 3){
      diffs->SetBinContent(i, diff.Z());
      diffs->SetBinError(i, fLayout->GetPadPitchZ()/2.);
    };
    i++;
  };
  float disPos;
  float disSize;

  FindDiscontinuity(disPos, disSize, size, diffs, threshold);

  delete diffs;
  if (disSize > 0.){
    return int(disPos); 
  };
  return 0;
}
void ND::TTPCVolGroupMan::FindDiscontinuity(float& pos, float& step, float size, TH1F* hist, float threshold){
  TF1* fitLine = new TF1("line", "[0]", 1., size);
  TF1* fitStep = new TF1("step", "[0] + [1]*(x<[2])", 1., size);

  hist->Fit(fitLine, "q");

  fitStep->SetParameters(fitLine->GetParameter(0), 0., size/2.);
  fitStep->SetParameter(0, fitLine->GetParameter(0));

  hist->Fit(fitStep, "q");

  bool isStep = (fitLine->GetChisquare() / fitStep->GetChisquare()) > threshold;
  if(isStep){
    pos = fitStep->GetParameter(2);
    step = fitStep->GetParameter(1);
  }
  else{
    pos = 0.;
    step = 0.;
  };
  delete fitLine;
  delete fitStep;
}

std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ND::TTPCVolGroupMan::BreakPathsAboutKinks(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > inPaths, bool tryChargeCut){
  std::vector< ND::THandle<ND::TTPCVolGroup> > inVertices = GetJunctionsFromPaths(inPaths);
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > outPaths = std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >();

  // loop over all input paths
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = inPaths.begin(); pathIt != inPaths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

    ND::THandle<ND::TTPCOrderedVolGroup> outPath1 = ND::THandle<ND::TTPCOrderedVolGroup> (new ND::TTPCOrderedVolGroup(fLayout));
    ND::THandle<ND::TTPCOrderedVolGroup> outPath2 = ND::THandle<ND::TTPCOrderedVolGroup> (new ND::TTPCOrderedVolGroup(fLayout));

    ND::THandle<ND::TTPCVolGroup> kinkGroup;
    if(fLayout->GetUseAltEdgeDetection() && path->GetDeltaCriteriaMet()){
      // don't try kink finding if path is a delta
      kinkGroup = ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));
    }
    else{
      kinkGroup = GetFarHitsPreGroup(path, tryChargeCut);
      if(!kinkGroup->size()){
        kinkGroup = GetFarHitsGroup(path, tryChargeCut);
      };
    };

    TVector3 kinkPosXYZ = kinkGroup->GetAveragePosXYZ();

    std::vector<ND::TTPCPathVolume*>::iterator path1EndIt = path->end();
    bool broken = false;

    if(kinkGroup->size()){
      double minDist2=99999999.;

      // find closest point to kink to break around
      for(std::vector<ND::TTPCPathVolume*>::iterator pointIt = path->begin(); pointIt != path->end(); ++pointIt){
        ND::TTPCPathVolume* point = *pointIt;

        TVector3 diff = point->GetPosXYZ() - kinkPosXYZ;
        double diffDist2 = diff.Mag2();
        if(diffDist2 < minDist2){
          minDist2 = diffDist2;
          path1EndIt = pointIt;
          broken = true;
        };
      };
    };

    // build copy of old hits so the code doesn't get confused when separating them
    ND::THandle<ND::TTPCVolGroup> newExtHits1 (new ND::TTPCVolGroup(fLayout));
    ND::THandle<ND::TTPCVolGroup> newExtHits2 (new ND::TTPCVolGroup(fLayout));
    newExtHits1->AddHits(path->GetExtendedHits());
    newExtHits2->AddHits(path->GetExtendedHits());

    outPath1->AddExtendedHits(newExtHits1);
    outPath1->AddBackHits(path->GetBackHits());
    outPath1->SetBackIsVertex(path->GetBackIsVertex());

    // if broken, tie up this path and start another
    if(broken){
      outPath1->AddFrontHits(kinkGroup);
      outPath1->SetFrontIsVertex(true);

      outPath2->AddExtendedHits(newExtHits2);
      outPath2->AddBackHits(kinkGroup);
      outPath2->SetBackIsVertex(true);

      outPath2->AddFrontHits(path->GetFrontHits());
      outPath2->SetFrontIsVertex(path->GetFrontIsVertex());
    }
    else{
      outPath1->AddFrontHits(path->GetFrontHits());
      outPath1->SetFrontIsVertex(path->GetFrontIsVertex());
    };

    // add hits to first path up to its end
    for(std::vector<ND::TTPCPathVolume*>::iterator pointIt = path->begin(); pointIt != path1EndIt; ++pointIt){
      ND::TTPCPathVolume* point = *pointIt;
      outPath1->AddCell(point->GetUnitVolume());
    };

    // add hits to first/same path up to its end
    if(path1EndIt != path->end())
    for(std::vector<ND::TTPCPathVolume*>::iterator pointIt = path1EndIt+1; pointIt != path->end(); ++pointIt){
      ND::TTPCPathVolume* point = *pointIt;
      outPath2->AddCell(point->GetUnitVolume());
    };

    outPaths.push_back(outPath1);
    if(outPath2->size()){
      outPaths.push_back(outPath2);

      // paths broken, so re-do hit finding
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > brokenPaths = std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >();

      brokenPaths.push_back(outPath1);
      brokenPaths.push_back(outPath2);

      SeparateHits(brokenPaths);
    };
  };
  std::vector< ND::THandle<ND::TTPCVolGroup> > outVertices = GetJunctionsFromPaths(outPaths);

  return outPaths;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetFilteredPaths(std::vector< ND::THandle<ND::TTPCVolGroup> > inPaths){
  // paths to return
  std::vector< ND::THandle<ND::TTPCVolGroup> > mergedPaths = std::vector< ND::THandle<ND::TTPCVolGroup> >();

  // add old paths to merged paths
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator inPath = inPaths.begin(); inPath != inPaths.end(); ++inPath){
    ND::THandle<ND::TTPCVolGroup> newGroup(new ND::TTPCVolGroup(fLayout));
    newGroup->AddHits(*inPath);
    mergedPaths.push_back(newGroup);
  };

  if(mergedPaths.size() < 1) return mergedPaths;
  // remove overlapping paths from merged paths
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator mergedPath1 = mergedPaths.begin(); mergedPath1 != mergedPaths.end(); ++mergedPath1)
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator mergedPath2 = mergedPath1+1; mergedPath2 != mergedPaths.end(); ++mergedPath2){
    // loop over all pairs of hits to determine if they are identical
    int nOverlaps=0;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator pt1 = (*mergedPath1)->begin(); pt1 != (*mergedPath1)->end(); ++pt1)
    for(std::map<long, ND::TTPCUnitVolume*>::iterator pt2 = (*mergedPath2)->begin(); pt2 != (*mergedPath2)->end(); ++pt2){
      if(pt1->first == pt2->first) nOverlaps++;
    };
  };

  return mergedPaths;
}

bool ND::TTPCVolGroupMan::GetPathVolOverlap(ND::THandle<ND::TTPCOrderedVolGroup> path, ND::TTPCUnitVolume* vol, ND::TTPCConnection::Type type){
  if(!path->size()) return false;

  // don't check when too close to the first or last point in the path
  ND::TTPCPathVolume* firstPoint = *(path->begin());
  ND::TTPCPathVolume* lastPoint = *(path->end()-1);
  ND::TTPCUnitVolume* firstVol = firstPoint->GetUnitVolume();
  ND::TTPCUnitVolume* lastVol = lastPoint->GetUnitVolume();
  for(std::vector<ND::TTPCPathVolume*>::iterator pnt = path->begin(); pnt != path->end(); ++pnt){
    ND::TTPCUnitVolume* curVol = (*pnt)->GetUnitVolume();

    // make sure you're not too close to the first point
    if(HitTest(curVol, firstVol, type)) continue;
    // make sure you're not too close to the last point
    if(HitTest(curVol, lastVol, type)) continue;
    // then check to see how close the volume is
    if(HitTest(curVol, vol, type)) return true;
  };
  return false;
}
bool ND::TTPCVolGroupMan::GetGroupGroupOverlap(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, ND::TTPCConnection::Type type, bool simple, bool checkLean){
  int distX;
  int distY;
  int distZ;
  fLayout->GetTypeDistances(distX, distY, distZ, type);

  if(simple){
    // avoid complicated checks

    ND::TTPCCellInfo3D grp1Pos = group1->GetAveragePad();
    ND::TTPCCellInfo3D grp2Pos = group2->GetAveragePad();

    float nx = (float)( grp2Pos.x - grp1Pos.x )/(float)distX;
    float ny = (float)( grp2Pos.y - grp1Pos.y )/(float)distY;
    float nz = (float)( grp2Pos.z - grp1Pos.z )/(float)distZ;

    return ( (nx*nx) + (ny*ny) + (nz*nz) <= 1.);
  };

  if(checkLean){
    // don't check in directions of opposing lean
    if ((group1->GetXLean() * group2->GetXLean()) < 0) distX = 0;
    if ((group1->GetYLean() * group2->GetYLean()) < 0) distY = 0;
    if ((group1->GetZLean() * group2->GetZLean()) < 0) distZ = 0;
  };

  // check if groups are in range of each other at all
  if( (group1->GetXMin()-group2->GetXMax())>distX || (group2->GetXMin()-group1->GetXMax())>distX ) return false;
  if( (group1->GetYMin()-group2->GetYMax())>distY || (group2->GetYMin()-group1->GetYMax())>distY ) return false;
  if( (group1->GetZMin()-group2->GetZMax())>distZ || (group2->GetZMin()-group1->GetZMax())>distZ ) return false;

  // check if any individual pairs of hits are in range
  for(std::map<long, ND::TTPCUnitVolume*>::iterator hit1El = group1->begin(); hit1El != group1->end(); ++hit1El)
  for(std::map<long, ND::TTPCUnitVolume*>::iterator hit2El = group2->begin(); hit2El != group2->end(); ++hit2El){
    ND::TTPCUnitVolume* vol1 = hit1El->second;
    ND::TTPCUnitVolume* vol2 = hit2El->second;

    if(std::abs(vol1->GetX() - vol2->GetX()) > distX) continue;
    if(std::abs(vol1->GetY() - vol2->GetY()) > distY) continue;
    if(std::abs(vol1->GetZ() - vol2->GetZ()) > distZ) continue;
    return true;
  };

  return false;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::MergeGroups(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2){
  ND::THandle<ND::TTPCVolGroup> newGroup (new ND::TTPCVolGroup(fLayout, group1->GetID()));

  newGroup->AddHits(group1);
  newGroup->AddHits(group2);

  return newGroup;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::ExperimentalMergeGroups(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2){
  ND::THandle<ND::TTPCVolGroup> newGroup = MergeGroups(group1, group2);
  BulkGroup(newGroup);

  return newGroup;
}
void ND::TTPCVolGroupMan::BulkGroup(ND::THandle<ND::TTPCVolGroup> group){
  int xMin = group->GetXMin();
  int xMax = group->GetXMax();
  int yMin = group->GetYMin();
  int yMax = group->GetYMax();
  int zMin = group->GetZMin();
  int zMax = group->GetZMax();

  // TODO: either add directly or by lookup; whichever is faster
  for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = fPrimaryHits->begin(); hitEl != fPrimaryHits->end(); ++hitEl){
    ND::TTPCUnitVolume* vol = hitEl->second;

    if(vol->GetX() >= xMin && vol->GetX() <= xMax){
      if(vol->GetY() >= yMin && vol->GetY() <= yMax){
        if(vol->GetZ() >= zMin && vol->GetZ() <= zMax){
          group->AddHit(vol);
        };
      };
    };
  };
}
void ND::TTPCVolGroupMan::BulkGroups(std::vector< ND::THandle<ND::TTPCVolGroup> > groups){
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator groupIt = groups.begin(); groupIt != groups.end(); ++groupIt){
    BulkGroup(*groupIt);
  };
}

ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetUnorderedGroup(ND::THandle<ND::TTPCOrderedVolGroup> in){
  ND::THandle<ND::TTPCVolGroup> group(new ND::TTPCVolGroup(fLayout));

  for(std::vector<ND::TTPCPathVolume*>::iterator it = in->begin(); it != in->end(); ++it){
    group->AddHit((*it)->GetUnitVolume());
  };

  return group;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetExtendedUnorderedGroup(ND::THandle<ND::TTPCOrderedVolGroup> in){
  // check if group has extended hits and produce them if not
  if(!in->GetHasExtendedHits()){
    BuildGroupFriends(in);
  };

  // return extended hits
  return in->GetExtendedHits();
}
void ND::TTPCVolGroupMan::BuildGroupFriends(ND::THandle<ND::TTPCOrderedVolGroup> in, ND::TTPCConnection::Type type){
  // max distance from path
  int typeX;
  int typeY;
  int typeZ;
  fLayout->GetTypeDistances(typeX, typeY, typeZ, type);

  // extended group to create
  ND::THandle<ND::TTPCVolGroup> group(new ND::TTPCVolGroup(fLayout));

  // add edge hits
  if(in->HasFrontHits()) group->AddHits(in->frontHitsBegin(), in->frontHitsEnd());
  if(in->HasBackHits()) group->AddHits(in->backHitsBegin(), in->backHitsEnd());

  // two different modes for searching path hits - use whichever involves fewest iterations
  int volIterations = (4./3.)*3.14 * typeX*typeY*typeZ;
  if(volIterations > fPrimaryHits->size()){
    for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = fPrimaryHits->begin(); hitEl != fPrimaryHits->end(); ++hitEl){
      bool nearEnough = false;
      for(std::vector<ND::TTPCPathVolume*>::iterator inputIt = in->begin(); inputIt != in->end(); inputIt++){
        ND::TTPCUnitVolume* inputVol = (*inputIt)->GetUnitVolume(); 

        int dx = std::abs( inputVol->GetX() - hitEl->second->GetX() );
        int dy = std::abs( inputVol->GetY() - hitEl->second->GetY() );
        int dz = std::abs( inputVol->GetZ() - hitEl->second->GetZ() );

        if(dx > typeX) continue;
        if(dy > typeY) continue;
        if(dz > typeZ) continue;

        float nx = (float)dx/(float)typeX;
        float ny = (float)dy/(float)typeY;
        float nz = (float)dz/(float)typeZ;
        if(( nx*nx + ny*ny + nz*nz) > 1. ) continue;

        nearEnough = true;
        break;
      };
      if(nearEnough) group->AddHit(hitEl->second);
    };
  }
  else{
    for(std::vector<ND::TTPCPathVolume*>::iterator inputIt = in->begin(); inputIt != in->end(); inputIt++){
      ND::TTPCUnitVolume* inputVol = (*inputIt)->GetUnitVolume(); 

      // hits near cell
      ND::THandle<ND::TTPCVolGroup> nearHits = GetNearHits(fPrimaryHits, inputVol->GetID(), type, true);
      // add all of them
      group->AddHits(nearHits);
    };
  };

  // add associated hits to the group
  in->AddExtendedHits(group);
}
void ND::TTPCVolGroupMan::ClusterGroupFriends(ND::THandle<ND::TTPCOrderedVolGroup> in, bool doClustering, bool checkX, bool partial){
  // note - path may be empty after this
  // check if group has extended hits and produce them if not
  if(!in->GetHasExtendedHits()){
    BuildGroupFriends(in);
  };
  // stop here if 'do clustering' is off
  if(!doClustering) return;
  // if needed, check whether path is an x-path before going further
  if(checkX){
    if(IsXPathCandidate(in)){
      in->SetIsXPath(true);
    };
  };

  // connect all assocaited hits
  in->DoClustering(partial);
}
std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >  ND::TTPCVolGroupMan::ClearEmptyPaths(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& paths){
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup > >::iterator pathDel = paths.begin();
  while(pathDel != paths.end()){
    if(!(*pathDel)->size()) paths.erase(pathDel);
    else pathDel ++;
  };
  return paths;
}
void ND::TTPCVolGroupMan::BuildAllFriends(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // build friends for paths
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt) BuildGroupFriends(*pathIt);

  // ensure that none share hits
  SeparateHits(paths);
}
void ND::TTPCVolGroupMan::SeparateHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // loop over all pairs of groups
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It)
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
    ND::THandle<ND::TTPCOrderedVolGroup> path1 = *path1It;
    ND::THandle<ND::TTPCOrderedVolGroup> path2 = *path2It;

    // take associated hits
    ND::THandle<ND::TTPCVolGroup> extendedHits1 = path1->GetExtendedHits();
    ND::THandle<ND::TTPCVolGroup> extendedHits2 = path2->GetExtendedHits();

    // loop over path associated hits
    for(std::map<long, ND::TTPCUnitVolume*>::iterator extHit1El = extendedHits1->begin(); extHit1El != extendedHits1->end(); ++extHit1El)
    for(std::map<long, ND::TTPCUnitVolume*>::iterator extHit2El = extendedHits2->begin(); extHit2El != extendedHits2->end(); ++extHit2El){
      // shared hit - break tie in favour of closest path
      if(extHit1El->first == extHit2El->first){
        // ignore if already marked for clearing
        if(!extHit1El->second) continue;

        double path1Dist = GetMinDistance(path1, extHit1El->second);
        double path2Dist = GetMinDistance(path2, extHit1El->second);

        if(path1Dist > path2Dist){
          // mark hit for removal from from path1
          extendedHits1->MarkHit(extHit1El->first);
        }
        else{
          // mark hit for removal from path2
          extendedHits2->MarkHit(extHit1El->first);
        };
      };
    };
    // clear hits marked for deletion
    extendedHits1->ClearMarked();
    extendedHits2->ClearMarked();

    // renew path associated hits
    path1->AddExtendedHits(extendedHits1);
    path2->AddExtendedHits(extendedHits2);
  };
}
void ND::TTPCVolGroupMan::SeparateXPathHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& paths){
  // loop over all paths
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

    if(path->GetIsXPath()){
      // get all hits in path
      ND::THandle<ND::TTPCVolGroup> pathHits (new ND::TTPCVolGroup(fLayout));
      for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
        ND::TTPCPathVolume* pathVol = *pathVolIt;
        for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
          ND::TTPCUnitVolume* vol = *volIt;
          pathHits->AddHit(vol);
        };
      };

      // add those hits to junctions either side; later code will take care of any overlap
      if(path->GetFrontIsVertex()){
        path->GetFrontHits()->AddHits(pathHits);
      };
      if(path->GetBackIsVertex()){
        path->GetBackHits()->AddHits(pathHits);
      };
    };
  };

  // now remove all x-paths
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator reaper = paths.begin();
  while(reaper != paths.end()){
    if((*reaper)->GetIsXPath()){
      paths.erase(reaper);
    }
    else{
      ++reaper;
    };
  };
}

void ND::TTPCVolGroupMan::SeparateEmptyClusters(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;
    path->Clean();
  };
}
void ND::TTPCVolGroupMan::SeparateClusterHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  int mergeX;
  int mergeY;
  int mergeZ;
  fLayout->GetTypeDistances(mergeX, mergeY, mergeZ, ND::TTPCConnection::clusterMerge);

  // get rid of any clusters whos seed is shared with other paths
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It)
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
    ND::THandle<ND::TTPCOrderedVolGroup> path1 = *path1It;
    ND::THandle<ND::TTPCOrderedVolGroup> path2 = *path2It;

    for(std::vector<ND::TTPCPathVolume*>::iterator cluster1It = path1->begin(); cluster1It != path1->end(); ++cluster1It)
    for(std::vector<ND::TTPCPathVolume*>::iterator cluster2It = path2->begin(); cluster2It != path2->end(); ++cluster2It){
      ND::TTPCPathVolume* cluster1 = *cluster1It;
      ND::TTPCPathVolume* cluster2 = *cluster2It;

      if(!cluster1) continue;
      if(!cluster2) continue;

      ND::TTPCUnitVolume* vol1 = cluster1->GetUnitVolume();
      ND::TTPCUnitVolume* vol2 = cluster2->GetUnitVolume();

      if(vol1 == vol2){
        cluster1->ClearFriends();
        cluster2->ClearFriends();
      };
    };
  };

  // loop over pairs of paths, checking for overlapping clusters (ignore x-paths for now)
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It){
    ND::THandle<ND::TTPCOrderedVolGroup> path1 = *path1It;
    if(path1->GetIsXPath()) continue;

    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
      ND::THandle<ND::TTPCOrderedVolGroup> path2 = *path2It;
      if(path2->GetIsXPath()) continue;

      for(std::vector<ND::TTPCPathVolume*>::iterator cluster1It = path1->begin(); cluster1It != path1->end(); ++cluster1It)
      for(std::vector<ND::TTPCPathVolume*>::iterator cluster2It = path2->begin(); cluster2It != path2->end(); ++cluster2It){
        ND::TTPCPathVolume* cluster1 = *cluster1It;
        ND::TTPCPathVolume* cluster2 = *cluster2It;
        if(!cluster1) continue;
        if(!cluster2) continue;

        // ensure clusters are the same orientation
        if(cluster1->GetIsVertical() && !cluster2->GetIsVertical()) continue;
        if(!cluster1->GetIsVertical() && cluster2->GetIsVertical()) continue;

        bool verticalMerge = cluster1->GetIsVertical();
        // check if merge is possible
        if((cluster1->GetXMin() - cluster2->GetXMax()) > mergeX) continue;
        if((cluster2->GetXMin() - cluster1->GetXMax()) > mergeX) continue;
        if(verticalMerge){
          if(cluster1->GetZ() != cluster2->GetZ()) continue;
          if((cluster1->GetYMin() - cluster2->GetYMax()) > mergeY) continue;
          if((cluster2->GetYMin() - cluster1->GetYMax()) > mergeY) continue;
        }
        else{
          if(cluster1->GetY() != cluster2->GetY()) continue;
          if((cluster1->GetZMin() - cluster2->GetZMax()) > mergeZ) continue;
          if((cluster2->GetZMin() - cluster1->GetZMax()) > mergeZ) continue;
        };

        // if clusters are close, check for a gap in one of them to re-break around
        std::vector<int> positions;
        for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = cluster1->GetFriendsBegin(); volIt != cluster2->GetFriendsEnd(); ++volIt){
          if(volIt == cluster1->GetFriendsEnd()) volIt = cluster2->GetFriendsBegin();
          ND::TTPCUnitVolume* vol = *volIt;

          if(verticalMerge) positions.push_back(vol->GetY());
          else positions.push_back(vol->GetZ());
        };

        // find biggest gap
        if(!positions.size()) continue;
        std::sort(positions.begin(), positions.end());
        int gapPos = -1;
        int gapSize = 0;
        int prevPos = *positions.begin();
        for(std::vector<int>::iterator positionIt = positions.begin(); positionIt != positions.end(); ++positionIt){
          int position = *positionIt;
          int gap = position-prevPos;
          prevPos = position;

          if(gap > gapSize){
            gapSize = gap;
            gapPos = position-1;
          };
        };

        // if a suitable gap is found, break into two new clusters around it
        if(( verticalMerge && gapSize>=mergeZ ) || ( !verticalMerge && gapSize>=mergeY )){
          bool normalHierarchy = (verticalMerge && (cluster1->GetY()>cluster2->GetY())) || (!verticalMerge && (cluster1->GetZ()>cluster2->GetZ()));

          // add hits above or below the gap to the appropriate cluster
          for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = cluster1->GetFriendsBegin(); volIt != cluster1->GetFriendsEnd(); ++volIt){
            ND::TTPCUnitVolume* vol = *volIt;

            if( normalHierarchy^((verticalMerge && (vol->GetY()>gapPos)) || (!verticalMerge && (vol->GetZ()>gapPos))) ){
              cluster2->AddFriend(vol);
              cluster1->MarkFriend(volIt);
            };
          };
          cluster1->ClearMarked();
          for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = cluster2->GetFriendsBegin(); volIt != cluster2->GetFriendsEnd(); ++volIt){
            ND::TTPCUnitVolume* vol = *volIt;

            if( !normalHierarchy^((verticalMerge && (vol->GetY()>gapPos)) || (!verticalMerge && (vol->GetZ()>gapPos))) ){
              cluster1->AddFriend(vol);
              cluster2->MarkFriend(volIt);
            };
          };
          cluster2->ClearMarked();
        }
        // otherwise just delete the initial clusters and let the nearest vertex eat their hits
        else{
          if(*cluster1It) delete *cluster1It;
          *cluster1It = 0;
          if(*cluster2It) delete *cluster2It;
          *cluster2It = 0;
        };
      };

      // clear any marked clusters
      path1->Clean();
      path2->Clean();
    };
  };
};

void ND::TTPCVolGroupMan::SeparateAnomHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // ready to add hits to junctions
  std::vector< ND::THandle<ND::TTPCVolGroup> > vertices = GetJunctionsFromPaths(paths);

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

    if(path->GetFrontIsVertex() || path->GetBackIsVertex()){
      if(path->GetFrontIsVertex()){
        unsigned int junctionID = path->GetFrontID();
        std::vector<ND::TTPCUnitVolume*> hitsToAdd = SeparateAnomHitsPath(path, -1);

        // now add discarded hits to junction
        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
          ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;
          if(junctionID == vertex->GetID()){
            vertex->AddHits(hitsToAdd);
          };
        };
      };
      if(path->GetBackIsVertex()){
        unsigned int junctionID = path->GetBackID();
        std::vector<ND::TTPCUnitVolume*> hitsToAdd = SeparateAnomHitsPath(path, 1);

        // now add discarded hits to junction
        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
          ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;
          if(junctionID == vertex->GetID()){
            vertex->AddHits(hitsToAdd);
          };
        };
      };
    };
  };
}
std::vector<ND::TTPCUnitVolume*> ND::TTPCVolGroupMan::SeparateAnomHitsPath(ND::THandle<ND::TTPCOrderedVolGroup> path, int checkDir){
  int checkDist = fLayout->GetAnomCheckDist();
  int projectDist = fLayout->GetAnomProjectDist();
  double maxOffs = fLayout->GetAnomMaxOffs();

  // be prepared to add hits to junction
  std::vector<ND::TTPCUnitVolume*> hitsToAdd;

  // return empty handed if path is too small
  if(path->size() <= (checkDist + projectDist)){
    return hitsToAdd;
  };

  int rangeZero;
  int rangeStart;
  int rangeEnd;
  if(checkDir > 0){
    rangeZero = 0;
    rangeStart = rangeZero + checkDist;
    rangeEnd = rangeStart + projectDist;
  }
  else if(checkDir < 0){
    rangeZero = path->size()-1;
    rangeStart = rangeZero - checkDist;
    rangeEnd = rangeStart - projectDist;
  }
  else{
    return hitsToAdd;
  };

  // look for hits with max time deviation from expectation, approximated by extrapolating a straight line between the start and end of test range
  ND::TTPCPathVolume* range1 = path->at(rangeStart);
  ND::TTPCPathVolume* range2 = path->at(rangeEnd);

  TVector3 startPos = GetAvgPosRep(range1);
  TVector3 endPos = GetAvgPosRep(range2);
  TVector3 parNorm = (endPos - startPos).Unit();
  bool xOnly = false;
  double maxDev = -1.;

  // find max in the control region
  for(int i=rangeStart; i != rangeEnd; i+= checkDir){
    ND::TTPCPathVolume* pathVol = path->at(i);
    for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
      ND::TTPCUnitVolume* vol = *volIt;
      for(int sign=-1; sign <= 1; sign += 2){
        TVector3 volPos = GetAvgPosRep(vol, sign);
        TVector3 diff = startPos - volPos;
        TVector3 parallel = parNorm*diff.Dot(parNorm);
        TVector3 dist = diff - parallel;

        if(xOnly){
          maxDev = std::max(maxDev, std::abs(dist.X()));
        }
        else{
          maxDev = std::max(maxDev, dist.Mag());
        };
      };
    };
  };
  // now remove anything too far above it in the check region
  maxDev += maxOffs;
  for(int i=rangeZero; i != rangeStart; i+= checkDir){
    ND::TTPCPathVolume* pathVol = path->at(i);
    bool removed = false;
    for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
      ND::TTPCUnitVolume* vol = *volIt;
      for(int sign=-1; sign <= 1; sign += 2){
        TVector3 volPos = GetAvgPosRep(vol, sign);
        TVector3 diff = startPos - volPos;
        TVector3 parallel = parNorm*diff.Dot(parNorm);
        TVector3 dist = diff - parallel;

        if(xOnly){
          if(std::abs(diff.X()) > maxDev){
            // friend is too far away; remove it and add to junction
            hitsToAdd.push_back(vol);
            pathVol->MarkFriend(volIt);
            removed = true;
            break;
          };
        }
        else{
          if(dist.Mag() > maxDev){
            // friend is too far away; remove it and add to junction
            hitsToAdd.push_back(vol);
            pathVol->MarkFriend(volIt);
            removed = true;
            break;
          };
        };
      };
    };
    if(removed){
      pathVol->ClearMarked();
    };
  };
  path->Clean();

  return hitsToAdd;
}

void ND::TTPCVolGroupMan::SeparateJunctionHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  int mergeX;
  int mergeY;
  int mergeZ;
  fLayout->GetTypeDistances(mergeX, mergeY, mergeZ, ND::TTPCConnection::clusterMerge);

  // first build list of potential vertices (path groups tagged as vertices)
  std::vector< ND::THandle<ND::TTPCVolGroup> > tempVertices = GetJunctionsFromPaths(paths);

  // select un-duplicated vertices and merge together any that overlap
  std::vector< ND::THandle<ND::TTPCVolGroup> > vertices;
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator tempVertexIt = tempVertices.begin(); tempVertexIt != tempVertices.end(); tempVertexIt++){
    ND::THandle<ND::TTPCVolGroup> tempVertex = *tempVertexIt;

    // check if each vertex overlaps with a current one and merge them if it does
    bool vertexFound = false;
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
      ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;

      // as a timesaver, ensure vertices might overlap before doing more intensive checks
      bool potentialOverlap = true;
      if(vertex->GetXMin() > tempVertex->GetXMax() || tempVertex->GetXMin() > vertex->GetXMax()) potentialOverlap = false;
      if(vertex->GetYMin() > tempVertex->GetYMax() || tempVertex->GetYMin() > vertex->GetYMax()) potentialOverlap = false;
      if(vertex->GetZMin() > tempVertex->GetZMax() || tempVertex->GetZMin() > vertex->GetZMax()) potentialOverlap = false;

      bool thisVertexFound = false;
      if(potentialOverlap){
        // overlaps between the two
        for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = tempVertex->begin(); volEl != tempVertex->end(); ++volEl){
          if(vertex->Contains(volEl->first)){
            vertexFound = thisVertexFound = true;
            break;
          };
        };
      };
      if(thisVertexFound){
        // merge
        vertex->AddHits(tempVertex);
      };
    };
    if(!vertexFound){
      // add vertex if it doesn't already exist
      vertices.push_back(tempVertex);
    };
  };

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;
    // check if any HV cluster overlaps with any vertex and swap those that do from the cluster to the vertex
    std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin();
    while(pathVolIt != path->end()){
      if(!*pathVolIt) path->erase(pathVolIt);
      else pathVolIt++;
    };
  };

  // now merge overlapping HV clusters from paths into vertices
  bool mightMerge = true;
  while(mightMerge){
    mightMerge = false;
    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

      // build sets of hits in horizontal clusters and vertical clusters
      std::set<ND::TTPCUnitVolume*> clusteredH;
      std::set<ND::TTPCUnitVolume*> clusteredV;
      if(!path->GetIsXPath()){
        // add hits in horizontal and vertical clusters to their respective sets
        for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
          ND::TTPCPathVolume* pathVol = *pathVolIt;
          if(!pathVol) continue;

          for(std::vector<ND::TTPCUnitVolume*>::iterator pathFriendIt = pathVol->GetFriendsBegin(); pathFriendIt != pathVol->GetFriendsEnd(); ++pathFriendIt){
            ND::TTPCUnitVolume* pathFriend = *pathFriendIt;

            if(pathVol->GetIsVertical()){
              clusteredV.insert(pathFriend);
            }
            else{
              clusteredH.insert(pathFriend);
            };
          };
        };
      };

      for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
        ND::TTPCPathVolume* pathVol = *pathVolIt;
        if(!pathVol) continue;

        bool verticalMerge = pathVol->GetIsVertical();

        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
          ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;

          int xMin = vertex->GetXMin();
          int xMax = vertex->GetXMax();
          int yMin = vertex->GetYMin();
          int yMax = vertex->GetYMax();
          int zMin = vertex->GetZMin();
          int zMax = vertex->GetZMax();

          xMin -= mergeX;
          xMax += mergeX;
          if(verticalMerge){
            yMin -= mergeY;
            yMax += mergeY;
          }
          else{
            zMin -= mergeZ;
            zMax += mergeZ;
          };

          // as a timesaver, ensure cluster is actually within range of vertex before doing more intensive checks
          if(!pathVol->GetIsXCluster()){
            if(pathVol->GetXMax() < xMin || pathVol->GetXMin() > xMax) continue;
            if(verticalMerge){
              if(pathVol->GetYMax() < yMin || pathVol->GetYMin() > yMax) continue;
              if(pathVol->GetZ() < zMin || pathVol->GetZ() > zMax) continue;
            }
            else{
              if(pathVol->GetY() < yMin || pathVol->GetY() > yMax) continue;
              if(pathVol->GetZMax() < zMin || pathVol->GetZMin() > zMax) continue;
            };
          };

          // loop over all hits in cluster and see if any are contained in the vertex
          bool clusterOverlaps = false;
          for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); volIt++){
            ND::TTPCUnitVolume* vol = *volIt;
            if(!vol) continue;

            // first check if vol might actually overlap to save time
            if(vol->GetX() < xMin || vol->GetX() > xMax) continue;
            if(vol->GetY() < yMin || vol->GetY() > yMax) continue;
            if(vol->GetZ() < zMin || vol->GetZ() > zMax) continue;

            // potential hit - do a more intense search
            for(std::map<long, ND::TTPCUnitVolume*>::iterator vertexEl = vertex->begin(); vertexEl != vertex->end(); ++vertexEl){
              ND::TTPCUnitVolume* vertexVol = vertexEl->second;

              int dX = std::abs(vol->GetX() - vertexVol->GetX());
              int dY = std::abs(vol->GetY() - vertexVol->GetY());
              int dZ = std::abs(vol->GetZ() - vertexVol->GetZ());

              if(dX > mergeX) continue;
              if(dY > (verticalMerge ? mergeY : 0)) continue;
              if(dZ > (verticalMerge ? 0 : mergeZ)) continue;

              // if all conditions are passed, an overlap is found
              clusterOverlaps = true;

              // add to the vertex if not spotted in another cluster
              if(verticalMerge){
                if(!clusteredH.count(vol)){
                  vertex->AddHit(vol);
                };
                clusteredV.erase(vol);
              }
              else{
                if(!clusteredV.count(vol)){
                  vertex->AddHit(vol);
                };
                clusteredH.erase(vol);
              };

              // remove from its group
              pathVol->MarkFriend(volIt);

              // exit the vertex loop
              break;
            };
          };

          // if an overlap was found, clean up tagged entries in path volume and continue merge at next iteration
          if(clusterOverlaps){
            mightMerge = true;
            pathVol->ClearMarked();

            if(!pathVol->GetHasCluster()){
              *pathVolIt = 0;
            };

            break;
          };
        };
      };
    };
  };

  // also merge in any isolated clusters
  int isoDist = fLayout->GetHVClusterMaxIso();
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

    if(path->GetIsXPath()) continue;
    if(path->size() < 2) continue;

    // look for gap at the back
    if(path->GetBackIsVertex()){
      bool possibleIso = true;
      while(possibleIso){
        possibleIso = false;
        ND::TTPCPathVolume* beginVol = 0;
        std::vector<ND::TTPCPathVolume*>::iterator beginVolIt;
        for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
          if(*pathVolIt){
            beginVolIt = pathVolIt;
            beginVol = *pathVolIt;
            break;
          };
        };

        if(beginVol){
          bool wasVertical = beginVol->GetIsVertical();
          int lastPos = wasVertical ? beginVol->GetZ() : beginVol->GetY();

          int i=0;
          for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = beginVolIt+1; pathVolIt != path->end(); ++pathVolIt){
            ND::TTPCPathVolume* pathVol = *pathVolIt;
            if(!pathVol) continue;

            if(wasVertical==pathVol->GetIsVertical()){
              int pos = pathVol->GetIsVertical() ? pathVol->GetZ() : pathVol->GetY();

              // check if there's a gap and add everything before it to the vertex if so
              if(std::abs(pos - lastPos) > 1){
                for(std::vector<ND::TTPCPathVolume*>::iterator gapVolIt = path->begin(); gapVolIt != pathVolIt; ++gapVolIt){
                  if(*gapVolIt){
                    // no need to add - vertex should sweep up these unused hits automatically
                    //vertex->AddHits(gapVol->GetFriendsBegin(), gapVol->GetFriendsEnd());
                    *gapVolIt = 0;
                    possibleIso = true;
                  };
                };
                if(possibleIso) break;
              };
            };

            wasVertical = pathVol->GetIsVertical();
            lastPos = wasVertical ? pathVol->GetZ() : pathVol->GetY();

            i++;
            if(i > isoDist) break;
          };
        };
      };
    };
    // look for gap at the front
    if(path->GetFrontIsVertex()){
      bool possibleIso = true;
      while(possibleIso){
        possibleIso = false;
        std::vector<ND::TTPCPathVolume*>::reverse_iterator endVolRIt;
        ND::TTPCPathVolume* endVol = 0;
        for(std::vector<ND::TTPCPathVolume*>::reverse_iterator pathVolRIt = path->rbegin(); pathVolRIt != path->rend(); ++pathVolRIt){
          if(*pathVolRIt){
            endVolRIt = pathVolRIt;
            endVol = *pathVolRIt;
            break;
          };
        };

        if(endVol){
          bool wasVertical = endVol->GetIsVertical();
          int lastPos = wasVertical ? endVol->GetZ() : endVol->GetY();

          int i=0;
          for(std::vector<ND::TTPCPathVolume*>::reverse_iterator pathVolRIt = endVolRIt; pathVolRIt != path->rend(); ++pathVolRIt){
            ND::TTPCPathVolume* pathVol = *pathVolRIt;
            if(!pathVol) continue;

            if(wasVertical==pathVol->GetIsVertical()){
              int pos = pathVol->GetIsVertical() ? pathVol->GetZ() : pathVol->GetY();

              // check if there's a gap and add everything before it to the vertex if so
              if(std::abs(pos - lastPos) > 1){
                for(std::vector<ND::TTPCPathVolume*>::reverse_iterator gapVolRIt = path->rbegin(); gapVolRIt != pathVolRIt; ++gapVolRIt){
                  if(*gapVolRIt){
                    // no need to add - vertex should sweep up these unused hits automatically
                    //vertex->AddHits(gapVol->GetFriendsBegin(), gapVol->GetFriendsEnd());
                    *gapVolRIt = 0;
                    possibleIso = true;
                  };
                };
                if(possibleIso) break;
              };
            };

            wasVertical = pathVol->GetIsVertical();
            lastPos = wasVertical ? pathVol->GetZ() : pathVol->GetY();

            i++;
            if(i > isoDist) break;
          };
        };
      };
    };
  };

  // clear superfluous HV clusters from paths into vertices
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;
    // check if any HV cluster overlaps with any vertex and swap those that do from the cluster to the vertex
    std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin();
    while(pathVolIt != path->end()){
      if(!*pathVolIt) path->erase(pathVolIt);
      else pathVolIt++;
    };
  };

  // copy expanded vertices to their initial positions
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
    ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;

    // replace any overlapping path beginning or end with this
    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

      // overlap with front hits
      bool frontOverlaps = path->GetFrontID() == vertex->GetID();
      bool backOverlaps = path->GetBackID() == vertex->GetID();

      // set as new vector
      if(frontOverlaps) path->AddFrontHits(vertex);
      if(backOverlaps) path->AddBackHits(vertex);
    };
  };
}
void ND::TTPCVolGroupMan::EnforceOrdering(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // loop over all paths
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt !=paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;
    // count turning points
    int turningPoints90=0;
    bool wasVertical=true;
    for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
      ND::TTPCPathVolume* pathVol = *pathVolIt;

      bool isVertical = pathVol->GetIsVertical();
      if(isVertical ^ wasVertical){
        if(pathVolIt != path->begin()) turningPoints90++;
      };
      wasVertical = isVertical;
    };
    // define turning points as half changes from horizontal to vertical or vice versa
    int turningPoints = turningPoints90/2;

    // case #1: positive z
    if(turningPoints<1){
      path->OrderForwardsDirection();
    }
    // case #2: negative curvature
    else if(turningPoints<2){
      path->OrderNegativeCurvature();
    }
    // case #3:
    else{
      //TODO: implement as starting from junction?
      path->OrderForwardsDirection();
    }
  };
}
bool ND::TTPCVolGroupMan::GetOverlaps(ND::THandle<ND::TTPCVolGroup> group1, ND::THandle<ND::TTPCVolGroup> group2, bool checkHits){
  // basic check
  if(group1->GetXMin() > group2->GetXMax() || group2->GetXMin() > group1->GetXMax()) return false;
  if(group1->GetYMin() > group2->GetYMax() || group2->GetYMin() > group1->GetYMax()) return false;
  if(group1->GetZMin() > group2->GetZMax() || group2->GetZMin() > group1->GetZMax()) return false;

  // if not flagged to check hits, just return here
  if(!checkHits) return true;

  // compare for individual hits for overlaps
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el1 = group1->begin(); el1 != group1->end(); ++el1)
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el2 = group2->begin(); el2 != group2->end(); ++el2)
  if(el1->first == el2->first) return true;

  return false;
}

void ND::TTPCVolGroupMan::ResetVertexStatuses(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths, bool partial){
  // first build list of potential vertices (path groups tagged as vertices)
  std::vector< ND::THandle<ND::TTPCVolGroup> > vertices = GetJunctionsFromPaths(paths);

  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
    ND::THandle<ND::TTPCVolGroup> vertex = *vertexIt;
    int connectedPaths = 0;

    // count connected paths
    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

      // overlap with front hits
      bool frontOverlaps = path->GetFrontID() == vertex->GetID();
      bool backOverlaps = path->GetBackID() == vertex->GetID();

      // increment connected paths
      if(frontOverlaps) connectedPaths++;
      if(backOverlaps) connectedPaths++;
    };

    // status as vertex
    bool isVertex = connectedPaths>1;

    // replace any overlapping path beginning or end with this
    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

      // overlap with front hits
      bool frontOverlaps = path->GetFrontID() == vertex->GetID();
      bool backOverlaps = path->GetBackID() == vertex->GetID();


      // set new status as vertex, adding hits if the status has changed
      bool redoClustering = false;
      if(frontOverlaps){
        if(path->GetFrontIsVertex() != isVertex){
          redoClustering = true;
        };
        path->SetFrontIsVertex(isVertex);
      }
      else if(backOverlaps){
        if(path->GetBackIsVertex() != isVertex){
          redoClustering = true;
        };
        path->SetBackIsVertex(isVertex);
      };

      // do NOT allow hits the same vertex to be duplicated
      if(path->GetFrontID() == path->GetBackID()){
        if(path->GetFrontIsVertex() && path->GetBackIsVertex()){
          path->SetBackIsVertex(false);
        };
      }
      else if(redoClustering){
        path->GetExtendedHits()->AddHits(vertex);
        path->DoClustering(partial);
      };
    };
  };
}

ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetGapFillHits(ND::THandle<ND::TTPCVolGroup> in){
  ND::THandle<ND::TTPCVolGroup> gapFillHits(new ND::TTPCVolGroup(fLayout));

  // identify hits at upper and lower edges of MM volumes in x, y and z
  ND::THandle<ND::TTPCVolGroup> hitsUX(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> hitsUY(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> hitsUZ(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> hitsDX(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> hitsDY(new ND::TTPCVolGroup(fLayout));
  ND::THandle<ND::TTPCVolGroup> hitsDZ(new ND::TTPCVolGroup(fLayout));

  // fill edge groups with found hits
  for(std::map<long, ND::TTPCUnitVolume*>::iterator it = in->begin(); it != in->end(); ++it){
    if (it->second->GetEdgeX() == 1) hitsUX->AddHit(it->second);
    if (it->second->GetEdgeY() == 1) hitsUY->AddHit(it->second);
    if (it->second->GetEdgeZ() == 1) hitsUZ->AddHit(it->second);
    if (it->second->GetEdgeX() == -1) hitsDX->AddHit(it->second);
    if (it->second->GetEdgeY() == -1) hitsDY->AddHit(it->second);
    if (it->second->GetEdgeZ() == -1) hitsDZ->AddHit(it->second);
  };

  // get extrapolations over the gap from those groups
  ND::THandle<ND::TTPCVolGroup> projectionsUX = GetGapProjections(hitsUX, 1, 0, 0, in->GetHitMap());
  ND::THandle<ND::TTPCVolGroup> projectionsUY = GetGapProjections(hitsUY, 0, 1, 0, in->GetHitMap());
  ND::THandle<ND::TTPCVolGroup> projectionsUZ = GetGapProjections(hitsUZ, 0, 0, 1, in->GetHitMap());
  ND::THandle<ND::TTPCVolGroup> projectionsDX = GetGapProjections(hitsDX, -1, 0, 0, in->GetHitMap());
  ND::THandle<ND::TTPCVolGroup> projectionsDY = GetGapProjections(hitsDY, 0, -1, 0, in->GetHitMap());
  ND::THandle<ND::TTPCVolGroup> projectionsDZ = GetGapProjections(hitsDZ, 0, 0, -1, in->GetHitMap());

  // add extrapolations to the main group as pseudo-hits
  gapFillHits->MergeHits(projectionsUX);
  gapFillHits->MergeHits(projectionsUY);
  gapFillHits->MergeHits(projectionsUZ);
  gapFillHits->MergeHits(projectionsDX);
  gapFillHits->MergeHits(projectionsDY);
  gapFillHits->MergeHits(projectionsDZ);

  return gapFillHits;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetGapProjections(ND::THandle<ND::TTPCVolGroup> in, int dirX, int dirY, int dirZ, std::map<long, ND::TTPCUnitVolume*> refHits){
  ND::THandle<ND::TTPCVolGroup> projGroup(new ND::TTPCVolGroup(fLayout));

  if(in->size() < 1) return projGroup;
  // break this group into connected sub-groups
  std::vector< ND::THandle<ND::TTPCVolGroup> > groups;
  for(int i=0; i<9999; i++){
    if(in->size() < 1) break;

    std::map<long, ND::TTPCUnitVolume*>::iterator el = in->begin();

    ND::THandle<ND::TTPCVolGroup> newGroup(new ND::TTPCVolGroup(fLayout));
    // build connected hits from hits in this group recursively
    RecursiveFriendBuild(el->first, newGroup, in, ND::TTPCConnection::path);

    // add to list of unique groups
    groups.push_back(newGroup);
  };

  // process each unique group individually
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator grp = groups.begin(); grp != groups.end(); ++grp){
    // work out parameters for extrapolation
    float q = (*grp)->GetCharge() / 2;
    q = 1000.;
    ND::TTPCCellInfo3D cell = (*grp)->GetAveragePad();
    int x = cell.x; 
    int y = cell.y;
    int z = cell.z;
    float sigmaX = (*grp)->GetSigmaPadX();
    float sigmaY = (*grp)->GetSigmaPadY();
    float sigmaZ = (*grp)->GetSigmaPadZ();

    ND::THandle<ND::TTPCVolGroup> refGroup(new TTPCVolGroup(fLayout));
    refGroup->AddHits(refHits);

    ND::THandle<ND::TTPCVolGroup> nearGroup;

    // find nearest group of cells in given direction to work out direction to project over
    if (dirX == 1 || dirX == -1) nearGroup = GetNearestGroup(refGroup, x + dirX*(fLayout->GetGapOffsetX() + 1), y, z, 1);
    else if (dirY == 1 || dirY == -1) nearGroup = GetNearestGroup(refGroup, x, y + dirY*(fLayout->GetGapOffsetY() + 1), z, 2);
    else if (dirZ == 1 || dirZ == -1) nearGroup = GetNearestGroup(refGroup, x, y, z + dirZ*(fLayout->GetGapOffsetZ() + 1), 3);
    else nearGroup = ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));

    if(nearGroup->size() < 1) continue;

    // work out direction to project into as difference between this group of cells and the nearest group
    TVector3 offset = nearGroup->GetAveragePosition() - (*grp)->GetAveragePosition();
    float diffX = offset.X() / fLayout->GetXCellSize();
    float diffY = offset.Y() / fLayout->GetPadPitchY();
    float diffZ = offset.Z() / fLayout->GetPadPitchZ();

    // extend across gap as a series of 2D gaussian distributions of charge in either zy, zx or yx planes
    if (dirX==-1 || dirX==1){
      float dY_dX = diffY / diffX;
      float dZ_dX = diffZ / diffX;
      for(int i=1; i<=fLayout->GetGapOffsetX(); i++){
        int i2 = dirX*i;
        int x2 = x + i2;
        int y2 = y + int(dY_dX*i2);
        int z2 = z + int(dZ_dX*i2);
        projGroup->AddNewPseudoGauss(x2, y2, z2, q, 1, sigmaX, sigmaY, sigmaZ);
      };
    };
    if (dirY==-1 || dirY==1){
      float dZ_dY = diffZ / diffY;
      float dX_dY = diffX / diffY;
      for(int j=1; j<=fLayout->GetGapOffsetY(); j++){
        int j2 = dirY*j;
        int x2 = x + int(dX_dY*j2);
        int y2 = y + j2;
        int z2 = z + int(dZ_dY*j2);
        projGroup->AddNewPseudoGauss(x2, y2, z2, q, 2, sigmaX, sigmaY, sigmaZ);
      };
    };
    if (dirZ==-1 || dirZ==1){
      float dX_dZ = diffX / diffZ;
      float dY_dZ = diffY / diffZ;
      for(int k=1; k<=fLayout->GetGapOffsetZ(); k++){
        int k2 = dirZ*k;
        int x2 = x + int(dX_dZ*k2);
        int y2 = y + int(dY_dZ*k2);
        int z2 = z + k2;
        projGroup->AddNewPseudoGauss(x2, y2, z2, q, 3, sigmaX, sigmaY, sigmaZ);
      };
    };
  };

  // return extrapolated charge
  return projGroup;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetNearestGroup(ND::THandle<ND::TTPCVolGroup> in, int x, int y, int z, int axis, int maxDist){
  ND::THandle<ND::TTPCVolGroup> newGroup(new ND::TTPCVolGroup(fLayout));

  // create group of cells in the next plane back
  ND::THandle<ND::TTPCVolGroup> planeHits(new ND::TTPCVolGroup(fLayout));
  for(std::map<long, ND::TTPCUnitVolume*>::iterator it = in->begin(); it != in->end(); ++it){
    if(axis == 1 && it->second->GetX() == x) planeHits->AddHit(it->second);
    if(axis == 2 && it->second->GetY() == y) planeHits->AddHit(it->second);
    if(axis == 3 && it->second->GetZ() == z) planeHits->AddHit(it->second);
  };

  long startID = -1;
  // find nearest hit in plane group
  if(planeHits->size() > 0) startID = planeHits->GetNearestHit(x, y, z, maxDist);
  // build whole group recursively from that
  if(startID >= 0) RecursiveFriendBuild(startID, newGroup, planeHits, ND::TTPCConnection::path);

  return newGroup;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetNearestHitGroup(ND::THandle<ND::TTPCVolGroup> in, int x, int y, int z, int maxDist){
  ND::THandle<ND::TTPCVolGroup> grp(new ND::TTPCVolGroup(fLayout));
  // find nearest hit
  long id = in->GetNearestHit(x, y, z, maxDist);
  // add it to empty group
  grp->AddHit(in->GetHit(id));
  return grp; 
}

ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetOverlap(ND::THandle<ND::TTPCVolGroup> grp1, ND::THandle<ND::TTPCVolGroup> grp2, bool swap){
  ND::THandle<ND::TTPCVolGroup> grp(new ND::TTPCVolGroup(fLayout));

  // determine cells which appear in both grp1 and grp2 (or one and not the other, if swap=true)
  for(std::map<long, ND::TTPCUnitVolume*>::iterator it = grp1->begin(); it != grp1->end(); ++it){
    if((swap)^(grp2->Contains(it->first))){
      grp->AddHit(it->second);
    };
  };

  return grp;
}

std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetConnectedHits(ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter, bool usabilityCheck){
  return GetConnectedHits(fPrimaryHits, type, typeFilter, usabilityCheck);
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetConnectedHits(ND::THandle<ND::TTPCVolGroup> in, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter, bool usabilityCheck){
  // make vector of exclusive, connected groups
  std::vector< ND::THandle<ND::TTPCVolGroup> > conGroups;

  // make dummy group of all cells
  ND::THandle<ND::TTPCVolGroup> conHits(new ND::TTPCVolGroup(fLayout));
  if(typeFilter == ND::TTPCHitGroupings::delta){
    for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = in->begin(); hitEl != in->end(); ++hitEl){
      if(hitEl->second->GetDeltaTagged()) conHits->AddHit(hitEl->second);
    };
  }
  else if(typeFilter == ND::TTPCHitGroupings::nonDelta){
    for(std::map<long, ND::TTPCUnitVolume*>::iterator hitEl = in->begin(); hitEl != in->end(); ++hitEl){
      if(!hitEl->second->GetDeltaTagged()) conHits->AddHit(hitEl->second);
    };
  }
  else{
    conHits->AddHitMap(in->GetHitMap());
  };

  for(int i=0; i<9999; i++){
    // break out of loop if dummy becomes empty
    if(conHits->size() < 1) break;
    // start at first element in dummy
    std::map<long, ND::TTPCUnitVolume*>::iterator el = conHits->begin();

    // define new group
    ND::THandle<ND::TTPCVolGroup> newGroup(new ND::TTPCVolGroup(fLayout));
    // recursively build set of connected cells, which are removed from the dummy and added to the new group
    RecursiveFriendBuild(el->first, newGroup, conHits, type, typeFilter);

    // check to make sure this group is useful
    if(usabilityCheck){
      if(CheckUsability(newGroup)){
        conGroups.push_back(newGroup);
      };
    }
    else{
      conGroups.push_back(newGroup);
    };
  };

  return conGroups;
}

bool ND::TTPCVolGroupMan::CheckUsability(ND::THandle<ND::TTPCVolGroup> inGroup){
  unsigned int minSize = fLayout->GetMinPatternPads();

  std::set< std::pair<int, int> > patternPads;
  for(std::map<long, ND::TTPCUnitVolume*>::iterator cellEl = inGroup->begin(); cellEl != inGroup->end(); ++cellEl){
    std::pair<int, int> patternPad (cellEl->second->GetY(), cellEl->second->GetZ());
    patternPads.insert(patternPad);
  };

  // count number of pads
  return (minSize <= patternPads.size());
}

ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetNearHits(ND::THandle<ND::TTPCVolGroup> in, long id, ND::TTPCConnection::Type type, bool inclusive, bool singular, float distFilter, ND::TTPCHitGroupings::Type typeFilter, bool square){
  // try to find id in map
  std::map<long, ND::TTPCUnitVolume*>::iterator el = in->find(id);
  // return empty group if it's not found
  if(el == in->end()) return ND::THandle<ND::TTPCVolGroup>(new ND::TTPCVolGroup(fLayout));
  // otherwise get hits near to its specified volume
  ND::TTPCUnitVolume* vol = el->second;

  return GetNearHits(in, vol, type, inclusive, singular, distFilter, typeFilter, square);
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetNearHits(ND::THandle<ND::TTPCVolGroup> in, ND::TTPCUnitVolume* vol, ND::TTPCConnection::Type type, bool inclusive, bool singular, float distFilter, ND::TTPCHitGroupings::Type typeFilter, bool square){
  // set up empty group of nearby hits
  ND::THandle<ND::TTPCVolGroup> nearHits(new ND::TTPCVolGroup(fLayout)); 

  // fetch x, y and z ids and edge statuses from provided volume
  int x = vol->GetX();
  int y = vol->GetY();
  int z = vol->GetZ();

  int edgeX = 0;
  int edgeY = 0;
  int edgeZ = 0;
  edgeX = vol->GetEdgeX();
  edgeY = vol->GetEdgeY();
  edgeZ = vol->GetEdgeZ();

  int distXP;
  int distYP;
  int distZP;
  int distXN;
  int distYN;
  int distZN;
  fLayout->GetTypeDistances(distXP, distYP, distZP, type);
  fLayout->GetTypeDistances(distXN, distYN, distZN, type);

  // expand possible connection distance if on an edge
  if(edgeX != 0 && fLayout->GetJumpX()){
    int incG = fLayout->GetGapOffsetX() * fLayout->GetGapOffsetAdjacent();
    if(edgeX == 1) distXP += fLayout->GetGapOffsetX();
    else if(edgeX == -1) distXN += fLayout->GetGapOffsetX();
    distYP += incG;
    distYN += incG;
    distZP += incG;
    distZN += incG;
  };
  if(edgeY != 0 && fLayout->GetJumpY()){
    int incG = fLayout->GetGapOffsetY() * fLayout->GetGapOffsetAdjacent();
    if(edgeY == 1) distYP += fLayout->GetGapOffsetY();
    else if(edgeY == -1) distYN += fLayout->GetGapOffsetY();
    distXP += incG;
    distXN += incG;
    distZP += incG;
    distZN += incG;
  };
  if(edgeZ != 0 && fLayout->GetJumpZ()){
    int incG = fLayout->GetGapOffsetZ() * fLayout->GetGapOffsetAdjacent();
    if(edgeZ == 1) distZP += fLayout->GetGapOffsetZ();
    else if(edgeZ == -1) distZN += fLayout->GetGapOffsetZ();
    distXP += incG;
    distXN += incG;
    distYP += incG;
    distYN += incG;
  };

  // iterate over ellipsoid if needed
  for(int i=-distXN; i<=distXP; i++){
    int distYP2;
    int distYN2;
    float factX = (i<0) ? float(i*i)/float(distXN*distXN) : float(i*i)/float(distXP*distXP);
    factX = std::min(factX, (float).99);

    if (square){
      distYP2 = distYP;
      distYN2 = distYN;
    }
    else{
      distYP2 = int( distYP * std::sqrt(1. - factX) + .6);
      distYN2 = int( distYN * std::sqrt(1. - factX) + .6);
    };
    for(int j=-distYN2; j<=distYP2; j++){
      int distZP2;
      int distZN2;
      int distY2Eff = (j<0) ? distYN2 : distYP2;
      float factXY = factX;
      if(distY2Eff!=0) factXY += (float(j*j)/float(distY2Eff*distY2Eff));
      factXY = std::min(factXY, (float).99);

      if (square){
        distZP2 = distZP;
        distZN2 = distZN;
      }
      else{
        distZP2 = int( distZP * std::sqrt(1. - factXY) + .6);
        distZN2 = int( distZN * std::sqrt(1. - factXY) + .6);
      };
      for(int k=-distZN2; k<=distZP2; k++){
        // veto current cell if inclusive is set to false
        if(i==0 && j==0 && k==0 && !inclusive) continue;

        // get unique id from x, y and z id
        long newID = fLayout->SafeMash(x+i, y+j, z+k);
        // skip if id is invalud
        if(newID < 0) continue;
        ND::TTPCUnitVolume* newHit = in->GetHit(newID);
        if (!newHit) continue;

        // filter hits based on input requirements
        if(typeFilter == ND::TTPCHitGroupings::diff){
          if(newHit->GetDeltaTagged() != vol->GetDeltaTagged()) continue;
        };

        // veto cells on far side of cathode if x jump is disabled
        if(!fLayout->GetJumpX()){
          if(vol->GetSegX() != newHit->GetSegX()) continue;
        };
        // veto cells on far side of y MM gap if y jump is disabled
        if(!fLayout->GetJumpY()){
          if(vol->GetSegY() != newHit->GetSegY()) continue;
        };
        // veto cells on far side of z MM gap if z jump is disabled
        if(!fLayout->GetJumpZ()){
          if(vol->GetSegZ() != newHit->GetSegZ()) continue;
        };

        if (distFilter>-1. && distFilter>=newHit->GetFriendDist()) continue;
        if (distFilter>-1. && newHit->GetFriendDist() < 0) continue;

        // if id is contained in hit map, add the corresponding cell to the list of near hits
        bool added = nearHits->AddHit(newHit);
        if(singular && added){
          return nearHits;
        };
      };
    };
  };
  return nearHits;
}
std::vector< ND::THandle<ND::TTPCVolGroup> > ND::TTPCVolGroupMan::GetJunctionsFromPaths(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  // count vertices in input
  std::vector< ND::THandle<ND::TTPCVolGroup> > tempVertices;
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;
    // find front junction
    if(path->GetFrontIsVertex()) tempVertices.push_back(path->GetFrontHits());
    // find back junction
    if(path->GetBackIsVertex()) tempVertices.push_back(path->GetBackHits());
  };

  // add unique vertices to output
  std::vector< ND::THandle<ND::TTPCVolGroup> > vertices;
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertex1It = tempVertices.begin(); vertex1It != tempVertices.end(); ++vertex1It){
    ND::THandle<ND::TTPCVolGroup> vertex1 = *vertex1It;
    bool found = false;

    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator vertex2It = vertices.begin(); vertex2It != vertices.end(); ++vertex2It){
      ND::THandle<ND::TTPCVolGroup> vertex2 = *vertex2It;
      if(vertex1->GetID() == vertex2->GetID()){
        found = true;
        break;
      };
    };
    if(!found){
      vertices.push_back(vertex1);
    };
  };
  return vertices;
}
ND::THandle<ND::TTPCVolGroup> ND::TTPCVolGroupMan::GetUnusedHits(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  std::vector< ND::THandle<ND::TTPCVolGroup> > junctions = GetJunctionsFromPaths(paths);

  // start full
  ND::THandle<ND::TTPCVolGroup> unusedHits (new ND::TTPCVolGroup(fLayout));
  unusedHits->AddHits(fPrimaryHits);

  // remove path hits
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;
    for(std::vector<ND::TTPCPathVolume*>::iterator clusterIt = path->begin(); clusterIt != path->end(); ++clusterIt){
      ND::TTPCPathVolume* cluster = *clusterIt;
      for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = cluster->GetFriendsBegin(); volIt != cluster->GetFriendsEnd(); ++volIt){
        ND::TTPCUnitVolume* vol = *volIt;

        unusedHits->RemoveHit(vol->GetID());
      };
    };
  };
  // remove junction hits
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    ND::THandle<ND::TTPCVolGroup> junction = *junctionIt;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = junction->begin(); volEl != junction->end(); ++volEl){
      unusedHits->RemoveHit(volEl->first);
    };
  };

  return unusedHits;
}
void ND::TTPCVolGroupMan::GetDeltaNonDelta(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& nonDelta, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& delta, ND::THandle<ND::TTPCOrderedVolGroup> input){
}

void ND::TTPCVolGroupMan::AssociateUnusedWithJunctions(ND::THandle<ND::TTPCVolGroup> unused, std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  std::vector< ND::THandle<ND::TTPCVolGroup> > junctions = GetJunctionsFromPaths(paths);
  if(!junctions.size()) return;

  // add all junction hits to unused
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    ND::THandle<ND::TTPCVolGroup> junction = *junctionIt;
    if(!junction->size()) continue;

    unused->AddHits(junction);
  };

  // new junctions to add to paths
  std::vector< ND::THandle<ND::TTPCVolGroup> > newJunctions;

  // try and associate each junction with new hits
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    ND::THandle<ND::TTPCVolGroup> junction = *junctionIt;

    // break out of loop if unused becomes empty, otherwise try to add hits starting from anywhere in junction
    if(unused->size() < 1) break;

    // new junction to hold hits
    ND::THandle<ND::TTPCVolGroup> newJunction(new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()));

    // seed from every element in old junction so nothing is missed
    for(std::map<long, ND::TTPCUnitVolume*>::iterator seedEl = junction->begin(); seedEl != junction->end(); ++seedEl){
      RecursiveFriendBuild(seedEl->first, newJunction, unused);
    };

    if(!newJunction->size()) continue;

    // add to list of junctions
    newJunctions.push_back(newJunction);
  };

  // now add expanded junctions to paths
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator newJunctionIt = newJunctions.begin(); newJunctionIt != newJunctions.end(); ++newJunctionIt){
    ND::THandle<ND::TTPCVolGroup> newJunction = *newJunctionIt;
    if(!newJunction->size()) continue;

    for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

      bool frontOverlaps = GetOverlaps(newJunction, path->GetFrontHits());
      bool backOverlaps = GetOverlaps(newJunction, path->GetBackHits());

      if(frontOverlaps) path->AddFrontHits(newJunction);
      if(backOverlaps) path->AddBackHits(newJunction);
    };
  };
}

void ND::TTPCVolGroupMan::BreakLongJunctions(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >& paths){
  std::vector< ND::THandle<ND::TTPCVolGroup> > junctions = GetJunctionsFromPaths(paths);

  int xSizeThreshold = fLayout->GetXSizeThreshold();
  int pathSizeThreshold = fLayout->GetPathSizeThreshold();

  // iterate over junctions, checking for any that are long in x
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    ND::THandle<ND::TTPCVolGroup> junction = *junctionIt;

    if(junction->GetXSize() >= xSizeThreshold){
      // get total x range from junction itself
      int xMin = junction->GetXMin();
      int xMax = junction->GetXMax();

      // get x range from adjacent paths
      std::vector< std::pair<int, int> > xRanges;

      // get desired vertex extent from connected paths
      for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
        ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

        bool frontOverlaps = GetOverlaps(junction, path->GetFrontHits());
        bool backOverlaps = GetOverlaps(junction, path->GetBackHits());

        ND::TTPCPathVolume* checkCluster;
        if(frontOverlaps) checkCluster = *(path->rbegin());
        else if(backOverlaps) checkCluster = *(path->begin());
        else continue;

        // grab mins and maxes
        int juncXMin = checkCluster->GetXMin();
        int juncXMax = checkCluster->GetXMax();

        // enforce finite size within junction
        if(juncXMax < xMin){
          juncXMax = xMin;
        };
        if(juncXMin > xMax){
          juncXMin = xMax;
        };

        // add to list
        xRanges.push_back( std::make_pair(juncXMin, juncXMax) );
      };

      // if any limits are within range of max or min, expand them to cover
      for(std::vector< std::pair<int, int> >::iterator xRange1It = xRanges.begin(); xRange1It != xRanges.end(); ++xRange1It){
        if(xRange1It->first < xMin + pathSizeThreshold) xRange1It->first = xMin;
        if(xRange1It->second > xMax - pathSizeThreshold) xRange1It->second = xMax;
      };

      // merge overlapping limits together
      bool canMerge = true;
      while(canMerge){
        canMerge = false;

        for(std::vector< std::pair<int, int> >::iterator xRange1It = xRanges.begin(); xRange1It != xRanges.end(); ++xRange1It){
          for(std::vector< std::pair<int, int> >::iterator xRange2It = xRange1It+1; xRange2It != xRanges.end(); ++xRange2It){
            if( xRange1It->first < (xRange2It->second + pathSizeThreshold) && xRange1It->second > (xRange2It->first - pathSizeThreshold) ){
              // overlap found - merge
              xRange1It->first = std::min(xRange1It->first, xRange2It->first);
              xRange1It->second = std::max(xRange1It->second, xRange2It->second);

              xRanges.erase(xRange2It);
              canMerge = true;
            };

            if(canMerge) break;
          };
          if(canMerge) break;
        };
      };

      // order remaining limits
      std::vector< std::pair<int, int> > pathRanges;
      std::vector< std::pair<int, int> > junctionRanges;

      int pathMin = xMin;
      bool needOrdering = true;
      while(needOrdering){
        needOrdering = false;

        int junctionMin = 9999999;
        std::vector< std::pair<int, int> >::iterator xRangeItId = xRanges.end();
        for(std::vector< std::pair<int, int> >::iterator xRange1It = xRanges.begin(); xRange1It != xRanges.end(); ++xRange1It){
          if(xRange1It->first < junctionMin){
            junctionMin = xRange1It->first;
            xRangeItId = xRange1It;
          };
        };

        if(xRangeItId != xRanges.end()){
          pathRanges.push_back( std::make_pair(pathMin, xRangeItId->first-1) );
          junctionRanges.push_back( std::make_pair(xRangeItId->first, xRangeItId->second) );

          pathMin = xRangeItId->second+1;
          xRanges.erase(xRangeItId);

          needOrdering = true;
        };
      };
      pathRanges.push_back( std::make_pair(pathMin, xMax) );

      // if breaks in middle are disabled, clear everything apart from the first and last junction and anything between them
      if(!fLayout->GetBreakInMiddle()){
        std::pair<int, int> pathPair1 = *(pathRanges.begin());
        std::pair<int, int> pathPair2 = *(pathRanges.rbegin());
        std::pair<int, int> junctionPair = std::make_pair(pathPair1.second+1, pathPair2.first-1);

        pathRanges.clear();
        pathRanges.push_back(pathPair1);
        pathRanges.push_back(pathPair2);

        junctionRanges.clear();
        junctionRanges.push_back(junctionPair);
      };

      // make new paths and new junctions from old big one
      std::vector< ND::THandle<ND::TTPCVolGroup> > newPathHitss;
      std::vector< ND::THandle<ND::TTPCVolGroup> > newJunctions;

      // get hits from original junction and add them to new paths or new junctions depending on where they are
      for(unsigned int i=0; i<pathRanges.size(); ++i){
        newPathHitss.push_back( ND::THandle<ND::TTPCVolGroup>( new ND::TTPCVolGroup(fLayout)) );
      };
      for(unsigned int i=0; i<junctionRanges.size(); ++i){
        newJunctions.push_back( ND::THandle<ND::TTPCVolGroup>( new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()) ) );
      };

      for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = junction->begin(); volEl != junction->end(); ++volEl){
        ND::TTPCUnitVolume* vol = volEl->second;
        bool pathHit = false;
        bool junctionHit = false;

        for(unsigned int i=0; i<pathRanges.size(); ++i){
          std::pair<int, int> pathRange = pathRanges.at(i);
          if(vol->GetX() >= pathRange.first && vol->GetX() <= pathRange.second){
            newPathHitss.at(i)->AddHit(vol);
            pathHit = true;
            break;
          };
        };
        if(pathHit) continue;

        for(unsigned int i=0; i<junctionRanges.size(); ++i){
          std::pair<int, int> junctionRange = junctionRanges.at(i);
          if(vol->GetX() >= junctionRange.first && vol->GetX() <= junctionRange.second){
            newJunctions.at(i)->AddHit(vol);
            junctionHit = true;
            break;
          };
        };
        if(!junctionHit){
          // some kind of warning or error?
        };
      };

      // create new x-paths from path hits
      std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > newPaths;
      for(unsigned int i=0; i<newPathHitss.size(); ++i){
        ND::THandle<ND::TTPCVolGroup> newPathHits = newPathHitss.at(i);

        // no new path
        if(newPathHits->size() < 2){
          continue;
        };

        ND::THandle<ND::TTPCOrderedVolGroup> newPath (new ND::TTPCOrderedVolGroup(fLayout));

        ND::TTPCUnitVolume* lowHit = 0;
        ND::TTPCUnitVolume* highHit = 0;

        for(std::map<long, ND::TTPCUnitVolume*>::iterator pathHitIt = newPathHits->begin(); pathHitIt != newPathHits->end(); ++pathHitIt){
          ND::TTPCUnitVolume* pathHit = pathHitIt->second;
          if(lowHit){
            if(pathHit->GetX() < lowHit->GetX()){
              lowHit = pathHit;
            };
          }
          else{
            lowHit = pathHit;
          };
          if(highHit){
            if(pathHit->GetX() > highHit->GetX()){
              highHit = pathHit;
            };
          }
          else{
            highHit = pathHit;
          };
        };

        if(lowHit == highHit){
          continue;
        };

        newPath->AddCell(lowHit);
        newPath->AddCell(highHit);
        newPath->AddExtendedHits(newPathHits);

        if(i > 0){
          newPath->AddBackHits(newJunctions.at(i-1));
          newPath->SetBackIsVertex(true);
        }
        else{
          ND::THandle<ND::TTPCVolGroup> emptyGroup( new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()) );
          newPath->AddBackHits(emptyGroup);
          newPath->SetBackIsVertex(false);
        };
        if(i < newPathHitss.size()-1){
          newPath->AddFrontHits(newJunctions.at(i));
          newPath->SetFrontIsVertex(true);
        }
        else{
          ND::THandle<ND::TTPCVolGroup> emptyGroup( new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()) );
          newPath->AddFrontHits(emptyGroup);
          newPath->SetFrontIsVertex(false);
        };
        // ensure facing away from junction
        if(i == 0){
          newPath->Flip();
        };
        newPath->SetIsXPath(true);
        newPath->DoClustering();

        newPaths.push_back(newPath);
      };

      for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator newPathIt = newPaths.begin(); newPathIt != newPaths.end(); ++newPathIt){
        ND::THandle<ND::TTPCOrderedVolGroup> newPath = *newPathIt;
        paths.push_back(newPath);
      };

      // now add new junctions to all other paths connected to them
      for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
        ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

        bool frontOverlaps = junction->GetID() == path->GetFrontID();
        bool backOverlaps = junction->GetID() == path->GetBackID();

        ND::TTPCPathVolume* checkCluster;
        if(frontOverlaps) checkCluster = *(path->rbegin());
        else if(backOverlaps) checkCluster = *(path->begin());
        else continue;

        // overlap found - find the closest new junction in x
        ND::THandle<ND::TTPCVolGroup> closestNewJunction = ND::THandle<ND::TTPCVolGroup>();
        int closestDist = 99999999;
        for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator newJunctionIt = newJunctions.begin(); newJunctionIt != newJunctions.end(); ++newJunctionIt){
          ND::THandle<ND::TTPCVolGroup> newJunction = *newJunctionIt;

          int distX = std::abs(junction->GetX() - checkCluster->GetX());
          if(distX < closestDist){
            closestDist = distX;
            closestNewJunction = newJunction;
          };
        };

        if((frontOverlaps && backOverlaps) || !closestNewJunction){
          // if for whatever reason a close junction isn't found, need to add empty groups to prevent bugs
          if(frontOverlaps){
            ND::THandle<ND::TTPCVolGroup> emptyGroup( new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()) );
            path->AddFrontHits(emptyGroup);
            path->SetFrontIsVertex(false);
          };
          if(backOverlaps){
            ND::THandle<ND::TTPCVolGroup> emptyGroup( new ND::TTPCVolGroup(fLayout, ND::TTPCVolGroup::GetFreeID()) );
            path->AddBackHits(emptyGroup);
            path->SetBackIsVertex(false);
          };
        };
        if(closestNewJunction){
          if(frontOverlaps){
            path->AddFrontHits(closestNewJunction);
            path->SetFrontIsVertex(true);
          }
          else if(backOverlaps){
            path->AddBackHits(closestNewJunction);
            path->SetBackIsVertex(true);
          };
        };
      };
    };
  };
}

std::pair<std::pair<int, int>, std::pair<int, int> > ND::TTPCVolGroupMan::CountDuplicates(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  std::vector< ND::THandle<ND::TTPCVolGroup> > junctions = GetJunctionsFromPaths(paths);

  std::set<ND::TTPCUnitVolume*> foundHVols;
  std::set<ND::TTPCUnitVolume*> foundVVols;
  std::set<ND::TTPCUnitVolume*> foundPathVols;
  std::set<ND::TTPCUnitVolume*> foundJunctionVols;
  int nPathHDupes=0;
  int nPathVDupes=0;
  int nJunctionJunctionDupes=0;
  int nJunctionPathDupes=0;

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

    for(std::vector<ND::TTPCPathVolume*>::iterator pathVolIt = path->begin(); pathVolIt != path->end(); ++pathVolIt){
      ND::TTPCPathVolume* pathVol = *pathVolIt;
      bool isVertical = pathVol->GetIsVertical();

      std::set<ND::TTPCUnitVolume*> pathVols = std::set<ND::TTPCUnitVolume*>();
      std::vector<ND::TTPCUnitVolume*> vols = pathVol->GetExtendedCell();
      for(std::vector<ND::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
        ND::TTPCUnitVolume* vol = *volIt;
        if(isVertical){
          if(foundVVols.count(vol)) nPathVDupes++;
          foundVVols.insert(vol);
        }
        else{
          if(foundHVols.count(vol)) nPathHDupes++;
          foundHVols.insert(vol);
        };
        foundPathVols.insert(vol);
      };
    };
  };

  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    ND::THandle<ND::TTPCVolGroup> junction = *junctionIt;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = junction->begin(); volEl != junction->end(); ++volEl){
      if(foundJunctionVols.count(volEl->second)){
        nJunctionJunctionDupes ++;
      };
      if(foundPathVols.count(volEl->second)){
        nJunctionPathDupes ++;
      };
      foundJunctionVols.insert(volEl->second);
    };
  };

  return std::make_pair(std::make_pair(nPathHDupes, nPathVDupes), std::make_pair(nJunctionJunctionDupes, nJunctionPathDupes));
}

std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > ND::TTPCVolGroupMan::SanityFilter(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > input){
  std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > output = std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >();

  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator pathIt = input.begin(); pathIt != input.end(); ++pathIt){
    ND::THandle<ND::TTPCOrderedVolGroup> path = *pathIt;

    // ensure that path contains no empty clusters
    std::vector<ND::TTPCPathVolume*>::iterator reaper = path->begin();
    while(reaper != path->end()){
      if(!(*reaper)->GetClusterSize()){
        path->erase(reaper);
      }
      else{
        ++reaper;
      };
    };

    // ensure path contains above threshold of clusters
    if(path->size() < fLayout->GetMinPathClusters()) continue;
    output.push_back(path);
  };

  return output;
}

void ND::TTPCVolGroupMan::RecursiveFriendBuild(long startID, ND::THandle<ND::TTPCVolGroup> target, ND::THandle<ND::TTPCVolGroup> source, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter){
  // initialise start cell by provided id
  ND::TTPCUnitVolume* startCell = fHitMap[startID];
  // check if cell exists in source
  if(!source->Contains(startID)) return;
  // attempt to remove start cell from the group being drained
  if(!source->RemoveHit(startID)) return;
  // attempt to add start cell to the group being filled
  if(!target->AddHit(startCell)) return;

  // get near hits from the group being drained and attempt to add all of them
  ND::THandle<ND::TTPCVolGroup> nearHits = GetNearHits(source, startCell, type, false,false,-1., typeFilter);
  for(std::map<long, ND::TTPCUnitVolume*>::iterator vol = nearHits->begin(); vol != nearHits->end(); ++vol){
    // repeat the process on every near hit found
    RecursiveFriendBuild(vol->first, target, source, type, typeFilter);
  };
}
void ND::TTPCVolGroupMan::FriendConnect(ND::THandle<ND::TTPCOrderedVolGroup> inList, ND::THandle<ND::TTPCVolGroup> inField){
  /*// find parents for all field points
  for(std::map<long, ND::TTPCUnitVolume*>::iterator fieldPnt = inField->begin(); fieldPnt != inField->end(); ++fieldPnt){
    if(fieldPnt->second->GetIsFocus()) continue;
    if(fieldPnt->second->GetFriendDist() < 9999.){
      fieldPnt->second->SetFriendDist(9999);
      for(std::vector<ND::TTPCPathVolume>::iterator pathPnt = inList->begin(); pathPnt != inList->end(); ++pathPnt){

        TVector3 diffVec = fieldPnt->second->GetPos() - pathPnt->GetPos();
        double diffMag = diffVec.Mag();
        if(fieldPnt->second->GetFriendDist() > diffMag){
          fieldPnt->second->SetFriendDist(diffMag);
          fieldPnt->second->SetFocusParent(*pathPnt);
        };
      };
    };
  };
  // add information on children to their parents
  for(std::map<long, ND::TTPCUnitVolume*>::iterator fieldPnt = inField->begin(); fieldPnt != inField->end(); ++fieldPnt){
    if(fieldPnt->second->GetIsFocus()) continue;
    if(fieldPnt->second->GetFriendDist() < 9999.){
      ND::TTPCUnitVolume* focusParent = fieldPnt->second->GetFocusParent();
      if(!focusParent) continue;
      focusParent->AddFocusFriend(fieldPnt->second);
    };
  };*/
}
void ND::TTPCVolGroupMan::RecursiveFriendSeek(ND::THandle<ND::TTPCOrderedVolGroup> inList, ND::THandle<ND::TTPCVolGroup> target, float dist, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter){
  RecursiveFriendSeek(GetUnorderedGroup(inList), target, dist, type, typeFilter);
}
void ND::TTPCVolGroupMan::RecursiveFriendSeek(ND::THandle<ND::TTPCVolGroup> inList, ND::THandle<ND::TTPCVolGroup> target, float dist, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter){
  for(std::map<long, ND::TTPCUnitVolume*>::iterator pnt = fPrimaryHits->begin(); pnt != fPrimaryHits->end(); ++pnt){
    pnt->second->SetFriendDist(9999);
  };
  for(std::map<long, ND::TTPCUnitVolume*>::iterator pnt = inList->begin(); pnt != inList->end(); ++pnt){
    pnt->second->SetFriendDist(0);
    RecursiveFriendListSeek(pnt->first, fPrimaryHits, 0, dist, type, typeFilter);
  };

  for(std::map<long, ND::TTPCUnitVolume*>::iterator pnt = fPrimaryHits->begin(); pnt != fPrimaryHits->end(); ++pnt){
    if(pnt->second->GetFriendDist() < 9998){
      target->AddHit(pnt->second);
    };
  };
}
void ND::TTPCVolGroupMan::RecursiveFriendListSeek(long startID, ND::THandle<ND::TTPCVolGroup> source, float curDist, float dist, ND::TTPCConnection::Type type, ND::TTPCHitGroupings::Type typeFilter){
  // break if over or at max distance
  if (curDist > dist) return;

  // set current cell's distance based on current input distance
  ND::TTPCUnitVolume* startCell = fHitMap[startID];
  startCell->SetFriendDist(curDist);

  // get near hits from the group being drained and attempt to add all of them
  ND::THandle<ND::TTPCVolGroup> nearHits = GetNearHits(source, startCell, type, false, false, dist, typeFilter);
  for(std::map<long, ND::TTPCUnitVolume*>::iterator vol = nearHits->begin(); vol != nearHits->end(); ++vol){
    // repeat the process on every near hit found
    RecursiveFriendListSeek(vol->first, source, curDist+1., dist, type, typeFilter);
  };
}

std::vector<long> ND::TTPCVolGroupMan::SpheroidIterable(ND::TTPCCell3D cell, ND::TTPCConnection::Type type){
  int sizeX;
  int sizeY;
  int sizeZ;
  fLayout->GetTypeDistances(sizeX, sizeY, sizeZ, type);

  return SpheroidIterable(cell, sizeX, sizeY, sizeZ);
}
std::vector<long> ND::TTPCVolGroupMan::SpheroidIterable(ND::TTPCCell3D cell, int diffX, int diffY, int diffZ){
  std::vector<long> spheroid;

  for(int i=-diffX; i<=diffX; i++){
    int diffY2 = int( std::sqrt(diffY*diffY - i*i) + .6);
    for(int j=-diffY2; j<=diffY2; j++){
      int diffZ2 = int( std::sqrt(diffZ*diffZ - i*i - j*j) + .6);
      for(int k=-diffZ2; k<=diffZ2; k++){
        long newID = fLayout->SafeMash(cell.x+i, cell.y+j, cell.z+k);
        if(newID < 0) continue;

        spheroid.push_back(newID);
      };
    };
  };
  return spheroid;
}
TVector3 ND::TTPCVolGroupMan::GetAvgPosRep(ND::TTPCPathVolume* vol){
  double avgTime = vol->GetAverageTime();
  TVector3 avgPos = vol->GetAveragePos();
  return TVector3(avgTime*fLayout->GetDriftSpeed(), avgPos.Y(), avgPos.Z());
}
TVector3 ND::TTPCVolGroupMan::GetAvgPosRep(ND::TTPCUnitVolume* vol, int sign){
  double avgTime;
  if(sign < 0){
    avgTime = vol->GetTimeMin();
  }
  else if(sign > 0){
    avgTime = vol->GetTimeMax();
  }
  else{
    avgTime = vol->GetTime();
  };
  TVector3 avgPos = vol->GetPos();

  return TVector3(avgTime*fLayout->GetDriftSpeed(), avgPos.Y(), avgPos.Z());
}

bool ND::TTPCVolGroupMan::IsInRange(ND::TTPCPathVolume* point1, ND::TTPCPathVolume* point2, int sizeX, int sizeY, int sizeZ){
  float dx = (float)(point1->GetX() - point2->GetX()) / float(sizeX);
  float dy = (float)(point1->GetY() - point2->GetY()) / float(sizeY);
  float dz = (float)(point1->GetZ() - point2->GetZ()) / float(sizeZ);

  return (dx*dx + dy*dy + dz*dz) < 1.;
  /*if( abs(point1->GetX() - point2->GetX()) > sizeX) return false;
  if( abs(point1->GetY() - point2->GetY()) > sizeY) return false;
  if( abs(point1->GetZ() - point2->GetZ()) > sizeZ) return false;
  return true;*/
}
bool ND::TTPCVolGroupMan::IsXPathCandidate(ND::THandle<ND::TTPCOrderedVolGroup> inPath){
  // get status as x cluster
  // first check; is cluster too small in y or z?
  int minClusters = fLayout->GetMinPathClusters();
  int ySize = (inPath->GetYMax() - inPath->GetYMin()) + 1;
  int zSize = (inPath->GetZMax() - inPath->GetZMin()) + 1;
  bool sizeCheck1 = (ySize < minClusters) && (zSize < minClusters);

  if(sizeCheck1) return true;

  int maxCheckClusters = fLayout->GetXPathMaxPads();
  double minEndRatio = fLayout->GetXPathMinEndRatio();

  bool sizeCheck2 = (ySize < maxCheckClusters) && (zSize < maxCheckClusters);

  ND::TTPCPathVolume* beginVol = *(inPath->begin());
  ND::TTPCPathVolume* endVol = *(inPath->rbegin());

  double xDiff = endVol->GetX() - beginVol->GetX();
  double yDiff = endVol->GetY() - beginVol->GetY();
  double zDiff = endVol->GetZ() - beginVol->GetZ();

  double  ratio2 = (xDiff*xDiff) / ( (yDiff*yDiff) + (zDiff*zDiff) );

  bool ratioCheck = ratio2 > (minEndRatio*minEndRatio);

  if(sizeCheck2 && ratioCheck) return true;

  return false;
}


ND::THandle<ND::THitSelection> ND::TTPCVolGroupMan::GetHits(){
  return fPrimaryHits->GetHits();
}

void ND::TTPCVolGroupMan::CheckExtendedHitOverlap(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> > paths){
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It)
  for(std::vector< ND::THandle<ND::TTPCOrderedVolGroup> >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
    ND::THandle<ND::TTPCOrderedVolGroup> path1 = *path1It;
    ND::THandle<ND::TTPCOrderedVolGroup> path2 = *path2It;

    for(std::map<long, ND::TTPCUnitVolume*>::iterator path1HitsEl = path1->extendedHitsBegin(); path1HitsEl != path1->extendedHitsEnd(); ++path1HitsEl)
    for(std::map<long, ND::TTPCUnitVolume*>::iterator path2HitsEl = path2->extendedHitsBegin(); path2HitsEl != path2->extendedHitsEnd(); ++path2HitsEl){
      if(path1HitsEl->second == path2HitsEl->second){
        std::cout << "overlap between extended hits!!!" << std::endl;
      };
    };
  };
};

bool ND::TTPCVolGroupMan::HitTest(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, ND::TTPCConnection::Type type){
  int distX;
  int distY;
  int distZ;
  fLayout->GetTypeDistances(distX, distY, distZ, type);

  return HitTest(vol1, vol2, distX, distY, distZ);
}
bool ND::TTPCVolGroupMan::HitTest(ND::TTPCUnitVolume* vol1, ND::TTPCUnitVolume* vol2, int xDist, int yDist, int zDist){
  float fx = (float)(vol2->GetX() - vol1->GetX())/(float)xDist;
  float fy = (float)(vol2->GetY() - vol1->GetY())/(float)yDist;
  float fz = (float)(vol2->GetZ() - vol1->GetZ())/(float)zDist;

  return (fx*fx + fy*fy + fz*fz <= 1.);
}
