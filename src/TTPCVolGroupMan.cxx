// eddy
#include "TTPCVolGroupMan.hxx"

trex::TTPCVolGroupMan::TTPCVolGroupMan(trex::TTPCLayout* layout):
  fLayout(layout),fPrimaryHits(layout){}


void trex::TTPCVolGroupMan::AddPrimaryHits(std::map<long, trex::TTPCUnitVolume*>& hitMap){
  fHitMap = hitMap;
  fPrimaryHits.AddHitMap(hitMap);
}

std::vector<trex::TTPCVolGroup> trex::TTPCVolGroupMan::GetAllEdgeGroups(){
  std::vector< trex::TTPCVolGroup > allGroups;
  std::vector< trex::TTPCVolGroup > nonDeltaGroups;
  GetConnectedHits(nonDeltaGroups,trex::TTPCConnection::path, trex::TTPCHitGroupings::nonDelta);

  // add the non-delta groups as one group
  allGroups.emplace_back(fLayout);
  trex::TTPCVolGroup& nonDeltaGroup=allGroups.back();
  for(std::vector< trex::TTPCVolGroup >::iterator grpIt = nonDeltaGroups.begin(); grpIt != nonDeltaGroups.end(); ++grpIt){
    nonDeltaGroup.AddHits(*grpIt);
  }

  std::vector<trex::TTPCVolGroup> output;

  for(std::vector< trex::TTPCVolGroup >::iterator groupIt = allGroups.begin(); groupIt != allGroups.end(); ++groupIt){
    std::vector< trex::TTPCVolGroup > edgeGroups = GetEdgeGroups(*groupIt, true);
    for(std::vector< trex::TTPCVolGroup >::iterator edgeIt = edgeGroups.begin(); edgeIt != edgeGroups.end(); ++edgeIt){
      trex::TTPCVolGroup& edge = *edgeIt;
      output.push_back(edge);
    };
  };

  for(std::vector< trex::TTPCVolGroup >::iterator grpCheckIt = output.begin(); grpCheckIt != output.end(); ++grpCheckIt){
    trex::TTPCVolGroup& grpCheck = *grpCheckIt;
  };

  return std::move(output);
}

std::vector< trex::TTPCVolGroup > trex::TTPCVolGroupMan::GetEdgeGroups(){
  std::vector< trex::TTPCVolGroup > out = GetEdgeGroups(fPrimaryHits);
  return std::move(out);
}

std::vector< trex::TTPCVolGroup > trex::TTPCVolGroupMan::GetEdgeGroups(trex::TTPCVolGroup& inGroup, bool fiddlyLeans){

  bool haveX=fLayout->GetHaveX();
  int layers = fLayout->GetEdgeLayers();
  trex::TTPCConnection::Type type = trex::TTPCConnection::path;

  int xLeanOverride = 0;
  int yLeanOverride = 0;
  int zLeanOverride = 0;
  // special checks for groups that're too small to get an edge for
  if(inGroup.GetXSize() < 2*layers){
    if(inGroup.GetXMax() == fPrimaryHits.GetXMax()) xLeanOverride += 1;
    if(inGroup.GetXMin() == fPrimaryHits.GetXMin()) xLeanOverride -= 1;
  };
  if(inGroup.GetYSize() < 2*layers){
    if(inGroup.GetYMax() == fPrimaryHits.GetYMax()) yLeanOverride += 1;
    if(inGroup.GetYMin() == fPrimaryHits.GetYMin()) yLeanOverride -= 1;
  };
  if(inGroup.GetZSize() < 2*layers){
    if(inGroup.GetZMax() == fPrimaryHits.GetZMax()) zLeanOverride += 1;
    if(inGroup.GetZMin() == fPrimaryHits.GetZMin()) zLeanOverride -= 1;
  };

  // get groups of all hits on x, y and z edges
  trex::TTPCVolGroup edgeHitsXHi(fLayout);
  trex::TTPCVolGroup edgeHitsXLo(fLayout);
  trex::TTPCVolGroup edgeHitsYHi(fLayout);
  trex::TTPCVolGroup edgeHitsYLo(fLayout);
  trex::TTPCVolGroup edgeHitsZHi(fLayout);
  trex::TTPCVolGroup edgeHitsZLo(fLayout);

  for(std::map<long, trex::TTPCUnitVolume*>::iterator vol = inGroup.begin(); vol != inGroup.end(); ++vol){
    // fetch x, y and z id and edge status for current cell
    int x = vol->second->GetX();
    int y = vol->second->GetY();
    int z = vol->second->GetZ();

    // add hits to relevant group, in x only in the case of deltas 
    if(haveX && (x > (inGroup.GetXMax() - layers)) ) edgeHitsXHi.AddHit(vol->second);
    if(haveX && (x < (inGroup.GetXMin() + layers)) ) edgeHitsXLo.AddHit(vol->second);
    if(y > (inGroup.GetYMax() - layers) ) edgeHitsYHi.AddHit(vol->second);
    if(y < (inGroup.GetYMin() + layers) ) edgeHitsYLo.AddHit(vol->second);
    if(z > (inGroup.GetZMax() - layers) ) edgeHitsZHi.AddHit(vol->second);
    if(z < (inGroup.GetZMin() + layers) ) edgeHitsZLo.AddHit(vol->second);
  };

  std::cout<<"Edge hit sizes: "<<edgeHitsXHi.size()<<", "<<edgeHitsYHi.size()<<", "<<edgeHitsZHi.size()<<", "
	   <<edgeHitsXLo.size()<<", "<<edgeHitsYLo.size()<<", "<<edgeHitsZLo.size()<<std::endl;

  std::vector< trex::TTPCVolGroup > edgeGroups;
  // recursively build groups of hits on x, y and z edges
  if(fiddlyLeans){
    FillWithSplitGroups(edgeGroups, edgeHitsXHi, type, 1*(inGroup.GetXMax() == fPrimaryHits.GetXMax()) ,0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsXLo, type, -1*(inGroup.GetXMin() == fPrimaryHits.GetXMin()),0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYHi, type, 0,1*(inGroup.GetYMax() == fPrimaryHits.GetYMax()) ,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYLo, type, 0,-1*(inGroup.GetYMin() == fPrimaryHits.GetYMin()),0);
    FillWithSplitGroups(edgeGroups, edgeHitsZHi, type, 0,0,1*(inGroup.GetZMax() == fPrimaryHits.GetZMax()) );
    FillWithSplitGroups(edgeGroups, edgeHitsZLo, type, 0,0,-1*(inGroup.GetZMin() == fPrimaryHits.GetZMin()));
  }
  else{
    FillWithSplitGroups(edgeGroups, edgeHitsXHi, type, 1 ,0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsXLo, type, -1,0,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYHi, type, 0,1 ,0);
    FillWithSplitGroups(edgeGroups, edgeHitsYLo, type, 0,-1,0);
    FillWithSplitGroups(edgeGroups, edgeHitsZHi, type, 0,0,1 );
    FillWithSplitGroups(edgeGroups, edgeHitsZLo, type, 0,0,-1);
  };

  std::cout<<"Ingroup min and max at: ("
	   <<inGroup.GetXMin()<<", "<<inGroup.GetYMin()<<", "<<inGroup.GetZMin()<<"), ("
	   <<inGroup.GetXMax()<<", "<<inGroup.GetYMax()<<", "<<inGroup.GetZMax()<<") "
	   <<std::endl;

  std::cout<<"GetEdgeGroups finds groups at: "<<std::endl;
  for(auto iGrp=edgeGroups.begin();iGrp!=edgeGroups.end();++iGrp){
    TVector3 pos=iGrp->GetAveragePosition();
    std::cout<<"  "<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<std::endl;
  }

  // loop over all found groups on x, y and z edges
  for(std::vector< trex::TTPCVolGroup >::iterator grp1It = edgeGroups.begin(); grp1It != edgeGroups.end(); ++grp1It){
    trex::TTPCVolGroup& grp1 = *grp1It;
    // if group overlaps with another, delete the biggest (break ties in favour of z, then y)
    // smallest group is most likely to be close to end, of e.g. low angle track
    for(std::vector< trex::TTPCVolGroup >::iterator grp2It = grp1It+1; grp2It != edgeGroups.end(); ++grp2It){
      trex::TTPCVolGroup& grp2 = *grp2It;
      if(grp2.empty()) continue;
      for(std::map<long, trex::TTPCUnitVolume*>::iterator grp2El = grp2.begin(); grp2El != grp2.end(); ++grp2El){
        if(grp1.Contains(grp2El->first)){
          if(grp2.size() < grp1.size()){
            grp1.SetXLean(grp2.GetXLean());
            grp1.SetYLean(grp2.GetYLean());
            grp1.SetZLean(grp2.GetZLean());
            grp1.clear();
          }
          else{
            grp2.SetXLean(grp1.GetXLean());
            grp2.SetYLean(grp1.GetYLean());
            grp2.SetZLean(grp1.GetZLean());
            grp2.clear();
          };
          break;
        };
      };
    };
  };
  // clean up empty groups
  std::vector< trex::TTPCVolGroup >::iterator grpDel = edgeGroups.begin();
  while(grpDel != edgeGroups.end()){
    if(grpDel->empty()) edgeGroups.erase(grpDel);
    else grpDel ++;
  };
  // override leans
  for(std::vector< trex::TTPCVolGroup >::iterator grpIt = edgeGroups.begin(); grpIt != edgeGroups.end(); ++grpIt){
    trex::TTPCVolGroup& grp = *grpIt;
    if(xLeanOverride) grp.SetXLean(xLeanOverride, true);
    if(yLeanOverride) grp.SetYLean(yLeanOverride, true);
    if(zLeanOverride) grp.SetZLean(zLeanOverride, true);
  };
  
  //MDH
  //This code looks like it doesn't do anything?!?
  //for(std::vector< trex::THandle<trex::TTPCVolGroup> >::iterator grp1It = edgeGroups.begin(); grp1It != edgeGroups.end(); ++grp1It){
  // trex::THandle<trex::TTPCVolGroup> grp1 = *grp1It;
  //};

  return std::move(edgeGroups);
}

void trex::TTPCVolGroupMan::FillWithSplitGroups(std::vector< trex::TTPCVolGroup >& container, trex::TTPCVolGroup& inputHits, trex::TTPCConnection::Type type, int maxFilterX, int maxFilterY, int maxFilterZ){
  //  std::vector< trex::THandle<trex::TTPCVolGroup> > splitGroups = GetSplitGroups(inputHits, type);
  std::vector< trex::TTPCVolGroup > containerTmp=GetSplitGroups(inputHits, type);
  for(std::vector< trex::TTPCVolGroup >::iterator grp = containerTmp.begin(); grp != containerTmp.end(); ++grp){
    if(maxFilterX) grp->SetXLean(maxFilterX);
    if(maxFilterY) grp->SetYLean(maxFilterY);
    if(maxFilterZ) grp->SetZLean(maxFilterZ);
    container.push_back(*grp);
  }
}

std::vector< trex::TTPCVolGroup > trex::TTPCVolGroupMan::GetSplitGroups(trex::TTPCVolGroup& inputHits, trex::TTPCConnection::Type type){
  std::vector< trex::TTPCVolGroup > groups;

  for(int i=0; i<9999; i++){
    if(inputHits.size() < 1) break;
    std::map<long, trex::TTPCUnitVolume*>::iterator el = inputHits.begin();

    groups.emplace_back(fLayout, trex::TTPCVolGroup::GetFreeID());
    trex::TTPCVolGroup& newGroup=groups.back();
    RecursiveFriendBuild(el->first, newGroup, inputHits, type);

    if(newGroup.size() == 0) groups.pop_back();
  };
  return std::move(groups);
}


//MDH
//Just changed sig for now - need to fix internals 
void trex::TTPCVolGroupMan::GetFoci(std::vector< trex::TTPCOrderedVolGroup >& paths, std::vector< trex::TTPCVolGroup >& groupsOut, float diffThreshold){
  // select pairs of groups to loop over
  std::vector< std::pair< trex::TTPCOrderedVolGroup*, trex::TTPCOrderedVolGroup* > > groupPairs;
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It){
    trex::TTPCOrderedVolGroup& path1 = *path1It;
    for(std::vector< trex::TTPCOrderedVolGroup >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
      trex::TTPCOrderedVolGroup& path2 = *path2It;

      // make sure paths are populated
      if(path1.size() < 1 || path2.size() < 1) continue;

      // make sure paths start in same place
      trex::TTPCUnitVolume* beginVol1 = (*(path1.begin()))->GetUnitVolume();
      trex::TTPCUnitVolume* beginVol2 = (*(path2.begin()))->GetUnitVolume();
      if(beginVol1 != beginVol2) continue;

      bool startAdded = false;
      // loop and make sure only the best pair is chosen for each start point
      for(std::vector< std::pair< trex::TTPCOrderedVolGroup*, trex::TTPCOrderedVolGroup* > >::iterator groupPairIt = groupPairs.begin(); groupPairIt != groupPairs.end(); ++groupPairIt){
        trex::TTPCOrderedVolGroup& pairGrp1 = *(groupPairIt->first);
        trex::TTPCOrderedVolGroup& pairGrp2 = *(groupPairIt->second);

        // do they start in the same place?
        if( (*(pairGrp1.begin()))->GetUnitVolume() != beginVol1) continue;
        startAdded = true;

        // end volumes
        trex::TTPCUnitVolume* endVol1 = (*(path1.rbegin()))->GetUnitVolume();
        trex::TTPCUnitVolume* endVol2 = (*(path2.rbegin()))->GetUnitVolume();
        trex::TTPCUnitVolume* endVol3 = (*(pairGrp1.rbegin()))->GetUnitVolume();
        trex::TTPCUnitVolume* endVol4 = (*(pairGrp2.rbegin()))->GetUnitVolume();

        // compare triangle areas to work out best one
        double area1 = GetTriangleArea(beginVol1, endVol1, endVol2);
        double area2 = GetTriangleArea(beginVol1, endVol3, endVol4);

        // replace old with new if new gives a bigger triangle
        if(area1 > area2){
          *groupPairIt = std::pair< trex::TTPCOrderedVolGroup*, trex::TTPCOrderedVolGroup* >(&path1, &path2);
        };
      };

      // push back if start not already found
      if(!startAdded) groupPairs.push_back( std::pair<trex::TTPCOrderedVolGroup*, trex::TTPCOrderedVolGroup* >(&path1, &path2) );
    };
  };

  // groups for holding hits at points of divergence
  trex::TTPCVolGroup grp(fLayout);
  // map for counting how many groups average at a given point
  std::map<long, int> avgPoints;
  for(std::vector< std::pair< trex::TTPCOrderedVolGroup*, trex::TTPCOrderedVolGroup* > >::iterator groupPairIt = groupPairs.begin(); groupPairIt != groupPairs.end(); ++groupPairIt){
    trex::TTPCOrderedVolGroup& path1 = *(groupPairIt->first);
    trex::TTPCOrderedVolGroup& path2 = *(groupPairIt->second);

    // get detatchment points
    trex::TTPCPathVolume* detatch1 = GetDetatchmentPoint(path1, path2, trex::TTPCConnection::vertexFind);
    trex::TTPCPathVolume* detatch2 = GetDetatchmentPoint(path2, path1, trex::TTPCConnection::vertexFind);

    // find the closest points between the two lines
    trex::TTPCUnitVolume* closestPoint = GetZero(path1, path2, detatch1, detatch2);

    // get closest point from projection points (defunct)
    //trex::TTPCPathVolume* project1 = GetProjectionPoint(path1, detatch1, fLayout->GetExtendedPathDist());
    //trex::TTPCPathVolume* project2 = GetProjectionPoint(path2, detatch2, fLayout->GetExtendedPathDist());
    //trex::TTPCUnitVolume* closestPoint = GetClosestPoint(project1, project2, detatch1, detatch2);

    // add those points and save number in each volume for future calculations
    if(closestPoint){
      grp.AddHit(closestPoint);
      if(!avgPoints.count(closestPoint->GetID())) avgPoints[closestPoint->GetID()] = 1;
      else avgPoints[closestPoint->GetID()]++;
    };
  };

  // break points of divergence into separate groups
  std::vector< trex::TTPCVolGroup > groups;
  GetConnectedHits(grp,groups, trex::TTPCConnection::vertexMerge);

  int maxSize = -999999;
  // find group of biggest size
  for(std::vector< trex::TTPCVolGroup >::iterator grpIt = groups.begin(); grpIt != groups.end(); ++grpIt){
    maxSize = std::max(maxSize, int(grpIt->size()));
  };
  // set cutoff for group size at fraction of maxSize defined by threshold
  float minSize = float(maxSize) * diffThreshold;

  // build list of groups of size above the threshold
  for(std::vector< trex::TTPCVolGroup >::iterator grpIt = groups.begin(); grpIt != groups.end(); ++grpIt){
    trex::TTPCVolGroup& grp = *grpIt;
    // find weighted average position of hit in group
    if(float(grp.size()) >= minSize){
      groupsOut.emplace_back(fLayout, trex::TTPCVolGroup::GetFreeID());
      trex::TTPCVolGroup& newGrp = groupsOut.back();
      int avgX = 0;
      int avgY = 0;
      int avgZ = 0;
      int weightSum = 0;
      for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = grp.begin(); volEl != grp.end(); ++volEl){
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
      long id = grp.GetNearestHit(avgX, avgY, avgZ);
      newGrp.AddHit(grp.GetEl(id));

      // also add nearby hits
      trex::TTPCVolGroup nearHits(fLayout);
      GetNearHits(fPrimaryHits, nearHits, id, trex::TTPCConnection::vertexHit);
      newGrp.AddHits(nearHits);

    };
  };

  // return those groups
}

trex::TTPCUnitVolume* trex::TTPCVolGroupMan::GetZero(trex::TTPCOrderedVolGroup& path1, trex::TTPCOrderedVolGroup& path2, trex::TTPCPathVolume* vol1, trex::TTPCPathVolume* vol2){
  trex::TTPCPathVolume* nearVol1 = GetPathZero(path1, path2, vol1);
  trex::TTPCPathVolume* nearVol2 = GetPathZero(path2, path1, vol2);

  if(nearVol1 && nearVol2){
    // average
    int ax = (int)((nearVol2->GetX() + nearVol1->GetX())*.5);
    int ay = (int)((nearVol2->GetY() + nearVol1->GetY())*.5);
    int az = (int)((nearVol2->GetZ() + nearVol1->GetZ())*.5);

    long id = fPrimaryHits.GetNearestHit(ax, ay, az);
    return fPrimaryHits.GetHit(id);
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

    long id = fPrimaryHits.GetNearestHit(ax, ay, az);
    return fPrimaryHits.GetHit(id);
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
trex::TTPCPathVolume* trex::TTPCVolGroupMan::GetPathZero(trex::TTPCOrderedVolGroup& path1, trex::TTPCOrderedVolGroup& path2, trex::TTPCPathVolume* vol){
  if(!vol) vol = *(path1.rbegin());

  // start checking after first points fed in
  bool go=false;

  // initial distances and number of iterations
  double initDist = -1.;
  int nIterations = 0;

  trex::TTPCPathVolume* nearVol = 0;

  // check backwards along first path
  for(std::vector<trex::TTPCPathVolume*>::reverse_iterator pathVolRit = path1.rbegin(); pathVolRit != path1.rend(); ++pathVolRit){
    trex::TTPCPathVolume* pathVol = *pathVolRit;
    if(go){
      if(nIterations<1){
        initDist = GetMinDistance(path2, pathVol);
      }
      else{
        double curDist = GetMinDistance(path2, pathVol);
        // if more than half way to the zero, project
        if( (curDist/initDist) / .5){
          while(nIterations > 0 && pathVolRit != (path1.rend()-1)){
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

double trex::TTPCVolGroupMan::GetTriangleArea(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, trex::TTPCUnitVolume* vol3){
  TVector3 sideA = vol2->GetPos() - vol1->GetPos();
  TVector3 sideB = vol3->GetPos() - vol1->GetPos();

  double a2 = sideA.Mag2();
  double b2 = sideB.Mag2();
  double ab = std::sqrt(a2*b2);
  double sinC = sideA.Angle(sideB);
  return .5*ab*sinC;
}
double trex::TTPCVolGroupMan::GetMinDistance(trex::TTPCOrderedVolGroup& path, trex::TTPCPathVolume* point){
  return GetMinDistance(path, point->GetUnitVolume());
}
double trex::TTPCVolGroupMan::GetMinDistance(trex::TTPCOrderedVolGroup& path, trex::TTPCUnitVolume* vol){
  double minDist2 = 99999999.;
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;

    double dist2 = (pathVol->GetPosXYZ() - vol->GetPosXYZ()).Mag2();
    minDist2 = std::min(minDist2, dist2);
  };
  return std::sqrt(minDist2);
}

trex::TTPCPathVolume* trex::TTPCVolGroupMan::GetDetatchmentPoint(trex::TTPCOrderedVolGroup& path1, trex::TTPCOrderedVolGroup& path2, trex::TTPCConnection::Type type){
  int sizeX;
  int sizeY;
  int sizeZ;
  fLayout->GetTypeDistances(sizeX, sizeY, sizeZ, trex::TTPCConnection::vertexFind);

  // check most likely point first to save time
  int prevPoint=0;
  int prevMax=path2.size();

  for(std::vector<trex::TTPCPathVolume*>::iterator point1It = path1.begin(); point1It != path1.end(); ++point1It){
    bool wasFound = false;

    trex::TTPCPathVolume* point1 = *point1It;
    trex::TTPCPathVolume* point2;

    int prevOld=prevPoint;
    int prevOff=0;

    // check same index first to save time, since this should usually be true if anything is
    while(prevOld+prevOff < prevMax || prevOld-prevOff >= 0){
      for (int sign=-1; sign<=1; sign+=2){
        int dummy = prevOld+sign*prevOff;
        if(dummy < 0) continue;
        if(dummy >= prevMax) continue;
        prevPoint = dummy;
        point2 = path2[prevPoint];
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
  return *(path1.end()-1);
}
trex::TTPCPathVolume* trex::TTPCVolGroupMan::GetProjectionPoint(trex::TTPCOrderedVolGroup& path, trex::TTPCPathVolume* pathVol, int dist){
  // get projection point
  bool foundDetatch = false;
  for(std::vector<trex::TTPCPathVolume*>::reverse_iterator pointRit = path.rbegin(); pointRit != path.rend(); ++pointRit){
    trex::TTPCPathVolume* point=*pointRit;
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
  return *(path.rend()-1);
}
trex::TTPCUnitVolume* trex::TTPCVolGroupMan::GetClosestPoint(trex::TTPCPathVolume* begin1, trex::TTPCPathVolume* end1, trex::TTPCPathVolume* begin2, trex::TTPCPathVolume* end2){
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

  long id = fPrimaryHits.GetNearestHit(ax, ay, az);
  return fPrimaryHits.GetHit(id);
}
void trex::TTPCVolGroupMan::CleanUpVertices(std::vector< trex::TTPCVolGroup >& edgeGroups, std::vector< trex::TTPCVolGroup >& vertices){
  std::cout<<"Cleaning vertices!"<<std::endl;
  while(vertices.size() > edgeGroups.size()-2){
    // get rid of spurious vertices, one at a time
    std::vector< trex::TTPCVolGroup >::iterator markedIt;
    double minDistance = 99999999.;

    // find closest vertex to any edge group and mark it as spurious
    for(std::vector< trex::TTPCVolGroup >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
      trex::TTPCVolGroup& vertex = *vertexIt;

      for(std::vector< trex::TTPCVolGroup >::iterator edgeIt = edgeGroups.begin(); edgeIt != edgeGroups.end(); ++edgeIt){
        trex::TTPCVolGroup& edge = *edgeIt;
        TVector3 diff = edge.GetAveragePosXYZ() - vertex.GetAveragePosXYZ(); 
        double dist = diff.Mag2();

        if(dist < minDistance){
          minDistance = dist;
          markedIt = vertexIt;
        };
      };
    };
    std::cout<<"Removing vertex!"<<std::endl;
    vertices.erase(markedIt);
  };
}

void trex::TTPCVolGroupMan::GetFarHitsPreGroup(trex::TTPCOrderedVolGroup& path, trex::TTPCVolGroup& hits){
  double edgePreDist = fLayout->GetEdgePreDist();
  double edgePreDist2 = edgePreDist*edgePreDist;
  double edgePreAng = fLayout->GetEdgePreAng();

  trex::TTPCPathVolume* kinkPoint = 0;
  double kinkAngle = 360.;

  // look for local kinks around each point
  for(int i=0; i<path.size(); ++i){
    trex::TTPCPathVolume* pathVol = path.at(i);

    // look for back and forward check points
    trex::TTPCPathVolume* backCheck = 0;
    trex::TTPCPathVolume* frontCheck = 0;
    for(int j=i; j>=0; --j){
      trex::TTPCPathVolume* backCheckVol = path.at(j);
      TVector3 posDiff = backCheckVol->GetPosXYZ() - pathVol->GetPosXYZ();

      if(posDiff.Mag2() >= edgePreDist2){
        backCheck = backCheckVol;
        break;
      };
    };
    for(int j=i; j<path.size(); ++j){
      trex::TTPCPathVolume* frontCheckVol = path.at(j);
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
    hits.AddHit(kinkPoint->GetUnitVolume());
  };

}

void trex::TTPCVolGroupMan::GetFarHitsGroup(trex::TTPCOrderedVolGroup& path, trex::TTPCVolGroup& hits){
  trex::TTPCPathVolume* startVol = *(path.begin());
  trex::TTPCPathVolume* endVol = *(path.rbegin());

  TVector3 startPos = startVol->GetPosXYZ();
  TVector3 endPos = endVol->GetPosXYZ();
  TVector3 parNorm = (startPos - endPos).Unit();
  float parPos = (startPos - endPos).Mag();

  // build group to look along when finding kinks
  trex::TTPCVolGroup inGroup(fLayout);
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    trex::TTPCUnitVolume* vol = pathVol->GetUnitVolume();

    inGroup.AddHit(vol);
  };
  // build group to look along when finding maximum hit
  trex::TTPCVolGroup extendedGroup(fLayout);
  extendedGroup = path.GetExtendedHits();

  // find single path max dist
  float maxDist2 = -1.;

  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = inGroup.begin(); el != inGroup.end(); ++el){
    TVector3 diff = startPos - el->second->GetPosXYZ();
    TVector3 parallel = parNorm*diff.Dot(parNorm);
    TVector3 dist = diff - parallel;

    float dist2 = dist.Mag2();
    if(dist2 > maxDist2){
      maxDist2 = dist2;
    };
  };

  if (maxDist2 < 0.) return;

  float threshDist2 = maxDist2 * fLayout->GetEdgeRange()*fLayout->GetEdgeRange();
  float maxDist = std::sqrt(maxDist2);

  // find extended group max dist
  float extendedMaxDist2 = -1.;
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = extendedGroup.begin(); el != extendedGroup.end(); ++el){
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
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = inGroup.begin(); el != inGroup.end(); ++el){
    TVector3 diff = startPos - el->second->GetPosXYZ();
    TVector3 parallel = parNorm*diff.Dot(parNorm);
    TVector3 dist = diff - parallel;

    double distMag2 = dist.Mag2();
    if ((float)distMag2 > threshDist2){
      nAbove ++;
    };
  };

  float threshFrac = (float)nAbove / (float)inGroup.size();
  if(!threshFrac) return;

  // vector for parallel and perpendicular positions of all hits above extended threshold
  std::vector<float> parPositions;
  std::vector<float> perpPositions;
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el = extendedGroup.begin(); el != extendedGroup.end(); ++el){
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
    //remember hits has to be created in the calling function with a unique ID
    //trex::THandle<trex::TTPCVolGroup> outGroup (new trex::TTPCVolGroup(fLayout, trex::TTPCVolGroup::GetFreeID()) );

    trex::TTPCUnitVolume* furthestHit = 0;
    double furthestDist = -1;
    for(std::map<long, trex::TTPCUnitVolume*>::iterator el = extendedGroup.begin(); el != extendedGroup.end(); ++el){
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
      hits.AddHit(furthestHit);
    };
    return;
  };

  // default
  return;
}

void trex::TTPCVolGroupMan::BreakPathsAboutKinks(std::vector< trex::TTPCOrderedVolGroup >& paths){

  std::vector< trex::TTPCOrderedVolGroup > outPaths;
  // loop over all input paths
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    std::vector< trex::TTPCOrderedVolGroup > tempOutPaths;
    tempOutPaths.emplace_back(fLayout);
    tempOutPaths.emplace_back(fLayout);
    trex::TTPCOrderedVolGroup& outPath1=tempOutPaths[0];
    trex::TTPCOrderedVolGroup& outPath2=tempOutPaths[1];
    trex::TTPCVolGroup kinkGroup(fLayout);
    GetFarHitsPreGroup(path, kinkGroup);
    if(!kinkGroup.size()){
      GetFarHitsGroup(path, kinkGroup);
    }
    TVector3 kinkPosXYZ = kinkGroup.GetAveragePosXYZ();

    std::vector<trex::TTPCPathVolume*>::iterator path1EndIt = path.end();
    bool broken = false;

    if(kinkGroup.size()){
      double minDist2=99999999.;

      // find closest point to kink to break around
      for(std::vector<trex::TTPCPathVolume*>::iterator pointIt = path.begin(); pointIt != path.end(); ++pointIt){
        trex::TTPCPathVolume* point = *pointIt;

        TVector3 diff = point->GetPosXYZ() - kinkPosXYZ;
        double diffDist2 = diff.Mag2();
        if(diffDist2 < minDist2){
          minDist2 = diffDist2;
          path1EndIt = pointIt;
          broken = true;
        }
      }
    }
    // build copy of old hits so the code doesn't get confused when separating them
    trex::TTPCVolGroup newExtHits1(fLayout);
    trex::TTPCVolGroup newExtHits2(fLayout);

    newExtHits1.AddHits(path.GetExtendedHits());
    newExtHits2.AddHits(path.GetExtendedHits());

    outPath1.AddExtendedHits(newExtHits1);
    outPath1.AddBackHits(path.GetBackHits());
    outPath1.SetBackIsVertex(path.GetBackIsVertex());
    // if broken, tie up this path and start another
    if(broken){
      outPath1.AddFrontHits(kinkGroup);
      outPath1.SetFrontIsVertex(true);

      outPath2.AddExtendedHits(newExtHits2);
      outPath2.AddBackHits(kinkGroup);
      outPath2.SetBackIsVertex(true);

      outPath2.AddFrontHits(path.GetFrontHits());
      outPath2.SetFrontIsVertex(path.GetFrontIsVertex());
    }
    else{
      outPath1.AddFrontHits(path.GetFrontHits());
      outPath1.SetFrontIsVertex(path.GetFrontIsVertex());
    };

    // add hits to first path up to its end
    for(std::vector<trex::TTPCPathVolume*>::iterator pointIt = path.begin(); pointIt != path1EndIt; ++pointIt){
      trex::TTPCPathVolume* point = *pointIt;
      outPath1.AddCell(point->GetUnitVolume());
    };

    // add hits to first/same path up to its end
    if(path1EndIt != path.end())
    for(std::vector<trex::TTPCPathVolume*>::iterator pointIt = path1EndIt+1; pointIt != path.end(); ++pointIt){
      trex::TTPCPathVolume* point = *pointIt;
      outPath2.AddCell(point->GetUnitVolume());
    };

    if(outPath2.size()){
      SeparateHits(tempOutPaths);
      outPaths.emplace_back(std::move(tempOutPaths[0]));
      outPaths.emplace_back(std::move(tempOutPaths[1]));
    }
    else{
      outPaths.emplace_back(std::move(tempOutPaths[0]));
    }
  }

  paths=std::move(outPaths);
}

bool trex::TTPCVolGroupMan::GetPathVolOverlap(trex::TTPCOrderedVolGroup& path, trex::TTPCUnitVolume* vol, trex::TTPCConnection::Type type){
  if(!path.size()) return false;

  // don't check when too close to the first or last point in the path
  trex::TTPCPathVolume* firstPoint = *(path.begin());
  trex::TTPCPathVolume* lastPoint = *(path.end()-1);
  trex::TTPCUnitVolume* firstVol = firstPoint->GetUnitVolume();
  trex::TTPCUnitVolume* lastVol = lastPoint->GetUnitVolume();
  for(std::vector<trex::TTPCPathVolume*>::iterator pnt = path.begin(); pnt != path.end(); ++pnt){
    trex::TTPCUnitVolume* curVol = (*pnt)->GetUnitVolume();

    // make sure you're not too close to the first point
    if(HitTest(curVol, firstVol, type)) continue;
    // make sure you're not too close to the last point
    if(HitTest(curVol, lastVol, type)) continue;
    // then check to see how close the volume is
    if(HitTest(curVol, vol, type)) return true;
  };
  return false;
}

bool trex::TTPCVolGroupMan::GetGroupGroupOverlap(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, trex::TTPCConnection::Type type, bool simple, bool checkLean){
  int distX;
  int distY;
  int distZ;
  fLayout->GetTypeDistances(distX, distY, distZ, type);

  if(simple){
    // avoid complicated checks

    trex::TTPCCellInfo3D grp1Pos = group1.GetAveragePad();
    trex::TTPCCellInfo3D grp2Pos = group2.GetAveragePad();

    float nx = (float)( grp2Pos.x - grp1Pos.x )/(float)distX;
    float ny = (float)( grp2Pos.y - grp1Pos.y )/(float)distY;
    float nz = (float)( grp2Pos.z - grp1Pos.z )/(float)distZ;

    return ( (nx*nx) + (ny*ny) + (nz*nz) <= 1.);
  };

  if(checkLean){
    // don't check in directions of opposing lean
    if ((group1.GetXLean() * group2.GetXLean()) < 0) distX = 0;
    if ((group1.GetYLean() * group2.GetYLean()) < 0) distY = 0;
    if ((group1.GetZLean() * group2.GetZLean()) < 0) distZ = 0;
  };

  // check if groups are in range of each other at all
  if( (group1.GetXMin()-group2.GetXMax())>distX || (group2.GetXMin()-group1.GetXMax())>distX ) return false;
  if( (group1.GetYMin()-group2.GetYMax())>distY || (group2.GetYMin()-group1.GetYMax())>distY ) return false;
  if( (group1.GetZMin()-group2.GetZMax())>distZ || (group2.GetZMin()-group1.GetZMax())>distZ ) return false;

  // check if any individual pairs of hits are in range
  for(std::map<long, trex::TTPCUnitVolume*>::iterator hit1El = group1.begin(); hit1El != group1.end(); ++hit1El)
  for(std::map<long, trex::TTPCUnitVolume*>::iterator hit2El = group2.begin(); hit2El != group2.end(); ++hit2El){
    trex::TTPCUnitVolume* vol1 = hit1El->second;
    trex::TTPCUnitVolume* vol2 = hit2El->second;

    if(std::abs(vol1->GetX() - vol2->GetX()) > distX) continue;
    if(std::abs(vol1->GetY() - vol2->GetY()) > distY) continue;
    if(std::abs(vol1->GetZ() - vol2->GetZ()) > distZ) continue;
    return true;
  };

  return false;
}

void trex::TTPCVolGroupMan::MergeGroups(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, trex::TTPCVolGroup& merged){

  merged.AddHits(group1);
  merged.AddHits(group2);

}

void trex::TTPCVolGroupMan::BulkGroup(trex::TTPCVolGroup& group){
  int xMin = group.GetXMin();
  int xMax = group.GetXMax();
  int yMin = group.GetYMin();
  int yMax = group.GetYMax();
  int zMin = group.GetZMin();
  int zMax = group.GetZMax();

  // TODO: either add directly or by lookup; whichever is faster
  for(std::map<long, trex::TTPCUnitVolume*>::iterator hitEl = fPrimaryHits.begin(); hitEl != fPrimaryHits.end(); ++hitEl){
    trex::TTPCUnitVolume* vol = hitEl->second;

    if(vol->GetX() >= xMin && vol->GetX() <= xMax){
      if(vol->GetY() >= yMin && vol->GetY() <= yMax){
        if(vol->GetZ() >= zMin && vol->GetZ() <= zMax){
          group.AddHit(vol);
        }
      }
    }
  }
}

void trex::TTPCVolGroupMan::BulkGroups(std::vector< trex::TTPCVolGroup >& groups){
  for(std::vector< trex::TTPCVolGroup >::iterator groupIt = groups.begin(); groupIt != groups.end(); ++groupIt){
    BulkGroup(*groupIt);
  };
}

void trex::TTPCVolGroupMan::GetUnorderedGroup(trex::TTPCOrderedVolGroup& in, trex::TTPCVolGroup& out){
 
  for(std::vector<trex::TTPCPathVolume*>::iterator it = in.begin(); it != in.end(); ++it){
    out.AddHit((*it)->GetUnitVolume());
  };
}

void trex::TTPCVolGroupMan::BuildGroupFriends(trex::TTPCOrderedVolGroup& in, trex::TTPCConnection::Type type){
  // max distance from path
  int typeX;
  int typeY;
  int typeZ;
  fLayout->GetTypeDistances(typeX, typeY, typeZ, type);

  // extended group to create
  trex::TTPCVolGroup group(fLayout);

  // add edge hits
  if(in.HasFrontHits()) group.AddHits(in.frontHitsBegin(), in.frontHitsEnd());
  if(in.HasBackHits()) group.AddHits(in.backHitsBegin(), in.backHitsEnd());

  // two different modes for searching path hits - use whichever involves fewest iterations
  int volIterations = (4./3.)*3.14 * typeX*typeY*typeZ;
  if(volIterations > fPrimaryHits.size()){
    for(std::map<long, trex::TTPCUnitVolume*>::iterator hitEl = fPrimaryHits.begin(); hitEl != fPrimaryHits.end(); ++hitEl){
      bool nearEnough = false;
      for(std::vector<trex::TTPCPathVolume*>::iterator inputIt = in.begin(); inputIt != in.end(); inputIt++){
        trex::TTPCUnitVolume* inputVol = (*inputIt)->GetUnitVolume(); 

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
      if(nearEnough) group.AddHit(hitEl->second);
    };
  }
  else{
    for(std::vector<trex::TTPCPathVolume*>::iterator inputIt = in.begin(); inputIt != in.end(); inputIt++){
      trex::TTPCUnitVolume* inputVol = (*inputIt)->GetUnitVolume(); 

      // hits near cell
      trex::TTPCVolGroup nearHits(fLayout);
      GetNearHits(fPrimaryHits,nearHits, inputVol->GetID(), type, true);
      // add all of them
      group.AddHits(nearHits);
    };
  };

  // add associated hits to the group
  in.AddExtendedHits(group);
}

void trex::TTPCVolGroupMan::ClusterGroupFriends(trex::TTPCOrderedVolGroup& in, bool doClustering, bool checkX, bool partial){
  // note - path may be empty after this
  // check if group has extended hits and produce them if not
  if(!in.GetHasExtendedHits()){
    BuildGroupFriends(in);
  };
  // stop here if 'do clustering' is off
  if(!doClustering) return;
  // if needed, check whether path is an x-path before going further
  if(checkX){
    if(IsXPathCandidate(in)){
      in.SetIsXPath(true);
    };
  };

  // connect all assocaited hits
  in.DoClustering(partial);
}

void trex::TTPCVolGroupMan::BuildAllFriends(std::vector< trex::TTPCOrderedVolGroup >& paths){
  // build friends for paths
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    BuildGroupFriends(*pathIt);
  }
  // ensure that none share hits
  SeparateHits(paths);
}

void trex::TTPCVolGroupMan::SeparateHits(std::vector< trex::TTPCOrderedVolGroup >& paths){
  // loop over all pairs of groups
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It)
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
    trex::TTPCOrderedVolGroup& path1 = *path1It;
    trex::TTPCOrderedVolGroup& path2 = *path2It;

    // take associated hits
    trex::TTPCVolGroup& extendedHits1=path1.GetExtendedHits();
    trex::TTPCVolGroup& extendedHits2=path2.GetExtendedHits();

    // loop over path associated hits
    for(std::map<long, trex::TTPCUnitVolume*>::iterator extHit1El = extendedHits1.begin(); extHit1El != extendedHits1.end(); ++extHit1El)
    for(std::map<long, trex::TTPCUnitVolume*>::iterator extHit2El = extendedHits2.begin(); extHit2El != extendedHits2.end(); ++extHit2El){
      // shared hit - break tie in favour of closest path
      if(extHit1El->first == extHit2El->first){
        // ignore if already marked for clearing
        if(!extHit1El->second) continue;

        double path1Dist = GetMinDistance(path1, extHit1El->second);
        double path2Dist = GetMinDistance(path2, extHit1El->second);

        if(path1Dist > path2Dist){
          // mark hit for removal from from path1
          extendedHits1.MarkHit(extHit1El->first);
        }
        else{
          // mark hit for removal from path2
          extendedHits2.MarkHit(extHit1El->first);
        };
      };
    };
    // clear hits marked for deletion
    extendedHits1.ClearMarked();
    extendedHits2.ClearMarked();

    // renew path associated hits
    path1.AddExtendedHits(extendedHits1);
    path2.AddExtendedHits(extendedHits2);
  };
}

void trex::TTPCVolGroupMan::SeparateXPathHits(std::vector< trex::TTPCOrderedVolGroup >& paths){
  // loop over all paths
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;

    if(path.GetIsXPath()){
      // get all hits in path
      trex::TTPCVolGroup pathHits(fLayout);
      for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
        trex::TTPCPathVolume* pathVol = *pathVolIt;
        for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
          trex::TTPCUnitVolume* vol = *volIt;
          pathHits.AddHit(vol);
        };
      };

      // add those hits to junctions either side; later code will take care of any overlap
      if(path.GetFrontIsVertex()){
        path.GetFrontHits().AddHits(pathHits);
      };
      if(path.GetBackIsVertex()){
        path.GetBackHits().AddHits(pathHits);
      };
    };
  };

  // now remove all x-paths
  std::vector<trex::TTPCOrderedVolGroup> pathsToKeep;

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator reaper = paths.begin(); reaper != paths.end(); ++reaper){
    if(!reaper->GetIsXPath()){
      pathsToKeep.emplace_back(std::move(*reaper));
    }
  }
  paths=std::move(pathsToKeep);

}

void trex::TTPCVolGroupMan::SeparateEmptyClusters(std::vector< trex::TTPCOrderedVolGroup >& paths){
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    path.Clean();
  };
}
void trex::TTPCVolGroupMan::SeparateClusterHits(std::vector< trex::TTPCOrderedVolGroup >& paths){
  int mergeX;
  int mergeY;
  int mergeZ;
  fLayout->GetTypeDistances(mergeX, mergeY, mergeZ, trex::TTPCConnection::clusterMerge);

  // get rid of any clusters whos seed is shared with other paths
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It)
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
    trex::TTPCOrderedVolGroup& path1 = *path1It;
    trex::TTPCOrderedVolGroup& path2 = *path2It;

    for(std::vector<trex::TTPCPathVolume*>::iterator cluster1It = path1.begin(); cluster1It != path1.end(); ++cluster1It)
    for(std::vector<trex::TTPCPathVolume*>::iterator cluster2It = path2.begin(); cluster2It != path2.end(); ++cluster2It){
      trex::TTPCPathVolume* cluster1 = *cluster1It;
      trex::TTPCPathVolume* cluster2 = *cluster2It;

      if(!cluster1) continue;
      if(!cluster2) continue;

      trex::TTPCUnitVolume* vol1 = cluster1->GetUnitVolume();
      trex::TTPCUnitVolume* vol2 = cluster2->GetUnitVolume();

      if(vol1 == vol2){
        cluster1->ClearFriends();
        cluster2->ClearFriends();
      };
    };
  };

  // loop over pairs of paths, checking for overlapping clusters (ignore x-paths for now)
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator path1It = paths.begin(); path1It != paths.end(); ++path1It){
    trex::TTPCOrderedVolGroup& path1 = *path1It;
    if(path1.GetIsXPath()) continue;

    for(std::vector< trex::TTPCOrderedVolGroup >::iterator path2It = path1It+1; path2It != paths.end(); ++path2It){
      trex::TTPCOrderedVolGroup& path2 = *path2It;
      if(path2.GetIsXPath()) continue;

      for(std::vector<trex::TTPCPathVolume*>::iterator cluster1It = path1.begin(); cluster1It != path1.end(); ++cluster1It)
      for(std::vector<trex::TTPCPathVolume*>::iterator cluster2It = path2.begin(); cluster2It != path2.end(); ++cluster2It){
        trex::TTPCPathVolume* cluster1 = *cluster1It;
        trex::TTPCPathVolume* cluster2 = *cluster2It;
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
        for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = cluster1->GetFriendsBegin(); volIt != cluster2->GetFriendsEnd(); ++volIt){
          if(volIt == cluster1->GetFriendsEnd()) volIt = cluster2->GetFriendsBegin();
          trex::TTPCUnitVolume* vol = *volIt;

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
          for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = cluster1->GetFriendsBegin(); volIt != cluster1->GetFriendsEnd(); ++volIt){
            trex::TTPCUnitVolume* vol = *volIt;

            if( normalHierarchy^((verticalMerge && (vol->GetY()>gapPos)) || (!verticalMerge && (vol->GetZ()>gapPos))) ){
              cluster2->AddFriend(vol);
              cluster1->MarkFriend(volIt);
            };
          };
          cluster1->ClearMarked();
          for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = cluster2->GetFriendsBegin(); volIt != cluster2->GetFriendsEnd(); ++volIt){
            trex::TTPCUnitVolume* vol = *volIt;

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
      path1.Clean();
      path2.Clean();
    };
  };
};

void trex::TTPCVolGroupMan::SeparateAnomHits(std::vector< trex::TTPCOrderedVolGroup >& paths){
  // ready to add hits to junctions
  std::vector< trex::TTPCVolGroup* > vertices = GetJunctionsFromPaths(paths);

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;

    if(path.GetFrontIsVertex() || path.GetBackIsVertex()){
      if(path.GetFrontIsVertex()){
        unsigned int junctionID = path.GetFrontID();
        std::vector<trex::TTPCUnitVolume*> hitsToAdd = SeparateAnomHitsPath(path, -1);

        // now add discarded hits to junction
        for(std::vector< trex::TTPCVolGroup* >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
          trex::TTPCVolGroup& vertex = **vertexIt;
          if(junctionID == vertex.GetID()){
            vertex.AddHits(hitsToAdd);
          };
        };
      };
      if(path.GetBackIsVertex()){
        unsigned int junctionID = path.GetBackID();
        std::vector<trex::TTPCUnitVolume*> hitsToAdd = SeparateAnomHitsPath(path, 1);

        // now add discarded hits to junction
        for(std::vector< trex::TTPCVolGroup* >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
          trex::TTPCVolGroup& vertex = **vertexIt;
          if(junctionID == vertex.GetID()){
            vertex.AddHits(hitsToAdd);
          };
        };
      };
    };
  };
}

std::vector<trex::TTPCUnitVolume*> trex::TTPCVolGroupMan::SeparateAnomHitsPath(trex::TTPCOrderedVolGroup& path, int checkDir){
  int checkDist = fLayout->GetAnomCheckDist();
  int projectDist = fLayout->GetAnomProjectDist();
  double maxOffs = fLayout->GetAnomMaxOffs();

  // be prepared to add hits to junction
  std::vector<trex::TTPCUnitVolume*> hitsToAdd;

  // return empty handed if path is too small
  if(path.size() <= (checkDist + projectDist)){
    return std::move(hitsToAdd);
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
    rangeZero = path.size()-1;
    rangeStart = rangeZero - checkDist;
    rangeEnd = rangeStart - projectDist;
  }
  else{
    return std::move(hitsToAdd);
  };

  // look for hits with max time deviation from expectation, approximated by extrapolating a straight line between the start and end of test range
  trex::TTPCPathVolume* range1 = path.at(rangeStart);
  trex::TTPCPathVolume* range2 = path.at(rangeEnd);

  TVector3 startPos = GetAvgPosRep(range1);
  TVector3 endPos = GetAvgPosRep(range2);
  TVector3 parNorm = (endPos - startPos).Unit();
  bool xOnly = false;
  double maxDev = -1.;

  // find max in the control region
  for(int i=rangeStart; i != rangeEnd; i+= checkDir){
    trex::TTPCPathVolume* pathVol = path.at(i);
    for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
      trex::TTPCUnitVolume* vol = *volIt;
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
    trex::TTPCPathVolume* pathVol = path.at(i);
    bool removed = false;
    for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); ++volIt){
      trex::TTPCUnitVolume* vol = *volIt;
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
  path.Clean();

  return std::move(hitsToAdd);
}

void trex::TTPCVolGroupMan::SeparateJunctionHits(std::vector< trex::TTPCOrderedVolGroup >& paths){
  int mergeX;
  int mergeY;
  int mergeZ;
  fLayout->GetTypeDistances(mergeX, mergeY, mergeZ, trex::TTPCConnection::clusterMerge);

  // first build list of potential vertices (path groups tagged as vertices)
  std::vector< trex::TTPCVolGroup* > tempVertices = GetJunctionsFromPaths(paths);

  // select un-duplicated vertices and merge together any that overlap
  std::vector< trex::TTPCVolGroup* > vertices;
  for(std::vector< trex::TTPCVolGroup* >::iterator tempVertexIt = tempVertices.begin(); tempVertexIt != tempVertices.end(); tempVertexIt++){
    trex::TTPCVolGroup& tempVertex = **tempVertexIt;

    // check if each vertex overlaps with a current one and merge them if it does
    bool vertexFound = false;
    for(std::vector< trex::TTPCVolGroup* >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
      trex::TTPCVolGroup& vertex = **vertexIt;

      // as a timesaver, ensure vertices might overlap before doing more intensive checks
      bool potentialOverlap = true;
      if(vertex.GetXMin() > tempVertex.GetXMax() || tempVertex.GetXMin() > vertex.GetXMax()) potentialOverlap = false;
      if(vertex.GetYMin() > tempVertex.GetYMax() || tempVertex.GetYMin() > vertex.GetYMax()) potentialOverlap = false;
      if(vertex.GetZMin() > tempVertex.GetZMax() || tempVertex.GetZMin() > vertex.GetZMax()) potentialOverlap = false;

      bool thisVertexFound = false;
      if(potentialOverlap){
        // overlaps between the two
        for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = tempVertex.begin(); volEl != tempVertex.end(); ++volEl){
          if(vertex.Contains(volEl->first)){
            vertexFound = thisVertexFound = true;
            break;
          };
        };
      };
      if(thisVertexFound){
        // merge
        vertex.AddHits(tempVertex);
      };
    };
    if(!vertexFound){
      // add vertex if it doesn't already exist
      vertices.push_back(&tempVertex);
    };
  };

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    // check if any HV cluster overlaps with any vertex and swap those that do from the cluster to the vertex
    std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin();
    while(pathVolIt != path.end()){
      if(!*pathVolIt) path.erase(pathVolIt);
      else pathVolIt++;
    };
  };

  // now merge overlapping HV clusters from paths into vertices
  bool mightMerge = true;
  while(mightMerge){
    mightMerge = false;
    for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      trex::TTPCOrderedVolGroup& path = *pathIt;

      // build sets of hits in horizontal clusters and vertical clusters
      std::set<trex::TTPCUnitVolume*> clusteredH;
      std::set<trex::TTPCUnitVolume*> clusteredV;
      if(!path.GetIsXPath()){
        // add hits in horizontal and vertical clusters to their respective sets
        for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
          trex::TTPCPathVolume* pathVol = *pathVolIt;
          if(!pathVol) continue;

          for(std::vector<trex::TTPCUnitVolume*>::iterator pathFriendIt = pathVol->GetFriendsBegin(); pathFriendIt != pathVol->GetFriendsEnd(); ++pathFriendIt){
            trex::TTPCUnitVolume* pathFriend = *pathFriendIt;

            if(pathVol->GetIsVertical()){
              clusteredV.insert(pathFriend);
            }
            else{
              clusteredH.insert(pathFriend);
            };
          };
        };
      };

      for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
        trex::TTPCPathVolume* pathVol = *pathVolIt;
        if(!pathVol) continue;

        bool verticalMerge = pathVol->GetIsVertical();

        for(std::vector< trex::TTPCVolGroup* >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
          trex::TTPCVolGroup& vertex = **vertexIt;

          int xMin = vertex.GetXMin();
          int xMax = vertex.GetXMax();
          int yMin = vertex.GetYMin();
          int yMax = vertex.GetYMax();
          int zMin = vertex.GetZMin();
          int zMax = vertex.GetZMax();

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
          for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathVol->GetFriendsBegin(); volIt != pathVol->GetFriendsEnd(); volIt++){
            trex::TTPCUnitVolume* vol = *volIt;
            if(!vol) continue;

            // first check if vol might actually overlap to save time
            if(vol->GetX() < xMin || vol->GetX() > xMax) continue;
            if(vol->GetY() < yMin || vol->GetY() > yMax) continue;
            if(vol->GetZ() < zMin || vol->GetZ() > zMax) continue;

            // potential hit - do a more intense search
            for(std::map<long, trex::TTPCUnitVolume*>::iterator vertexEl = vertex.begin(); vertexEl != vertex.end(); ++vertexEl){
              trex::TTPCUnitVolume* vertexVol = vertexEl->second;

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
                  vertex.AddHit(vol);
		  std::cout<<"Vertical merge!"<<std::endl;
                };
                clusteredV.erase(vol);
              }
              else{
                if(!clusteredV.count(vol)){
                  vertex.AddHit(vol);
		  std::cout<<"Horizontal merge!"<<std::endl;
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
	    std::cout<<"Cluster goes from "<<pathVol->GetClusterSize();
            pathVol->ClearMarked();
	    std::cout<<" to "<<pathVol->GetClusterSize()<<" hits"<<std::endl;

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
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;

    if(path.GetIsXPath()) continue;
    if(path.size() < 2) continue;

    // look for gap at the back
    if(path.GetBackIsVertex()){
      bool possibleIso = true;
      while(possibleIso){
        possibleIso = false;
        trex::TTPCPathVolume* beginVol = 0;
        std::vector<trex::TTPCPathVolume*>::iterator beginVolIt;
        for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
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
          for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = beginVolIt+1; pathVolIt != path.end(); ++pathVolIt){
            trex::TTPCPathVolume* pathVol = *pathVolIt;
            if(!pathVol) continue;

            if(wasVertical==pathVol->GetIsVertical()){
              int pos = pathVol->GetIsVertical() ? pathVol->GetZ() : pathVol->GetY();

              // check if there's a gap and add everything before it to the vertex if so
              if(std::abs(pos - lastPos) > 1){
                for(std::vector<trex::TTPCPathVolume*>::iterator gapVolIt = path.begin(); gapVolIt != pathVolIt; ++gapVolIt){
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
    if(path.GetFrontIsVertex()){
      bool possibleIso = true;
      while(possibleIso){
        possibleIso = false;
        std::vector<trex::TTPCPathVolume*>::reverse_iterator endVolRIt;
        trex::TTPCPathVolume* endVol = 0;
        for(std::vector<trex::TTPCPathVolume*>::reverse_iterator pathVolRIt = path.rbegin(); pathVolRIt != path.rend(); ++pathVolRIt){
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
          for(std::vector<trex::TTPCPathVolume*>::reverse_iterator pathVolRIt = endVolRIt; pathVolRIt != path.rend(); ++pathVolRIt){
            trex::TTPCPathVolume* pathVol = *pathVolRIt;
            if(!pathVol) continue;

            if(wasVertical==pathVol->GetIsVertical()){
              int pos = pathVol->GetIsVertical() ? pathVol->GetZ() : pathVol->GetY();

              // check if there's a gap and add everything before it to the vertex if so
              if(std::abs(pos - lastPos) > 1){
                for(std::vector<trex::TTPCPathVolume*>::reverse_iterator gapVolRIt = path.rbegin(); gapVolRIt != pathVolRIt; ++gapVolRIt){
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
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    // check if any HV cluster overlaps with any vertex and swap those that do from the cluster to the vertex
    std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin();
    while(pathVolIt != path.end()){
      if(!*pathVolIt) path.erase(pathVolIt);
      else pathVolIt++;
    };
  };

  // copy expanded vertices to their initial positions
  for(std::vector< trex::TTPCVolGroup* >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
    trex::TTPCVolGroup& vertex = **vertexIt;

    // replace any overlapping path beginning or end with this
    for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      trex::TTPCOrderedVolGroup& path = *pathIt;

      // overlap with front hits
      bool frontOverlaps = path.GetFrontID() == vertex.GetID();
      bool backOverlaps = path.GetBackID() == vertex.GetID();

      // set as new vector
      if(frontOverlaps) path.AddFrontHits(vertex);
      if(backOverlaps) path.AddBackHits(vertex);
    };
  };
}
void trex::TTPCVolGroupMan::EnforceOrdering(std::vector< trex::TTPCOrderedVolGroup >& paths){
  // loop over all paths
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt !=paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    // count turning points
    int turningPoints90=0;
    bool wasVertical=true;
    for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = path.begin(); pathVolIt != path.end(); ++pathVolIt){
      trex::TTPCPathVolume* pathVol = *pathVolIt;

      bool isVertical = pathVol->GetIsVertical();
      if(isVertical ^ wasVertical){
        if(pathVolIt != path.begin()) turningPoints90++;
      };
      wasVertical = isVertical;
    };
    // define turning points as half changes from horizontal to vertical or vice versa
    int turningPoints = turningPoints90/2;

    // case #1: positive z
    if(turningPoints<1){
      path.OrderForwardsDirection();
    }
    // case #2: negative curvature
    else if(turningPoints<2){
      path.OrderNegativeCurvature();
    }
    // case #3:
    else{
      //TODO: implement as starting from junction?
      path.OrderForwardsDirection();
    }
  };
}
bool trex::TTPCVolGroupMan::GetOverlaps(trex::TTPCVolGroup& group1, trex::TTPCVolGroup& group2, bool checkHits){
  // basic check
  if(group1.GetXMin() > group2.GetXMax() || group2.GetXMin() > group1.GetXMax()) return false;
  if(group1.GetYMin() > group2.GetYMax() || group2.GetYMin() > group1.GetYMax()) return false;
  if(group1.GetZMin() > group2.GetZMax() || group2.GetZMin() > group1.GetZMax()) return false;

  // if not flagged to check hits, just return here
  if(!checkHits) return true;

  // compare for individual hits for overlaps
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el1 = group1.begin(); el1 != group1.end(); ++el1)
  for(std::map<long, trex::TTPCUnitVolume*>::iterator el2 = group2.begin(); el2 != group2.end(); ++el2)
  if(el1->first == el2->first) return true;

  return false;
}

void trex::TTPCVolGroupMan::ResetVertexStatuses(std::vector< trex::TTPCOrderedVolGroup >& paths, bool partial){
  // first build list of potential vertices (path groups tagged as vertices)
  std::vector< trex::TTPCVolGroup* > vertices = GetJunctionsFromPaths(paths);

  for(std::vector< trex::TTPCVolGroup* >::iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt){
    trex::TTPCVolGroup& vertex = **vertexIt;
    int connectedPaths = 0;

    // count connected paths
    for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      trex::TTPCOrderedVolGroup& path = *pathIt;

      // overlap with front hits
      bool frontOverlaps = path.GetFrontID() == vertex.GetID();
      bool backOverlaps = path.GetBackID() == vertex.GetID();

      // increment connected paths
      if(frontOverlaps) connectedPaths++;
      if(backOverlaps) connectedPaths++;
    };

    // status as vertex
    bool isVertex = connectedPaths>1;

    // replace any overlapping path beginning or end with this
    for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      trex::TTPCOrderedVolGroup& path = *pathIt;

      // overlap with front hits
      bool frontOverlaps = path.GetFrontID() == vertex.GetID();
      bool backOverlaps = path.GetBackID() == vertex.GetID();


      // set new status as vertex, adding hits if the status has changed
      bool redoClustering = false;
      if(frontOverlaps){
        if(path.GetFrontIsVertex() != isVertex){
          redoClustering = true;
        };
        path.SetFrontIsVertex(isVertex);
      }
      else if(backOverlaps){
        if(path.GetBackIsVertex() != isVertex){
          redoClustering = true;
        };
        path.SetBackIsVertex(isVertex);
      };

      // do NOT allow hits the same vertex to be duplicated
      if(path.GetFrontID() == path.GetBackID()){
        if(path.GetFrontIsVertex() && path.GetBackIsVertex()){
          path.SetBackIsVertex(false);
        };
      }
      else if(redoClustering){
        path.GetExtendedHits().AddHits(vertex);
        path.DoClustering(partial);
      }
    }
  }
}

void trex::TTPCVolGroupMan::GetConnectedHits(std::vector<trex::TTPCVolGroup>& out, trex::TTPCConnection::Type type, trex::TTPCHitGroupings::Type typeFilter, bool usabilityCheck){
  GetConnectedHits(fPrimaryHits, out, type, typeFilter, usabilityCheck);
}

void trex::TTPCVolGroupMan::GetConnectedHits(trex::TTPCVolGroup& in, std::vector<trex::TTPCVolGroup>& out, trex::TTPCConnection::Type type, trex::TTPCHitGroupings::Type typeFilter, bool usabilityCheck){

  // make dummy group of all cells
  trex::TTPCVolGroup conHits(fLayout);
    conHits.AddHitMap(in.GetHitMap());

  for(int i=0; i<9999; i++){
    // break out of loop if dummy becomes empty
    if(conHits.size() < 1) break;
    // start at first element in dummy
    std::map<long, trex::TTPCUnitVolume*>::iterator el = conHits.begin();

    // define new group
    out.emplace_back(fLayout);
    // recursively build set of connected cells, which are removed from the dummy and added to the new group
    RecursiveFriendBuild(el->first, out.back(), conHits, type, typeFilter);

    // check to make sure this group is useful
    if(usabilityCheck){
      if(!CheckUsability(out.back())){
        out.pop_back();
      }
    }
  }
}

bool trex::TTPCVolGroupMan::CheckUsability(trex::TTPCVolGroup& inGroup){
  unsigned int minSize = fLayout->GetMinPatternPads();

  std::set< std::pair<int, int> > patternPads;
  for(std::map<long, trex::TTPCUnitVolume*>::iterator cellEl = inGroup.begin(); cellEl != inGroup.end(); ++cellEl){
    std::pair<int, int> patternPad (cellEl->second->GetY(), cellEl->second->GetZ());
    patternPads.insert(patternPad);
  };

  // count number of pads
  return (minSize <= patternPads.size());
}

void trex::TTPCVolGroupMan::GetNearHits(trex::TTPCVolGroup& in, trex::TTPCVolGroup& cellsOut, long id, trex::TTPCConnection::Type type, bool inclusive, bool singular, float distFilter, trex::TTPCHitGroupings::Type typeFilter, bool square){
  // try to find id in map
  std::map<long, trex::TTPCUnitVolume*>::iterator el = in.find(id);
  // return empty group if it's not found
  if(el == in.end()) return;
  // otherwise get hits near to its specified volume
  trex::TTPCUnitVolume* vol = el->second;

  GetNearHits(in, cellsOut,vol, type, inclusive, singular, distFilter, typeFilter, square);
}
void trex::TTPCVolGroupMan::GetNearHits(trex::TTPCVolGroup& in, trex::TTPCVolGroup& cellsOut, trex::TTPCUnitVolume* vol, trex::TTPCConnection::Type type, bool inclusive, bool singular, float distFilter, trex::TTPCHitGroupings::Type typeFilter, bool square){
  // set up empty group of nearby hits

  // fetch x, y and z ids and edge statuses from provided volume
  int x = vol->GetX();
  int y = vol->GetY();
  int z = vol->GetZ();

  int distXP;
  int distYP;
  int distZP;
  int distXN;
  int distYN;
  int distZN;
  fLayout->GetTypeDistances(distXP, distYP, distZP, type);
  fLayout->GetTypeDistances(distXN, distYN, distZN, type);


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
        trex::TTPCUnitVolume* newHit = in.GetHit(newID);
        if (!newHit) continue;

        if (distFilter>-1. && distFilter>=newHit->GetFriendDist()) continue;
        if (distFilter>-1. && newHit->GetFriendDist() < 0) continue;

        // if id is contained in hit map, add the corresponding cell to the list of near hits
        bool added = cellsOut.AddHit(newHit);
        if(singular && added){
          return;
        };
      };
    };
  };
  return;
}

std::vector< trex::TTPCVolGroup* > trex::TTPCVolGroupMan::GetJunctionsFromPaths(std::vector< trex::TTPCOrderedVolGroup >& paths){
  // count vertices in input

  std::vector< trex::TTPCVolGroup* > tempVertices;
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    // find front junction
    if(path.GetFrontIsVertex()) tempVertices.push_back(&(path.GetFrontHits()));
    // find back junction
    if(path.GetBackIsVertex()) tempVertices.push_back(&(path.GetBackHits()));
  };

  // add unique vertices to output
  std::vector< trex::TTPCVolGroup* > vertices;
  for(std::vector< trex::TTPCVolGroup* >::iterator vertex1It = tempVertices.begin(); vertex1It != tempVertices.end(); ++vertex1It){
    trex::TTPCVolGroup* vertex1 = *vertex1It;
    bool found = false;

    for(std::vector< trex::TTPCVolGroup* >::iterator vertex2It = vertices.begin(); vertex2It != vertices.end(); ++vertex2It){
      trex::TTPCVolGroup* vertex2 = *vertex2It;
      if(vertex1->GetID() == vertex2->GetID()){
        found = true;
        break;
      };
    };
    if(!found){
      vertices.push_back(vertex1);
    };
  };
  return std::move(vertices);
}

void trex::TTPCVolGroupMan::GetUnusedHits(std::vector< trex::TTPCOrderedVolGroup >& paths, TTPCVolGroup& unusedHits){
  
  std::vector< trex::TTPCVolGroup* > junctions = GetJunctionsFromPaths(paths);

  // start full

  unusedHits.AddHits(fPrimaryHits);

  // remove path hits
  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;
    for(std::vector<trex::TTPCPathVolume*>::iterator clusterIt = path.begin(); clusterIt != path.end(); ++clusterIt){
      trex::TTPCPathVolume* cluster = *clusterIt;
      for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = cluster->GetFriendsBegin(); volIt != cluster->GetFriendsEnd(); ++volIt){
        trex::TTPCUnitVolume* vol = *volIt;

        unusedHits.RemoveHit(vol->GetID());
      }
    }
  }
  // remove junction hits
  for(std::vector< trex::TTPCVolGroup* >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    trex::TTPCVolGroup& junction = **junctionIt;
    for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = junction.begin(); volEl != junction.end(); ++volEl){
      unusedHits.RemoveHit(volEl->first);
    }
  }
}


//MDH TODO: This is building new junctions in a stack container then passing out their pointers!!! :-S
void trex::TTPCVolGroupMan::AssociateUnusedWithJunctions(trex::TTPCVolGroup& unused, 
							 std::vector< trex::TTPCOrderedVolGroup >& paths,
							 std::vector< trex::TTPCVolGroup >& junctionsToAddTo){

  std::vector< trex::TTPCVolGroup* > junctions = GetJunctionsFromPaths(paths);
  if(!junctions.size()) return;

  // add all junction hits to unused
  for(std::vector< trex::TTPCVolGroup* >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    trex::TTPCVolGroup& junction = **junctionIt;
    if(!junction.size()) continue;

    unused.AddHits(junction);
  };

  std::vector<trex::TTPCVolGroup*> newJuncPtrs;

  // try and associate each junction with new hits
  for(std::vector< trex::TTPCVolGroup* >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
    trex::TTPCVolGroup& junction = **junctionIt;

    // break out of loop if unused becomes empty, otherwise try to add hits starting from anywhere in junction
    if(unused.size() < 1) break;

    // new junction to hold hits
    junctionsToAddTo.emplace_back(fLayout, trex::TTPCVolGroup::GetFreeID());

    // seed from every element in old junction so nothing is missed
    for(std::map<long, trex::TTPCUnitVolume*>::iterator seedEl = junction.begin(); seedEl != junction.end(); ++seedEl){
      RecursiveFriendBuild(seedEl->first, junctionsToAddTo.back(), unused);
    }

    if(!junctionsToAddTo.back().size()){
      junctionsToAddTo.pop_back();
    }
    else{
      newJuncPtrs.push_back(&(junctionsToAddTo.back()));
    }
  }

  // now add expanded junctions to paths
  for(auto newJunctionIt = newJuncPtrs.begin(); newJunctionIt != newJuncPtrs.end(); ++newJunctionIt){
    trex::TTPCVolGroup& newJunction = **newJunctionIt;
    if(!newJunction.size()) continue;

    for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      trex::TTPCOrderedVolGroup& path = *pathIt;

      bool frontOverlaps = GetOverlaps(newJunction, path.GetFrontHits());
      bool backOverlaps = GetOverlaps(newJunction, path.GetBackHits());

      if(frontOverlaps) path.AddFrontHits(newJunction);
      if(backOverlaps) path.AddBackHits(newJunction);
    }
  }
}

void trex::TTPCVolGroupMan::SanityFilter(std::vector< trex::TTPCOrderedVolGroup >& input){
  std::vector< trex::TTPCOrderedVolGroup > output;

  for(std::vector< trex::TTPCOrderedVolGroup >::iterator pathIt = input.begin(); pathIt != input.end(); ++pathIt){
    trex::TTPCOrderedVolGroup& path = *pathIt;

    // ensure that path contains no empty clusters
    std::vector<trex::TTPCPathVolume*>::iterator reaper = path.begin();
    while(reaper != path.end()){
      if(!(*reaper)->GetClusterSize()){
        path.erase(reaper);
      }
      else{
        ++reaper;
      };
    };

    // ensure path contains above threshold of clusters
    if(path.size() < fLayout->GetMinPathClusters()) continue;
    output.emplace_back(std::move(path));
  };

  input=std::move(output);
}

void trex::TTPCVolGroupMan::RecursiveFriendBuild(long startID, trex::TTPCVolGroup& target, trex::TTPCVolGroup& source, trex::TTPCConnection::Type type, trex::TTPCHitGroupings::Type typeFilter){
  // initialise start cell by provided id
  trex::TTPCUnitVolume* startCell = fHitMap[startID];
  // check if cell exists in source
  if(!source.Contains(startID)) return;
  // attempt to remove start cell from the group being drained
  if(!source.RemoveHit(startID)) return;
  // attempt to add start cell to the group being filled
  if(!target.AddHit(startCell)) return;

  // get near hits from the group being drained and attempt to add all of them
  trex::TTPCVolGroup nearHits(fLayout);
  GetNearHits(source,nearHits, startCell, type, false,false,-1., typeFilter);
  for(std::map<long, trex::TTPCUnitVolume*>::iterator vol = nearHits.begin(); vol != nearHits.end(); ++vol){
    // repeat the process on every near hit found
    RecursiveFriendBuild(vol->first, target, source, type, typeFilter);
  };
}

void trex::TTPCVolGroupMan::RecursiveFriendSeek(trex::TTPCOrderedVolGroup& inList, trex::TTPCVolGroup& target, float dist, trex::TTPCConnection::Type type, trex::TTPCHitGroupings::Type typeFilter){
  trex::TTPCVolGroup unordered(fLayout);
  GetUnorderedGroup(inList,unordered);
  RecursiveFriendSeek(unordered, target, dist, type, typeFilter);
}
void trex::TTPCVolGroupMan::RecursiveFriendSeek(trex::TTPCVolGroup& inList, trex::TTPCVolGroup& target, float dist, trex::TTPCConnection::Type type, trex::TTPCHitGroupings::Type typeFilter){
  for(std::map<long, trex::TTPCUnitVolume*>::iterator pnt = fPrimaryHits.begin(); pnt != fPrimaryHits.end(); ++pnt){
    pnt->second->SetFriendDist(9999);
  };
  for(std::map<long, trex::TTPCUnitVolume*>::iterator pnt = inList.begin(); pnt != inList.end(); ++pnt){
    pnt->second->SetFriendDist(0);
    RecursiveFriendListSeek(pnt->first, fPrimaryHits, 0, dist, type, typeFilter);
  };

  for(std::map<long, trex::TTPCUnitVolume*>::iterator pnt = fPrimaryHits.begin(); pnt != fPrimaryHits.end(); ++pnt){
    if(pnt->second->GetFriendDist() < 9998){
      target.AddHit(pnt->second);
    };
  };
}

void trex::TTPCVolGroupMan::RecursiveFriendListSeek(long startID, trex::TTPCVolGroup& source, float curDist, float dist, trex::TTPCConnection::Type type, trex::TTPCHitGroupings::Type typeFilter){
  // break if over or at max distance
  if (curDist > dist) return;

  // set current cell's distance based on current input distance
  trex::TTPCUnitVolume* startCell = fHitMap[startID];
  startCell->SetFriendDist(curDist);

  // get near hits from the group being drained and attempt to add all of them
  trex::TTPCVolGroup nearHits(fLayout);
  GetNearHits(source, nearHits, startCell, type, false, false, dist, typeFilter);
  for(std::map<long, trex::TTPCUnitVolume*>::iterator vol = nearHits.begin(); vol != nearHits.end(); ++vol){
    // repeat the process on every near hit found
    RecursiveFriendListSeek(vol->first, source, curDist+1., dist, type, typeFilter);
  };
}

TVector3 trex::TTPCVolGroupMan::GetAvgPosRep(trex::TTPCPathVolume* vol){
  TVector3 avgPos = vol->GetAveragePos();
  return avgPos;
}
TVector3 trex::TTPCVolGroupMan::GetAvgPosRep(trex::TTPCUnitVolume* vol, int sign){
  TVector3 avgPos = vol->GetPos();
  return avgPos;
}

bool trex::TTPCVolGroupMan::IsInRange(trex::TTPCPathVolume* point1, trex::TTPCPathVolume* point2, int sizeX, int sizeY, int sizeZ){
  float dx = (float)(point1->GetX() - point2->GetX()) / float(sizeX);
  float dy = (float)(point1->GetY() - point2->GetY()) / float(sizeY);
  float dz = (float)(point1->GetZ() - point2->GetZ()) / float(sizeZ);

  return (dx*dx + dy*dy + dz*dz) < 1.;
  /*if( abs(point1->GetX() - point2->GetX()) > sizeX) return false;
  if( abs(point1->GetY() - point2->GetY()) > sizeY) return false;
  if( abs(point1->GetZ() - point2->GetZ()) > sizeZ) return false;
  return true;*/
}
bool trex::TTPCVolGroupMan::IsXPathCandidate(trex::TTPCOrderedVolGroup& inPath){
  // get status as x cluster
  // first check; is cluster too small in y or z?
  int minClusters = fLayout->GetMinPathClusters();
  int ySize = (inPath.GetYMax() - inPath.GetYMin()) + 1;
  int zSize = (inPath.GetZMax() - inPath.GetZMin()) + 1;
  bool sizeCheck1 = (ySize < minClusters) && (zSize < minClusters);

  if(sizeCheck1) return true;

  int maxCheckClusters = fLayout->GetXPathMaxPads();
  double minEndRatio = fLayout->GetXPathMinEndRatio();

  bool sizeCheck2 = (ySize < maxCheckClusters) && (zSize < maxCheckClusters);

  trex::TTPCPathVolume* beginVol = *(inPath.begin());
  trex::TTPCPathVolume* endVol = *(inPath.rbegin());

  double xDiff = endVol->GetX() - beginVol->GetX();
  double yDiff = endVol->GetY() - beginVol->GetY();
  double zDiff = endVol->GetZ() - beginVol->GetZ();

  double  ratio2 = (xDiff*xDiff) / ( (yDiff*yDiff) + (zDiff*zDiff) );

  bool ratioCheck = ratio2 > (minEndRatio*minEndRatio);

  if(sizeCheck2 && ratioCheck) return true;

  return false;
}


std::vector<trex::TTPCHitPad*> trex::TTPCVolGroupMan::GetHits(){
  return fPrimaryHits.GetHits();
}

bool trex::TTPCVolGroupMan::HitTest(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, trex::TTPCConnection::Type type){
  int distX;
  int distY;
  int distZ;
  fLayout->GetTypeDistances(distX, distY, distZ, type);

  return HitTest(vol1, vol2, distX, distY, distZ);
}
bool trex::TTPCVolGroupMan::HitTest(trex::TTPCUnitVolume* vol1, trex::TTPCUnitVolume* vol2, int xDist, int yDist, int zDist){
  float fx = (float)(vol2->GetX() - vol1->GetX())/(float)xDist;
  float fy = (float)(vol2->GetY() - vol1->GetY())/(float)yDist;
  float fz = (float)(vol2->GetZ() - vol1->GetZ())/(float)zDist;

  return (fx*fx + fy*fy + fz*fz <= 1.);
}
