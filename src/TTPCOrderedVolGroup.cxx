#include "TTPCOrderedVolGroup.hxx"

//ClassImp(trex::TTPCOrderedVolGroup);
//trex::TTPCOrderedVolGroup::TTPCOrderedVolGroup(){
//  SetUp();
//}
trex::TTPCOrderedVolGroup::TTPCOrderedVolGroup(trex::TTPCLayout* layout):
  fFrontHits(layout),fBackHits(layout),fExtendedHits(layout){
  fLayout = layout;
  SetUp();
}

trex::TTPCOrderedVolGroup::TTPCOrderedVolGroup(trex::TTPCOrderedVolGroup&& in):
  fFrontHits(in.fFrontHits),fBackHits(in.fBackHits), fExtendedHits(in.fExtendedHits){

    fLayout=in.fLayout;
    fHits=std::move(in.fHits);
    in.fHits.clear();
    fAddedFrontHits=in.fAddedFrontHits;
    fAddedBackHits=in.fAddedBackHits;
    fFrontIsVertex=in.fFrontIsVertex;
    fBackIsVertex=in.fBackIsVertex;
    fHasExtendedHits=in.fHasExtendedHits;
    fXMin=in.fXMin;
    fXMax=in.fXMax;
    fYMin=in.fYMin;
    fYMax=in.fYMax;
    fZMin=in.fZMin;
    fZMax=in.fZMax;
    fIsXPath=in.fIsXPath;
    fClosed=in.fClosed;
}

trex::TTPCOrderedVolGroup::~TTPCOrderedVolGroup(){
  for(auto hit=fHits.begin();hit!=fHits.end();++hit){
    delete *hit;
  }
}

void trex::TTPCOrderedVolGroup::SetUp(){

  fAddedFrontHits = false;
  fAddedBackHits = false;
  fHasExtendedHits = false;

  fFrontIsVertex = false;
  fBackIsVertex = false;
  fIsXPath = false;

  fClosed = false;
}

void trex::TTPCOrderedVolGroup::AddFrontHits(const trex::TTPCVolGroup& frontHits){
  fFrontHits = frontHits;
  fAddedFrontHits = true;
}
void trex::TTPCOrderedVolGroup::AddBackHits(const trex::TTPCVolGroup& backHits){
  fBackHits = backHits;
  fAddedBackHits = true;
}
void trex::TTPCOrderedVolGroup::AddExtendedHits(const trex::TTPCVolGroup& extendedHits){
  fHasExtendedHits = true;
  fExtendedHits = extendedHits;
}

trex::TTPCPathVolume* trex::TTPCOrderedVolGroup::AddCell(trex::TTPCUnitVolume* cell, bool isXCluster){
  // create new path with input cell
  trex::TTPCPathVolume* path = new trex::TTPCPathVolume(cell);
  path->SetIsXCluster(isXCluster);
  // add to list of hits
  fHits.push_back(path);
  fClosed = false;

  return path;
}
void trex::TTPCOrderedVolGroup::DoClustering(bool partial){
  if(fIsXPath){
    DoXClustering();
  }
  else{
    DoStandardClustering(partial);
  };
}
void trex::TTPCOrderedVolGroup::DoStandardClustering(bool partial){
  if(!partial){
    StripX();
    FindClusterAnglesByDichotomy();
    FindHVClustersFromAngles();
  };
  InterpolateHVClusters();
  ExtrapolateHVClusters();
  FillHVClusters();
  MergeHVClusters();
}
void trex::TTPCOrderedVolGroup::DoGreedyClustering(bool partial){
  if(!partial){
    StripX();
    FindClusterAnglesByDichotomy();
    FindHVClustersFromAngles();
  };
  InterpolateHVClusters();
  ExtrapolateHVClusters();
  GreedyFillHVClusters();
}
void trex::TTPCOrderedVolGroup::DoXClustering(){
  ExpandXClusters();
  SetXClusterHVAngles();
}
void trex::TTPCOrderedVolGroup::StripX(){
  // kill adjacent nearby volumes repeated only in z
  int lastY = -1;
  int lastZ = -1;
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    int curY = pathVol->GetY();
    int curZ = pathVol->GetZ();

    if((curY == lastY) && (curZ == lastZ)){
      delete *pathVolIt;
      *pathVolIt = 0;
    };

    lastY = curY;
    lastZ = curZ;
  };
  // erase dead path volumes
  std::vector<trex::TTPCPathVolume*>::iterator pathDel = fHits.begin();
  while(pathDel != fHits.end()){
    if(!*pathDel) fHits.erase(pathDel);
    else pathDel ++;
  };
  fClosed = false;
}
void trex::TTPCOrderedVolGroup::FindClusterAnglesByDichotomy(){
  int pathSize = fHits.size();

  // define vector establishing angles along track
  std::vector<float> angles (pathSize, 9999.);
  if(pathSize < 1){
    return;
  }
  else if(pathSize < 2){
    angles = std::vector<float>(pathSize, 0);
  }
  else if(pathSize < 3){
    TVector3 tDiff = fHits[0]->GetPos() - fHits[pathSize-1]->GetPos();
    float tAngle = std::abs( std::atan(tDiff.Y() / tDiff.Z()) * 180/TMath::Pi() );
    angles = std::vector<float>(pathSize, tAngle);
  }
  else{
    // first, look for turning points to break over
    std::vector<int> turningPoints = std::vector<int>();

    // first point
    turningPoints.push_back(0);

    // turning points
    int hSteps = 0;
    int hSense = 0;
    const int threshold = fLayout->GetDichotomyCutoff();

    int prevZ = -1;
    for(int i=0; i<pathSize; i++){
      trex::TTPCPathVolume* pathVol = fHits[i];

      int z = pathVol->GetZ();

      int curHSense = 0;

      if(prevZ >= 0){
        hSteps += (z>prevZ) ? 1 : -1;
      };

      if(std::abs(hSteps) > threshold){
        curHSense = (hSteps > 0) ? 1 : -1;
        hSteps = 0;
      };

      if(curHSense){
        if((hSense * curHSense) > 0 && i>threshold) turningPoints.push_back(i-threshold); 
        hSense = curHSense;
      };

      prevZ = z;
    };

    // last point
    turningPoints.push_back(pathSize-1);

    // connect between all pairs of turning points
    for(std::vector<int>::iterator point = turningPoints.begin(); point != turningPoints.end()-1; ++point){
      RecursiveDichotomy(angles, *(point), *(point+1), 0);
    };

    // set first and last angles to match nearest defined values
    if(angles[0] > 360.){
      for(std::vector<float>::iterator angle=angles.begin(); angle != angles.end(); ++angle){
        if(*angle <= 360.){
          angles[0] = *angle;
          break;
        };
      };
    };
    if(angles[pathSize-1] > 360.){
      for(std::vector<float>::reverse_iterator angle=angles.rbegin(); angle != angles.rend(); ++angle){
        if(*angle <= 360.){
          angles[pathSize-1] = *angle;
          break;
        };
      };
    };
    if(angles[0] > 360. && angles[pathSize-1] > 360.) angles = std::vector<float>(pathSize, 0);
  };

  // loop through angles, defining unknown values to interpolate between known values
  int undefAngles=0;
  float fAngle=9999.;
  float lAngle=9999.;
  for(std::vector<float>::iterator angle=angles.begin(); angle != angles.end(); ++angle){
    bool redef=false;

    if(*angle > 360.){
      undefAngles++;

      if(angle == angles.end()-1){
        redef=true;
        lAngle = fAngle;
      };
    }
    else{
      lAngle = *angle;
      if(fAngle > 360.) fAngle = *angle;

      redef=true;
    };
    if(redef){
      if(undefAngles > 0){
        float step=(lAngle-fAngle) / float(undefAngles+1);
        float val=fAngle;
        for(std::vector<float>::iterator edAngle=angle-undefAngles; edAngle!=angle; ++edAngle){
          val += step;
          *edAngle = val;
        };

        if(angle==angles.end()-1) *angle=fAngle;

        undefAngles=0;
        fAngle = *angle;
      };
    };
  };
  // set angles
  for(int i=0; i<pathSize; i++){
    float angle = angles[i];
    trex::TTPCPathVolume* pathVol = fHits[i];

    pathVol->SetPatRecAngle(angle);
  };
}
void trex::TTPCOrderedVolGroup::RecursiveDichotomy(std::vector<float>& angles, int firstID, int lastID, float prevAngDiff){
  int midPoint = int( (firstID+lastID)/2 );
  if(((lastID-firstID) < fLayout->GetDichotomyCutoff()) && (firstID != 0 || lastID != (int)angles.size()-1)) return;
  if(midPoint == firstID || midPoint == lastID) return;

  TVector3 pDiff1 = fHits[firstID]->GetPos() - fHits[midPoint]->GetPos();
  float angle1 = std::atan2(pDiff1.Y(), pDiff1.Z());

  TVector3 pDiff2 = fHits[lastID]->GetPos() - fHits[midPoint]->GetPos();
  float angle2 = std::atan2(pDiff2.Y(), pDiff2.Z());

  float angDiff = std::abs(angle2 - angle1);
  if(angDiff > TMath::Pi()) angDiff = ( 2*TMath::Pi() ) - angDiff;

  if(angDiff >= prevAngDiff){
    TVector3 pDiff = fHits[lastID]->GetPos() - fHits[firstID]->GetPos();
    float angle = std::abs( std::atan(pDiff.Y() / pDiff.Z()) * 180/TMath::Pi() );

    angles[midPoint] = angle;

    RecursiveDichotomy(angles, firstID, midPoint, angDiff);
    RecursiveDichotomy(angles, midPoint, lastID, angDiff);
  };
}
void trex::TTPCOrderedVolGroup::FindHVClustersFromAngles(){
  int pathSize = fHits.size();

  // relevant clustering parameters
  int chainThreshold = fLayout->GetDichotomyCutoff();
  int edgeThreshold = fLayout->GetHVEdgeDist();
  float thresholdAngle = fLayout->GetThresholdAngle();

  // make list of all runs of horizontal or vertical clusters
  std::vector< std::pair<char, int> > clusterRuns;

  char curOr = '\0';
  int curStep = 1;
  float maxAngl = 0.0;
  float minAngl = 0.0;
  if(pathSize>0)
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    if(pathVolIt == fHits.begin()){
      maxAngl = pathVol->GetPatRecAngle();
      minAngl = pathVol->GetPatRecAngle();
    }

    maxAngl = std::max(maxAngl, pathVol->GetPatRecAngle());
    minAngl = std::min(minAngl, pathVol->GetPatRecAngle());
    float diffAngl = maxAngl - minAngl;
    char nextOr = (std::abs(pathVol->GetPatRecAngle()) > thresholdAngle ) ? 'h' : 'v';
    if(!curOr) curOr = nextOr;

    if( (nextOr != curOr && fabs(diffAngl) > 5.0) || pathVolIt == fHits.end()-1){
      clusterRuns.push_back(std::make_pair(curOr, curStep));

      curOr = nextOr;
      curStep = 1;

      maxAngl = pathVol->GetPatRecAngle();
      minAngl = pathVol->GetPatRecAngle();
    } else{
      curStep ++;
    };
  };

  // remove short runs, starting with the shortest
  bool possibleShort = true;
  if(clusterRuns.size() > 1)
  while(possibleShort){
    possibleShort = false;

    int shortestLength = 9999;
    std::vector< std::pair<char, int> >::iterator shortestIt;
    for(std::vector< std::pair<char, int> >::iterator runIt = clusterRuns.begin(); runIt != clusterRuns.end(); ++runIt){
      if(runIt->second < shortestLength){
        shortestLength = runIt->second;
        shortestIt = runIt;
      };
    };

    // if short track found, merge into neighbors
    bool mergeMe = false;
    if(shortestLength < chainThreshold){
      mergeMe = true;
    }
    else if(shortestLength < edgeThreshold){
      if(shortestIt == clusterRuns.begin() || shortestIt == (clusterRuns.end()-1)){
        mergeMe = true;
      };
    };

    if (mergeMe){
      if(shortestIt == clusterRuns.begin()){
        (shortestIt+1)->second = shortestIt->second + (shortestIt+1)->second;
        clusterRuns.erase(shortestIt);
      }
      else if(shortestIt == clusterRuns.end()-1){
        (shortestIt-1)->second = (shortestIt-1)->second + shortestIt->second;
        clusterRuns.erase(shortestIt);
      }
      else{
        (shortestIt-1)->second = (shortestIt-1)->second + shortestIt->second + (shortestIt+1)->second;
        clusterRuns.erase(shortestIt, shortestIt+2);
      };
      possibleShort = true;
    };
  };

  std::vector< std::pair<char, int> >::iterator shortestIt;
  int i=0;
  for(std::vector< std::pair<char, int> >::iterator runIt = clusterRuns.begin(); runIt != clusterRuns.end(); ++runIt){
    for(int j=0; j<runIt->second; j++){
      if(runIt->first == 'v') fHits[i]->SetIsVertical(true);
      else if(runIt->first == 'h') fHits[i]->SetIsVertical(false);
      i++;
    };
  };
}
void trex::TTPCOrderedVolGroup::ExtrapolateHVClusters(){
  if(fHits.size() < 2) return;

  // front and back hits
  std::vector<trex::TTPCPathVolume*> backAddHits = GetExtrapolatedClusters(fHits.begin(), -1);
  std::vector<trex::TTPCPathVolume*> frontAddHits = GetExtrapolatedClusters(fHits.end()-1, 1);

  // reverse front hits and insert to front
  std::vector<trex::TTPCPathVolume*> orderedBackHits;
  for(std::vector<trex::TTPCPathVolume*>::reverse_iterator pathVolRIt = backAddHits.rbegin(); pathVolRIt != backAddHits.rend(); ++pathVolRIt){
    orderedBackHits.push_back(*pathVolRIt);
  };
  fHits.insert(fHits.begin(), orderedBackHits.begin(), orderedBackHits.end());

  // insert end hits to the back
  fHits.insert(fHits.end(), frontAddHits.begin(), frontAddHits.end());

  fClosed = false;
}
void trex::TTPCOrderedVolGroup::InterpolateHVClusters(){
  if(fHits.size() < 2) return;

  int closestDist2 = fLayout->GetClusterConnectDist()*fLayout->GetClusterConnectDist();
  // look for pairs of hits where one is skipped in between
  for(unsigned int i = 0; i < fHits.size()-1; ++i){
    trex::TTPCPathVolume* pathVol1 = fHits.at(i);
    trex::TTPCPathVolume* pathVol2 = fHits.at(i+1);

    if(!pathVol1) continue;
    if(!pathVol2) continue;

    // make sure the paths are of the same orientation
    if(pathVol1->GetIsVertical() == pathVol2->GetIsVertical()){
      bool isVertical = pathVol1->GetIsVertical();

      // check if there's a gap in between them
      int diff;
      if(isVertical){
        diff = std::abs(pathVol2->GetZ() - pathVol1->GetZ());
      }
      else{
        diff = std::abs(pathVol2->GetY() - pathVol1->GetY());
      };

      if(diff > 1){
        // fill vector of points to insert
        std::vector<trex::TTPCPathVolume*> toInsert;

        // try and fill in points in between
        int dX = (pathVol2->GetX() - pathVol1->GetX()) / diff;
        int dY = (pathVol2->GetY() - pathVol1->GetY()) / diff;
        int dZ = (pathVol2->GetZ() - pathVol1->GetZ()) / diff;

        for(int check=1; check<diff; ++check){
          int targX = pathVol1->GetX() + (dX * check);
          int targY = pathVol1->GetY() + (dY * check);
          int targZ = pathVol1->GetZ() + (dZ * check);

          // look for a hit within this distance
          trex::TTPCUnitVolume* closestVol = 0;
          int checkClosestDist2 = closestDist2;
          for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = fExtendedHits.begin(); volEl != fExtendedHits.end(); ++volEl){
            trex::TTPCUnitVolume* vol = volEl->second;

            int distX = targX - vol->GetX();
            int distY = targY - vol->GetY();
            int distZ = targZ - vol->GetZ();
            int dist2;
            if(isVertical){
              if(distZ == 0){
                dist2 = distX*distX + distY*distY;
                if(dist2 <= checkClosestDist2){
                  closestVol = vol;
                  checkClosestDist2 = closestDist2;
                };
              };
            }
            else{
              if(distY == 0){
                dist2 = distX*distX + distZ*distZ;
                if(dist2 <= checkClosestDist2){
                  closestVol = vol;
                  checkClosestDist2 = closestDist2;
                };
              };
            };
          };
          if(closestVol){
            trex::TTPCPathVolume* newPathVol = new trex::TTPCPathVolume(closestVol);
            newPathVol->SetIsVertical(isVertical);

            toInsert.push_back(newPathVol);
          };
        };
        if(toInsert.size()){
          // now add all elements from the insertable
          fHits.insert(fHits.begin()+(i+1), toInsert.begin(), toInsert.end());

          // update index to continue iterating
          i += toInsert.size();
        };
      };
    };
  };
  fClosed = false;
}

void trex::TTPCOrderedVolGroup::FillHVClusters(){
  // loop over all points in path
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    bool isVertical = pathVol->GetIsVertical();

    // grab co-ordinates
    int pathX = pathVol->GetX();
    int pathY = pathVol->GetY();
    int pathZ = pathVol->GetZ();

    // set up cluster connection size and scaling factor in x direction
    int incSizeYZ = fLayout->GetClusterConnectDist();
    float xFact = (float)(fLayout->GetConnectDistX()) / (float)(fLayout->GetConnectDistY());
    int incSizeX = int(incSizeYZ*xFact + 1);

    // size of bounding box for adding hits
    int xMin = pathX - incSizeX;
    int xMax = pathX + incSizeX;
    int yzMin = (isVertical ? pathY : pathZ) - incSizeYZ;
    int yzMax = (isVertical ? pathY : pathZ) + incSizeYZ;

    // keep adding hits and expanding box to check until no more are found
    bool hitsFound = true;
    while(hitsFound){
      hitsFound = false;

      for(int i=xMin; i<=xMax; i++)
      for(int jk=yzMin; jk<=yzMax; jk++){
        int x = i;
        int y = (isVertical ? jk : pathY);
        int z = (isVertical ? pathZ : jk);

        // check if hit exists
        long cellID = fLayout->SafeMash(x, y, z);
        if(cellID < 0) continue;

        trex::TTPCUnitVolume* vol = fExtendedHits.GetEl(cellID);
        if(!vol) continue;

        // check that his isn't already added
        if(pathVol->GetFriendsContains(vol)) continue;

        // add hit and continue while loop over larger range if something is found
        pathVol->AddFriend(vol);

        xMin = std::min(pathX-incSizeX, vol->GetX()-incSizeX);
        xMax = std::max(pathX+incSizeX, vol->GetX()+incSizeX);
        if(isVertical){
          yzMin = std::min(pathY-incSizeYZ, vol->GetY()-incSizeYZ);
          yzMax = std::max(pathY+incSizeYZ, vol->GetY()+incSizeYZ);
        }
        else{
          yzMin = std::min(pathZ-incSizeYZ, vol->GetZ()-incSizeYZ);
          yzMax = std::max(pathZ+incSizeYZ, vol->GetZ()+incSizeYZ);
        };
        hitsFound = true;
      };
    };
  };
  fClosed = false;
}
void trex::TTPCOrderedVolGroup::MergeHVClusters(){
  // first, clear any empty clusters
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    if(!pathVol->GetHasCluster()){
      delete pathVol;
      *pathVolIt = 0;
    };
  };
  // loop until no more merges can be made
  bool canMerge = true;
  while(canMerge){
    canMerge = false;
    // compare pairs of path volumes 
    for(std::vector<trex::TTPCPathVolume*>::iterator pathVol1It = fHits.begin(); pathVol1It != fHits.end(); ++pathVol1It)
    for(std::vector<trex::TTPCPathVolume*>::iterator pathVol2It = pathVol1It+1; pathVol2It != fHits.end(); ++pathVol2It){
      trex::TTPCPathVolume* pathVol1 = *pathVol1It;
      trex::TTPCPathVolume* pathVol2 = *pathVol2It;
      // break out if either path fails to exist
      if(!pathVol1 || !pathVol2) continue;

      // only proceed if both are oriented same way in same layer
      if(pathVol1->GetIsVertical() && pathVol2->GetIsVertical()){
        if (pathVol1->GetZ() != pathVol2->GetZ()) continue;
      }
      else if(!pathVol1->GetIsVertical() && !pathVol2->GetIsVertical()){
        if (pathVol1->GetY() != pathVol2->GetY()) continue;
      }
      else{
        continue;
      };

      bool olFound = false;
      for(std::vector<trex::TTPCUnitVolume*>::iterator vol1It = pathVol1->GetFriendsBegin(); vol1It != pathVol1->GetFriendsEnd(); ++vol1It){
        for(std::vector<trex::TTPCUnitVolume*>::iterator vol2It = pathVol2->GetFriendsBegin(); vol2It != pathVol2->GetFriendsEnd(); ++vol2It){
          if((*vol1It) == (*vol2It)){
            olFound = true;
            break;
          };
        };
        if(olFound) break;
      };
      // if overlap found, merge clusters into one
      if(olFound){
        canMerge = true;
        for(std::vector<trex::TTPCUnitVolume*>::iterator volIt = pathVol2->GetFriendsBegin(); volIt != pathVol2->GetFriendsEnd(); ++volIt){
          pathVol1->SafeAddFriend(*volIt);
        };
        delete pathVol2;
        *pathVol2It = 0;
      };
    };
  };

  // erase dead path volumes
  std::vector<trex::TTPCPathVolume*>::iterator pathDel = fHits.begin();
  while(pathDel != fHits.end()){
    if(!*pathDel) fHits.erase(pathDel);
    else pathDel ++;
  };
  fClosed = false;
}
void trex::TTPCOrderedVolGroup::GreedyFillHVClusters(){
  // fill maps of clusters based on their positions
  std::map<int, std::vector<trex::TTPCPathVolume*> > vClusters;
  std::map<int, std::vector<trex::TTPCPathVolume*> > hClusters;

  // loop over all points in path
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); pathVolIt++){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    bool isVertical = pathVol->GetIsVertical();

    // store h clusters and v clusters separately to avoid duplicates
    if(isVertical){
      int bin = pathVol->GetZ();
      if(vClusters.find(bin) == vClusters.end()){
        vClusters[bin] = std::vector<trex::TTPCPathVolume*>();
      };
      vClusters[bin].push_back(pathVol);
    }
    else{
      int bin = pathVol->GetY();
      if(hClusters.find(bin) == hClusters.end()){
        hClusters[bin] = std::vector<trex::TTPCPathVolume*>();
      };
      hClusters[bin].push_back(pathVol);
    };
  };

  // merge adjacent clusters
  GreedyMergeHVClusters(vClusters);
  GreedyMergeHVClusters(hClusters);

  // try and add all hits to horizontal and vertical cluster
  for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = fExtendedHits.begin(); volEl != fExtendedHits.end(); ++volEl){
    GreedyAddVol(volEl->second, vClusters, true);
    GreedyAddVol(volEl->second, hClusters, false);
  };
}
void trex::TTPCOrderedVolGroup::GreedyMergeHVClusters(std::map<int, std::vector<trex::TTPCPathVolume*> >& clusters){
  int clusterMergeThreshold = fLayout->GetClusterConnectDist();
  int clusterMergeThreshold2 = clusterMergeThreshold*clusterMergeThreshold;

  for(std::map<int, std::vector<trex::TTPCPathVolume*> >::iterator clEl = clusters.begin(); clEl != clusters.end(); ++clEl){
    std::vector<trex::TTPCPathVolume*> oldCl = clEl->second;
    std::vector<trex::TTPCPathVolume*> newCl;

    // merge everything in a cell together if in range
    while(oldCl.size()){

      std::vector<trex::TTPCPathVolume*> mergeSet;
      bool canMerge = true;
      while(canMerge){
        canMerge = false;
        for(std::vector<trex::TTPCPathVolume*>::iterator oldVolIt = oldCl.begin(); oldVolIt != oldCl.end(); ++oldVolIt){
          trex::TTPCPathVolume* oldVol = *oldVolIt;
          if(!oldVol) continue;

          // if there's nothing to check yet, add this
          if(!mergeSet.size()){
            mergeSet.push_back(oldVol);
            *oldVolIt = 0;
          }
          // otherwise, add if it's connected to anything in the merge set
          else{
            for(std::vector<trex::TTPCPathVolume*>::iterator mergeIt = mergeSet.begin(); mergeIt != mergeSet.end(); ++mergeIt){
              trex::TTPCPathVolume* merge = *mergeIt;

              int dX = merge->GetX() - oldVol->GetX();
              int dY = merge->GetY() - oldVol->GetY();
              int dZ = merge->GetZ() - oldVol->GetZ();

              // merge if within range of threshold
              if(dX*dX + dY*dY + dZ*dZ <= clusterMergeThreshold2){
                mergeSet.push_back(oldVol);
                *oldVolIt = 0;
                canMerge = true;

                break;
              };
            };
          };
        };
      };
      // if there's stuff to merge, merge it into one vol and return
      if(mergeSet.size()){
        newCl.push_back(GetAverageVol(mergeSet));
      };

      // clear dead vols
      std::vector<trex::TTPCPathVolume*>::iterator clDel = oldCl.begin();
      while(clDel != oldCl.end()){
        if(!*clDel) oldCl.erase(clDel);
        else clDel ++;
      };
    };
    clEl->second = newCl;
  };
}
trex::TTPCPathVolume* trex::TTPCOrderedVolGroup::GetAverageVol(std::vector<trex::TTPCPathVolume*> clusters){
  trex::TTPCPathVolume* avgVol = 0;

  // find average
  float avgX = 0.;
  float avgY = 0.;
  float avgZ = 0.;
  float weight = 0.;
  for(std::vector<trex::TTPCPathVolume*>::iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt){
    trex::TTPCPathVolume* cluster = *clusterIt;

    avgX += cluster->GetX();
    avgY += cluster->GetY();
    avgZ += cluster->GetZ();
    weight += 1.;
  };
  avgX /= weight;
  avgY /= weight;
  avgZ /= weight;

  // get path vol closest to average
  float minR2 = 99999999.;
  for(std::vector<trex::TTPCPathVolume*>::iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt){
    trex::TTPCPathVolume* cluster = *clusterIt;

    float dX = cluster->GetX() - avgX;
    float dY = cluster->GetY() - avgY;
    float dZ = cluster->GetZ() - avgZ;
    float dR2 = dX*dX + dY*dY + dZ*dZ;

    if(dR2 < minR2){
      minR2 = dR2;
      avgVol = cluster;
    };
  };

  return avgVol;
}
bool trex::TTPCOrderedVolGroup::GreedyAddVol(trex::TTPCUnitVolume* vol, std::map<int, std::vector<trex::TTPCPathVolume*> >& clusters, bool isVertical){
  // add to nearest cluster in col
  int bin;
  if(isVertical){
    bin = vol->GetZ();
  }
  else{
    bin = vol->GetY();
  };

  if(clusters.find(bin) != clusters.end()){
    std::vector<trex::TTPCPathVolume*> cluster = clusters[bin];

    // add to nearest cluster in bin
    float minDist2 = 99999999.;
    trex::TTPCPathVolume* closestVol = 0;
    for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = cluster.begin(); pathVolIt != cluster.end(); ++pathVolIt){
      trex::TTPCPathVolume* pathVol = *pathVolIt;

      float dX = vol->GetX() - pathVol->GetX();
      float dY = vol->GetY() - pathVol->GetY();
      float dZ = vol->GetZ() - pathVol->GetZ();
      float dR2 = dX*dX + dY*dY + dZ*dZ;
      if(dR2 < minDist2){
        minDist2 = dR2;
        closestVol = pathVol;
      };
    };

    if(closestVol){
      closestVol->AddFriend(vol);
      return true;
    };
  };
  return false;
}

void trex::TTPCOrderedVolGroup::ExpandXClusters(){
  if(!fHits.size()) return;

  int xMin = GetXMin();
  int xMax = GetXMax();
  int zMin = fExtendedHits.GetZMin();
  int zMax = fExtendedHits.GetZMax();

  // determing if pointing up (normal order) or down (inverted)
  bool normalOrder = (*fHits.begin())->GetX() < (*fHits.rbegin())->GetX();

  // delete old clusters
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVol = fHits.begin(); pathVol != fHits.end(); ++pathVol){
    delete *pathVol;
  };
  fHits.clear();

  // add hits at desired x positions
  for(int x=xMin; x<=xMax; x++){
    for(int z=zMin; z<=zMax; z++){
      trex::TTPCPathVolume* hit = 0;
      // either add hits to current path volumes or create new ones
      for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = fExtendedHits.begin(); volEl != fExtendedHits.end(); ++volEl){
        trex::TTPCUnitVolume* vol = volEl->second;

        if(vol->GetX() != x) continue;
        if(vol->GetZ() != z) continue;

        if(!hit){
          hit = AddCell(vol, true);
        };
        hit->AddFriend(vol);
      };
    };
  };

  // flip to order away from path
  if(!normalOrder){
    Flip();
  };
};
void trex::TTPCOrderedVolGroup::SetXClusterHVAngles(){
  for(std::vector<trex::TTPCPathVolume*>::iterator hitIt = fHits.begin(); hitIt != fHits.end(); ++hitIt){
    trex::TTPCPathVolume* hit = *hitIt;
    hit->SetPatRecAngle(0);
    hit->SetIsVertical(true);
  };
}

void trex::TTPCOrderedVolGroup::Clean(){
  std::vector<trex::TTPCPathVolume*>::iterator reaper;
  reaper = begin();
  while(reaper != end()){
    if(!(*reaper)) erase(reaper);
    else if(!(*reaper)->GetHasCluster()) erase(reaper);
    else(reaper++);
  };
}

void trex::TTPCOrderedVolGroup::Flip(){
  // temporary copies

  trex::TTPCVolGroup newBackHits = fFrontHits;
  bool newFrontIsVertex = fBackIsVertex;
  bool newBackIsVertex = fFrontIsVertex;
  std::vector<trex::TTPCPathVolume*> newHits;

  // fill new hits
  for(std::vector<trex::TTPCPathVolume*>::reverse_iterator pathVolRit = fHits.rbegin(); pathVolRit != fHits.rend(); ++pathVolRit){
    trex::TTPCPathVolume* pathVol = *pathVolRit;
    newHits.push_back(pathVol);
  };

  trex::TTPCVolGroup hitsTemp = fFrontHits;
  fFrontHits = fBackHits;
  fBackHits = hitsTemp;

  fFrontIsVertex = newFrontIsVertex;
  fBackIsVertex = newBackIsVertex;
  fHits = std::move(newHits);
}
void trex::TTPCOrderedVolGroup::OrderForwardsDirection(){
  if(!fHits.size()) return;

  // check if front hit is ahead of back hit and flip if not
  trex::TTPCPathVolume* frontVol = *(fHits.rbegin());
  trex::TTPCPathVolume* backVol = *(fHits.begin());

  if(frontVol->GetZ() < backVol->GetZ()){
    Flip();
  };
}
void trex::TTPCOrderedVolGroup::OrderNegativeCurvature(){
  if(!fHits.size()) return;

  // check if curvature corresponds to forwards going negative and flip if not
  trex::TTPCPathVolume* frontVol = *(fHits.rbegin());
  trex::TTPCPathVolume* backVol = *(fHits.begin());

  TVector3 startPos = frontVol->GetPos();
  TVector3 endPos = backVol->GetPos();
  TVector3 parNorm = (startPos - endPos).Unit();

  double maxDist = -1;
  int distCurvature = 0;
  // determine if furthest hit in perpendicular direction is in positive or negative direction
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;

    TVector3 diff = startPos - pathVol->GetPos();
    TVector3 dist = diff - (parNorm*diff.Dot(parNorm));

    // curvature from first hit
    if(dist.Mag() > maxDist){
      maxDist = dist.Mag();
      distCurvature = (dist.Y() > 0.) ? -1 : 1;
    };
  };

  // if curving positively relative to first hit, flip
  if(distCurvature > 0){
    Flip();
  };
}

void trex::TTPCOrderedVolGroup::OrderFromJunction(trex::TTPCVolGroup& junction){
  OrderFromPosition(junction.GetAveragePosition());
}

void trex::TTPCOrderedVolGroup::OrderFromPosition(TVector3 pos){
  if(!fHits.size()) return;

  // check if curvature corresponds to forwards going negative and flip if not
  trex::TTPCPathVolume* frontVol = *(fHits.end()-1);
  trex::TTPCPathVolume* backVol = *(fHits.begin());

  TVector3 frontDist = frontVol->GetPos() - pos;
  TVector3 backDist = backVol->GetPos() - pos;

  if(frontDist.Mag() > backDist.Mag()){
    Flip();
  };
}


std::vector<trex::TTRExHVCluster*> trex::TTPCOrderedVolGroup::GetClusters(){

  if(!fClusters.size()){

    for(std::vector<trex::TTPCPathVolume*>::iterator id = fHits.begin(); id != fHits.end(); ++id){
      
      std::vector<trex::TTPCHitPad*> hits = (*id)->GetHits();
      
      //bool hasCluster = (*id)->GetHasCluster();
      //std::cout << "This hit has clustered friends" << std::endl;
      bool IsVertical = (*id)->GetIsVertical();
      TVector3 pos = (*id)->GetPos();
      double charge = (*id)->GetQ();
      
      fClusters.emplace_back();
      
      trex::TTRExHVCluster& HVcluster=fClusters.back();
      HVcluster.SetHits(hits);
      HVcluster.SetIsVertical(IsVertical);
      HVcluster.SetPosition(pos);
      HVcluster.SetCharge(charge);

    }
  }

  std::vector<trex::TTRExHVCluster*> retClusters;

  for(auto iter=fClusters.begin();iter!=fClusters.end();++iter){
    retClusters.push_back(&*iter);
  }

  return std::move(retClusters);
}

std::string trex::TTPCOrderedVolGroup::GetOrientations(){
  std::string orientationsList = "";
  if(fHits.size() > 0)
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    orientationsList += (*pathVolIt)->GetIsVertical() ? 'v' : 'h';
    if(pathVolIt != fHits.end()-1) orientationsList += ' ';
  };
  return orientationsList;
}
void trex::TTPCOrderedVolGroup::PrintOrientations(){
  std::cout << "  " << GetOrientations() << std::endl;
}
void trex::TTPCOrderedVolGroup::PrintPositions(bool size, bool clusterHits){
  std::cout << "  ";
  if(fHits.size() > 0)
  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    std::cout << "(" << pathVol->GetX() << ", " << pathVol->GetY() << ", " << pathVol->GetZ() << ")";
    if(clusterHits) pathVol->PrintPositions(true);
    else if(size) std::cout << "{" << pathVol->GetClusterSize() << "} ";
    else std::cout << " ";
  };
  std::cout << std::endl;
}
void trex::TTPCOrderedVolGroup::CheckClusters(){
  if(fHasExtendedHits){
    for(std::map<long, trex::TTPCUnitVolume*>::iterator volEl = fExtendedHits.begin(); volEl != fExtendedHits.end(); ++volEl){
      trex::TTPCUnitVolume* vol = volEl->second;
      if(!vol){
        std::cout << "WARNING:  empty trex::TTPCUnitVolume* found in TTPCPathVolume extended hits" << std::endl;
      }
      else{
        std::cout << "  attempting to access trex::TTPCOrderedVolGroup trex::TTPCUnitVolume* hits…" << std::endl;
        std::vector<trex::TTPCHitPad*> hits = vol->GetHits();
        std::cout << "  …got hits!" << std::endl;
        for(std::vector< trex::TTPCHitPad* >::iterator hitIt = vol->GetHitsBegin(); hitIt != vol->GetHitsEnd(); ++hitIt){
          trex::TTPCHitPad* nhit = *hitIt;
          if (!nhit){
            std::cout << "WARNING:  non trex::TTPCHitPad hit found in TTPCOrderedVolGroup extended hits" << std::endl;
          }
        }
      }
    }
  }

  for(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt = fHits.begin(); pathVolIt != fHits.end(); ++pathVolIt){
    trex::TTPCPathVolume* pathVol = *pathVolIt;
    std::cout << "checking trex::TTPCPathVol* hits…" << std::endl;
    pathVol->CheckHits();
    std::cout << "…done checking trex::TTPCPathVol* hits" << std::endl;
  }
}

std::vector<trex::TTPCPathVolume*> trex::TTPCOrderedVolGroup::GetExtrapolatedClusters(std::vector<trex::TTPCPathVolume*>::iterator pathVolIt, int dir){
  std::vector<trex::TTPCPathVolume*>::iterator pathVol2It;
  trex::TTPCPathVolume* pathVol = *pathVolIt;
  trex::TTPCPathVolume* pathVol2 = 0;

  std::vector<trex::TTPCPathVolume*> outClusters;

  // don't try to extrapolate into vertex
  if((dir == 1 && fFrontIsVertex) || (dir == -1 && fBackIsVertex)) return outClusters;

  // get rough extrapolation direction from first hit to be in different plane
  int diff = 0;
  int i = 0;
  while (diff == 0){
    if(dir == 1){
      pathVol2It = pathVolIt-i;
      if (pathVol2It == fHits.begin()) return outClusters;
    }
    else if(dir == -1){
      pathVol2It = pathVolIt+i;
      if (pathVol2It == (fHits.end()-1)) return outClusters;
    };
    i++;

    pathVol2 = *pathVol2It;
    if(pathVol->GetIsVertical()) diff = pathVol->GetZ() - pathVol2->GetZ();
    else diff = pathVol->GetY() - pathVol2->GetY();
  };

  int extDir = (diff < 0) ? -1 : 1;

  // loop, setting limit on maximum number of iterations
  trex::TTPCPathVolume* curVol = pathVol;
  float xFact = (float)(fLayout->GetConnectDistX()) / (float)(fLayout->GetConnectDistY());
  int unfoundCount = fLayout->GetHVClusterExtrapolateDist();
  int limit = fLayout->GetHVClusterExtrapolateLimit();

  bool isVertical = curVol->GetIsVertical();

  int extX = curVol->GetX();
  int extY = curVol->GetY();
  int extZ = curVol->GetZ();
  for(int i=0; i<limit; ++i){

    if(isVertical) extZ += extDir;
    else extY += extDir;

    // look for a hit within range of extrapolated position
    trex::TTPCUnitVolume* closestHit = 0;
    int closestDist2 = fLayout->GetClusterConnectDist()*fLayout->GetClusterConnectDist();
    for(std::map<long, trex::TTPCUnitVolume*>::iterator extHitEl = fExtendedHits.begin(); extHitEl != fExtendedHits.end(); ++extHitEl){
      trex::TTPCUnitVolume* extHit = extHitEl->second;
      int diffX = (extHit->GetX() - extX)/xFact;
      int diffY = extHit->GetY() - extY;
      int diffZ = extHit->GetZ() - extZ;

      int dist2;
      if(isVertical){
        if(diffZ) continue;
        dist2 = diffX*diffX + diffY*diffY;
      }
      else {
        if(diffY) continue;
        dist2 = diffX*diffX + diffZ*diffZ;
      };

      if(dist2 < closestDist2){
        closestHit = extHit;
        closestDist2 = dist2;
      };
    };

    // if a hit is found, add it as a new path volume and continue the search
    if(closestHit){
      trex::TTPCPathVolume* newVol = new trex::TTPCPathVolume(closestHit);
      newVol->SetIsVertical(isVertical);

      outClusters.push_back(newVol);

      curVol = newVol;

      // reset other variables
      isVertical = curVol->GetIsVertical();

      extX = curVol->GetX();
      extY = curVol->GetY();
      extZ = curVol->GetZ();

      unfoundCount = fLayout->GetHVClusterExtrapolateDist();
    }
    else{
      if(unfoundCount > 0){
        unfoundCount --;
      }
      else{
        break;
      };
    };
  };

  return outClusters;
}

void trex::TTPCOrderedVolGroup::Close(){
  if(fClosed) return;

  // get overall extents
  fXMin = 99999;
  fXMax = -99999;
  fYMin = 99999;
  fYMax = -99999;
  fZMin = 99999;
  fZMax = -99999;
  for(std::vector<trex::TTPCPathVolume*>::iterator hitIt = fHits.begin(); hitIt != fHits.end(); ++hitIt){
    trex::TTPCPathVolume* hit = *hitIt;
    if(!hit) continue;

    fXMin = std::min(fXMin, hit->GetX());
    fXMax = std::max(fXMax, hit->GetX());
    fYMin = std::min(fYMin, hit->GetY());
    fYMax = std::max(fYMax, hit->GetY());
    fZMin = std::min(fZMin, hit->GetZ());
    fZMax = std::max(fZMax, hit->GetZ());
  };

  fClosed = true;
}
