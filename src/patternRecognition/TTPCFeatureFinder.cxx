// eddy
#include "TTPCFeatureFinder.hxx"

ND::TTPCFeatureFinder::TTPCFeatureFinder(ND::TTPCLayout* layout){
  fUse3D = true;
  fUseCharge = false;

  if(fUse3D){
    fFeatures3D = new ND::TTPCFeatureFinder3D(layout);
  }
  else{
    fFeaturesZY = new ND::TTPCFeatureFinder2D(layout, 1);
    fFeaturesZX = new ND::TTPCFeatureFinder2D(layout, 2);
    fFeaturesYX = new ND::TTPCFeatureFinder2D(layout, 3);
  };

  fLayout = layout;
  fLayout->GetRanges(fSizeX,fMinX,fMaxX, fSizeY,fMinY,fMaxY, fSizeZ,fMinZ,fMaxZ);
}
ND::TTPCFeatureFinder::~TTPCFeatureFinder(){
  if(fUse3D){
    delete fFeatures3D;
  }
  else {
    delete fFeaturesZY;
    delete fFeaturesZX;
    delete fFeaturesYX;
  };
}

void ND::TTPCFeatureFinder::AddHits(std::map<long, ND::TTPCUnitVolume*> hitMap){
  for(std::map<long, ND::TTPCUnitVolume*>::iterator el = hitMap.begin(); el != hitMap.end(); ++ el){
    double charge = 1.;
    if(fUseCharge) charge = el->second->GetQ();
    else if(el->second->GetQ() < .1) charge = 0.;

    if(fUse3D){
      fFeatures3D->AddHit(el->second->GetX(), el->second->GetY(), el->second->GetZ(), charge);
    }
    else{
      fFeaturesZY->AddHit(el->second->GetZ(), el->second->GetY(), charge);
      fFeaturesZX->AddHit(el->second->GetZ(), el->second->GetX(), charge);
      fFeaturesYX->AddHit(el->second->GetY(), el->second->GetX(), charge);
    };
  };
}
void ND::TTPCFeatureFinder::Process(){
  if(fUse3D){
    fFeatures3D->Process();
  }
  else{
    fFeaturesZY->Process();
    fFeaturesZX->Process();
    fFeaturesYX->Process();
  };
}

std::vector<ND::TTPCCell3D> ND::TTPCFeatureFinder::GetFeatures(int closeness, int window){
  std::vector<ND::TTPCCell3D> foundCells;

  /*std::vector<ND::TTPCCell2D> poisZY = featuresZY->GetLocalMaxima(window);
  std::vector<ND::TTPCCell2D> poisZX = featuresZX->GetLocalMaxima(window);
  std::vector<ND::TTPCCell2D> poisYX = featuresYX->GetLocalMaxima(window);

  std::vector<ND::TTPCCell3D> matchedZYZX;
  for(std::vector<ND::TTPCCell2D>::iterator cell1 = poisZY.begin(); cell1 != poisZY.end(); ++ cell1)
  for(std::vector<ND::TTPCCell2D>::iterator cell2 = poisZX.begin(); cell2 != poisZX.end(); ++ cell2)
  if(abs(cell1->x - cell2->x) <= closeness){
  ND::TTPCCell3D matchedCell;
  matchedCell.x = cell2->x;
  matchedCell.y = cell1->y;
  matchedCell.z = (cell1->x + cell2->x) / 2;
  matchedZYZX.push_back(matchedCell);
  };

  for(std::vector<ND::TTPCCell3D>::iterator cell1 = matchedZYZX.begin(); cell1 != matchedZYZX.end(); ++ cell1)
  for(std::vector<ND::TTPCCell2D>::iterator cell2 = poisYX.begin(); cell2 != poisYX.end(); ++ cell2)
  if(abs(cell1->x - cell2->y) <= closeness && abs(cell1->y - cell2->x) <= closeness){
  ND::TTPCCell3D foundCell;
  foundCell.x = (cell1->x + cell2->y) / 2;
  foundCell.y = (cell1->y + cell2->x) / 2;
  foundCell.z = cell1->z;
  foundCells.push_back(foundCell);
  };*/

  return foundCells;
}

float** ND::TTPCFeatureFinder::GetHitArray(int view){
  if(fUse3D){
    return fFeatures3D->GetHitArray(view);
  }
  else{
    switch(view){
      case 3:
      return fFeaturesYX->GetHitArray();
      break;
      case 2:
      return fFeaturesZX->GetHitArray();
      break;
      default:
      return fFeaturesZY->GetHitArray();
      break;
    };
  };
}
void ND::TTPCFeatureFinder::GetSizes(int& xSize, int& ySize, int& zSize){
  if(fUse3D){
    fFeatures3D->GetSizes(xSize, ySize, zSize);
  }
  else{
    fFeaturesZX->GetSizes(zSize, xSize);
    fFeaturesZY->GetSizes(zSize, ySize);
  };
}
