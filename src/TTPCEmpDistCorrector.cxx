#include "TTPCEmpDistCorrector.hxx"
#include "TOARuntimeParameters.hxx"
#include "TGeomInfo.hxx"
#include "HEPUnits.hxx"
#include <fstream>
#include <string>
#include <iostream>

using namespace std;

// *****************************************************************************
TTPCEmpDistCorrector::TTPCEmpDistCorrector() {
  fLaserHitsAndDistortions.clear();
  fDistMapIsSet = false;
  // Get a reference to the singleton instance of TPC geometry information
  ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();
  fMaxDriftInX = tpcGeom.GetMaxDriftDistance();
  fGlobalXOfPadPlane = tpcGeom.GetMaxDriftDistance() + tpcGeom.GetCathodeWidth() / 2.;
  //Flag to use four neighbors, weighted average -> Disable as default
  fFourHits = false;
}

// *****************************************************************************
TTPCEmpDistCorrector::~TTPCEmpDistCorrector() {
}


// *****************************************************************************
bool TTPCEmpDistCorrector::SetDistMap(TString distMap) {
  if (fDistMapIsSet) {
    cout<<"ERROR in TTPCEmpDistCorrector::SetDistMap, distortion map already set!"<<endl;
    return false;
  }
  // Read in distortion file
  string line;
  double xAl, yAl, zAl, xDist, yDist, zDist, tDist;
  ifstream inputFile(distMap);
  if (!inputFile.is_open()) {
    cout<<"ERROR in TTPCEmpDistCorrector::SetDistMap, distortion map could not open!"<<endl;
    return false;
  }
  fLaserHitsAndDistortions.clear();
  //Get the header line
  getline(inputFile,line);
  //Read in all lines
  while (!inputFile.eof()) {
    inputFile >> xAl;
    inputFile >> yAl;
    inputFile >> zAl;
    inputFile >> xDist;
    inputFile >> yDist;
    inputFile >> zDist;
    inputFile >> tDist;
    TVector3 laserHit(xAl, yAl, zAl);
    TLorentzVector distortion(xDist, yDist, zDist, tDist);
    pair<TVector3, TLorentzVector> laserHitAndDistortion(laserHit, distortion);
    fLaserHitsAndDistortions.push_back(laserHitAndDistortion);
  }
  inputFile.close();
  fDistMapIsSet = true;
  return true;
}
//lorena

// *****************************************************************************
bool TTPCEmpDistCorrector::SetAllDistMap() {
  if (fDistMapIsSet) {
    cout<<"ERROR in TTPCEmpDistCorrector::SetDistMap, distortion map already set!"<<endl;
    return false;
  }
  // Read in distortion files
  string line;
  double xAl, yAl, zAl, xDist, yDist, zDist, tDist;
  fLaserHitsAndDistortions.clear();
  char distmap[128];

  for (int TPCnum = 1; TPCnum < 4 ; TPCnum++){
    for (int RPnum = 0; RPnum < 2 ; RPnum++){
      
      sprintf(distmap,"DistortionMap_TPC%d_RP%d.txt",TPCnum,RPnum);
      TString filename  = getenv("TREXRECONROOT") + TString("/parameters/") + TString(distmap);
      std::cout << "Distortion Map: file name " << filename << std::endl;
      ifstream inputFile(filename);
      if (!inputFile.is_open()) {
      cout<<"ERROR in TTPCEmpDistCorrector::SetAllDistMap, distortion map could not open!"<<endl;
      return false;
      }
      
      //Get the header line
      getline(inputFile,line);
      //Read in all lines
      while (!inputFile.eof()) {
	inputFile >> xAl;
	inputFile >> yAl;
	inputFile >> zAl;
	inputFile >> xDist;
	inputFile >> yDist;
	inputFile >> zDist;
	inputFile >> tDist;
	TVector3 laserHit(xAl, yAl, zAl);
	TLorentzVector distortion(xDist, yDist, zDist, tDist);
	pair<TVector3, TLorentzVector> laserHitAndDistortion(laserHit, distortion);
	fLaserHitsAndDistortions.push_back(laserHitAndDistortion);
      }
      inputFile.close();
    
      
    }
  }
  fDistMapIsSet = true;
  return true;
}






// *****************************************************************************
bool TTPCEmpDistCorrector::GetDistortionCorrectedPoint(const TVector3& recPos, TVector3 &corrPos) {
  if(fFourHits){
    //Get the four closest laser hits                                                                                                                                                                   
    vector<pair<TVector3,TLorentzVector> > closestFourHits;
    if (! GetClosestFourHits(recPos, closestFourHits) ) return false;

    //Nos correct the point
    if (! CorrectWeightedAverage(recPos, closestFourHits, corrPos) ) return false;
  }
 else{
   //Get to what laser hit this rec pos is closest to (in 2D)
  pair<TVector3,TLorentzVector> closestHit;
  if (! GetClosestHit(recPos, closestHit) ) return false;

  //We do not want to apply the distortions in the case of single point fit, we want to 
  //"correct" the distortion, so that means we will correct with the negative distortions
  TVector3 fullDriftDistCorr = -1 * closestHit.second.Vect();
  if (! Correct(recPos, fullDriftDistCorr, corrPos) ) return false;
 }
  return true;
}

// *****************************************************************************
bool TTPCEmpDistCorrector::ApplyDistOnPadPoint(const TVector3& trackPos, TLorentzVector &padPlanePoint) {
  if(fFourHits){
    //Get the four closest laser hits
    vector<pair<TVector3,TLorentzVector> > closestFourHits;
    if (! GetClosestFourHits(trackPos, closestFourHits) ) return false;

    //Apply the distortions                                                                                                                                                                            
    if (! ApplyWeightedAverage(trackPos, closestFourHits, padPlanePoint) ) return false;

  }
  else{
  //Get to what laser hit this rec pos is closest to (in 2D)
  pair<TVector3,TLorentzVector> closestHit;
  if (! GetClosestHit(trackPos, closestHit) ) return false;

  //Apply the distortions
  TLorentzVector fullDriftDistortion = closestHit.second;
  if (! Apply(trackPos, fullDriftDistortion, padPlanePoint) ) return false;
  }
  return true;
}

// *****************************************************************************
bool TTPCEmpDistCorrector::GetClosestHit(const TVector3& recPos, pair<TVector3,TLorentzVector>& closestHit) {
  double closestDist = 99999 * unit::m;
  vector<pair<TVector3,TLorentzVector> >::iterator hitIt;
  for (hitIt = fLaserHitsAndDistortions.begin(); hitIt != fLaserHitsAndDistortions.end(); hitIt++) {
    TVector3 laserHit = hitIt->first;
    //Check that the laser hit is on the same side of the central cathode (x=0) as recPos
    if (((laserHit.X() > 0)&&(recPos.X() > 0))||
	((laserHit.X() < 0)&&(recPos.X() < 0))) { 
      //If this laser hit is the closest (in z-y plane) to the recPos
      if (sqrt((recPos.Z() - laserHit.Z())*(recPos.Z() - laserHit.Z()) + 
	       (recPos.Y() - laserHit.Y())*(recPos.Y() - laserHit.Y())) < closestDist) {
	//Then save this laser-pad hit pair
	closestDist = sqrt((recPos.Z() - laserHit.Z())*(recPos.Z() - laserHit.Z()) + 
			   (recPos.Y() - laserHit.Y())*(recPos.Y() - laserHit.Y()));
	closestHit = (*hitIt);
      }
    }
  }
  return closestDist < 99998 * unit::m;
}

//Get the four closest laser hits in 2D
// *****************************************************************************
bool TTPCEmpDistCorrector::GetClosestFourHits(const TVector3& recPos, vector<pair<TVector3,TLorentzVector> >& closestFourHits) {
  double closestDist = 99999 * unit::m;
  vector<pair<TVector3,TLorentzVector> >::iterator hitIt;
  pair<TVector3,TLorentzVector> closestHit;

  //first find the closest laser hit
  for (hitIt = fLaserHitsAndDistortions.begin(); hitIt != fLaserHitsAndDistortions.end(); hitIt++) {
    TVector3 laserHit = hitIt->first;
    //Check that the laser hit is on the same side of the central cathode (x=0) as recPos                                                                                                               
    if (((laserHit.X() > 0)&&(recPos.X() > 0))||
        ((laserHit.X() < 0)&&(recPos.X() < 0))) {
      //If this laser hit is the closest (in z-y plane) to the recPos                                                                                                                                   
      if (sqrt((recPos.Z() - laserHit.Z())*(recPos.Z() - laserHit.Z()) +
               (recPos.Y() - laserHit.Y())*(recPos.Y() - laserHit.Y())) < closestDist) {
        //Then save this laser-pad hit pair                                                                                                                                                            
        closestDist = sqrt((recPos.Z() - laserHit.Z())*(recPos.Z() - laserHit.Z()) +
                           (recPos.Y() - laserHit.Y())*(recPos.Y() - laserHit.Y()));
        closestHit = (*hitIt);
      }
    }
  }
  //save the closest hit
  closestFourHits.push_back(closestHit);

  closestDist = 99999 * unit::m;
  //now find the other three neighbors
  for(int i =0 ; i < 3; i++){
    for (hitIt = fLaserHitsAndDistortions.begin(); hitIt != fLaserHitsAndDistortions.end(); hitIt++) {
      TVector3 laserHit = hitIt->first;
      //Check that the laser hit is on the same side of the central cathode (x=0) as recPos                                                                                                               
      if (((laserHit.X() > 0)&&(recPos.X() > 0))||
	  ((laserHit.X() < 0)&&(recPos.X() < 0))) {
	//If this laser hit is the closest (in z-y plane) to the recPos 
	if (sqrt((recPos.Z() - laserHit.Z())*(recPos.Z() - laserHit.Z()) +
		 (recPos.Y() - laserHit.Y())*(recPos.Y() - laserHit.Y())) < closestDist) {
	  //and it is not the same we previously found 
	  vector<pair<TVector3,TLorentzVector> >::iterator it;
	  bool diff = true;
	  for (it = closestFourHits.begin(); it != closestFourHits.end(); it++){
	    TVector3 posIt = it->first;
	    if((posIt.Y()== laserHit.Y()) && (posIt.Z()== laserHit.Z())) diff = false;
	  }
	  if(diff){
	    //Then save this laser-pad hit pair                                                                                                                                                      
	    closestDist = sqrt((recPos.Z() - laserHit.Z())*(recPos.Z() - laserHit.Z()) +
			       (recPos.Y() - laserHit.Y())*(recPos.Y() - laserHit.Y()));
	    closestHit = (*hitIt);
	    
	  }
	}
      }
    }
    closestFourHits.push_back(closestHit);
  }
  
  return (closestDist < 99998 * unit::m);
}





// *****************************************************************************
bool TTPCEmpDistCorrector::Apply(const TVector3& trackPos, const TLorentzVector& fullDriftDist, TLorentzVector& padPlanePoint) {
  double xDistToPadPlane = fGlobalXOfPadPlane - fabs(trackPos.X()); 
  double yDist = xDistToPadPlane / fMaxDriftInX * fullDriftDist.Y();
  double zDist = xDistToPadPlane / fMaxDriftInX * fullDriftDist.Z();
  double tDist = xDistToPadPlane / fMaxDriftInX * fullDriftDist.T();
  padPlanePoint.SetY(padPlanePoint.Y() + yDist); 
  padPlanePoint.SetZ(padPlanePoint.Z() + zDist); 
  padPlanePoint.SetT(padPlanePoint.T() + tDist); 
  //So far I see no reason to fail
  return true;
}



// *****************************************************************************
bool TTPCEmpDistCorrector::Correct(const TVector3& recPos, const TVector3& fullDriftDistCorr, TVector3& corrPos) {
  double xDistToPadPlane = fGlobalXOfPadPlane - fabs(recPos.X()); 
  double xCorr = xDistToPadPlane / fMaxDriftInX * fullDriftDistCorr.X();
  double yCorr = xDistToPadPlane / fMaxDriftInX * fullDriftDistCorr.Y();
  double zCorr = xDistToPadPlane / fMaxDriftInX * fullDriftDistCorr.Z();
  corrPos.SetX(recPos.X() + xCorr); 
  corrPos.SetY(recPos.Y() + yCorr); 
  corrPos.SetZ(recPos.Z() + zCorr); 
  //So far I see no reason to fail
  return true;
}

// *****************************************************************************
bool TTPCEmpDistCorrector::CorrectWeightedAverage(const TVector3& recPos,vector<pair<TVector3,TLorentzVector> >& closestFourHits , TVector3& corrPos) {

  //first calculate the distance in 2D to the four neighbors
  double d[4];
  vector<pair<TVector3,TLorentzVector> >::iterator it = closestFourHits.begin(); 
  for (int i =0; i<4; i++){
    TVector3 posIt = it->first;
    d[i]= sqrt((recPos.Z() - posIt.Z())*(recPos.Z() - posIt.Z()) + (recPos.Y() - posIt.Y())*(recPos.Y() - posIt.Y()));
    it++;
  }
  //now calculate the weights, inversely proportional to the distance
  double k;
  k=(d[0]*d[1]*d[2]*d[3])/(d[0]*d[1]*d[2]+d[0]*d[1]*d[3]+d[0]*d[2]*d[3]+d[1]*d[2]*d[3]) ;

  double w[4];
  for (int i=0; i<4; i++){
    w[i]=k/d[i];
  } 
  //calculate the full drift distortion
  double xDist =0;
  double yDist =0;
  double zDist =0;
  it = closestFourHits.begin();
  for (int i =0; i<4; i++){
    TVector3 fullDriftDistCorr = -1 * it->second.Vect();
    xDist += w[i]*fullDriftDistCorr.X();
    yDist += w[i]*fullDriftDistCorr.Y();
    zDist += w[i]*fullDriftDistCorr.Z();
    it++;
  }

  //finally calculate the correction linearly
  double xDistToPadPlane = fGlobalXOfPadPlane - fabs(recPos.X());
  double xCorr = xDistToPadPlane / fMaxDriftInX * xDist;
  double yCorr = xDistToPadPlane / fMaxDriftInX * yDist;
  double zCorr = xDistToPadPlane / fMaxDriftInX * zDist;
  corrPos.SetX(recPos.X() + xCorr);
  corrPos.SetY(recPos.Y() + yCorr);
  corrPos.SetZ(recPos.Z() + zCorr);

  return true;
}

// *****************************************************************************
bool TTPCEmpDistCorrector::ApplyWeightedAverage(const TVector3& trackPos,vector<pair<TVector3,TLorentzVector> >& closestFourHits , TLorentzVector& padPlanePoint) {
  //first calculate the distance in 2D to the four neighbors                                                                                                                                            
  double d[4];
  vector<pair<TVector3,TLorentzVector> >::iterator it = closestFourHits.begin();
  for (int i =0; i<4; i++){
    TVector3 posIt = it->first;
    d[i]= sqrt((trackPos.Z() - posIt.Z())*(trackPos.Z() - posIt.Z()) + (trackPos.Y() - posIt.Y())*(trackPos.Y() - posIt.Y()));
    it++;
  }
  //now calculate the weights, inversely proportional to the distance                                                                                                                                    
  double k;
  k=(d[0]*d[1]*d[2]*d[3])/(d[0]*d[1]*d[2]+d[0]*d[1]*d[3]+d[0]*d[2]*d[3]+d[1]*d[2]*d[3]) ;

  double w[4];
  for (int i=0; i<4; i++){
    w[i]=k/d[i];
  }
  //calculate the full drift distortion                                                                                                                                                                 
  double tFullDist =0;
  double yFullDist =0;
  double zFullDist =0;
  it = closestFourHits.begin();
  for (int i =0; i<4; i++){
    TLorentzVector fullDriftDist = it->second;
    tFullDist += w[i]*fullDriftDist.T();
    yFullDist += w[i]*fullDriftDist.Y();
    zFullDist += w[i]*fullDriftDist.Z();
    it++;
  }

  //finally calculate the correction linearly                                                                                                                                                           
  double xDistToPadPlane = fGlobalXOfPadPlane - fabs(trackPos.X());
  double tDist = xDistToPadPlane / fMaxDriftInX * tFullDist;
  double yDist = xDistToPadPlane / fMaxDriftInX * yFullDist;
  double zDist = xDistToPadPlane / fMaxDriftInX * zFullDist;
  padPlanePoint.SetT(padPlanePoint.T() + tDist);
  padPlanePoint.SetY(padPlanePoint.Y() + yDist);
  padPlanePoint.SetZ(padPlanePoint.Z() + zDist);

  return true;
}
