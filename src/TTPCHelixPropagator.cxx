#include "TTPCHelixPropagator.hxx"
#include "TTPCUtils.hxx"

/// The static member pointer to the singleton.
trex::TTPCHelixPropagator* trex::TTPCHelixPropagator::_helixPropagator = NULL;

//*****************************************************************************
trex::TTPCHelixPropagator& trex::helixPropagator(){

  return trex::TTPCHelixPropagator::Get();

}

//*****************************************************************************
trex::TTPCHelixPropagator& trex::TTPCHelixPropagator::Get(void) {
  if (!_helixPropagator) {
    _helixPropagator = new trex::TTPCHelixPropagator();
  }

  return *_helixPropagator;
}

//*****************************************************************************
trex::TTPCHelixPropagator::TTPCHelixPropagator(){
  Reset();
}


//*****************************************************************************
void trex::TTPCHelixPropagator::Reset(){
  fX0 = 0.0;
  fY0 = 0.0;
  fZ0 = 0.0;
  fDirX0 = 0.0;
  fDirY0 = 0.0;
  fDirZ0 = 0.0;
  fPhiQuad = 0.0;
  fPhiZY = 0.0;
  fRho0 = 0.0;
  fZc = 0.0;
  fYc = 0.0;
  fQuadrant = 0;
  fInitSense = 0;
}


//*****************************************************************************
void trex::TTPCHelixPropagator::CalculatePhisAndCenters(){
  // Use DirY and DirZ to get absolute phi angle wrt the ZY axis, fPhiZY.
  if (fInitSense > 0)
    fPhiZY = TMath::ATan2(fDirZ0, fDirY0);
  fPhiZY = TMath::ATan2(fDirZ0, fInitSense*fDirY0);
  if ( fQuadrant == 2){
    if (fInitSense < 0.0)
      fPhiZY = TMath::ATan2(-1.*fDirZ0, fDirY0);
    if (fInitSense > 0.0)
      fPhiZY = TMath::ATan2(fDirZ0, -1.*fDirY0);
    fPhiQuad = fPhiZY - TMath::PiOver2();
  } else if ( fQuadrant == 4){
    if (fInitSense < 0.0)
      fPhiZY = TMath::ATan2(fDirZ0, -1.*fDirY0);
    if (fInitSense > 0.0)
      fPhiZY = TMath::ATan2(-1.*fDirZ0, fDirY0);
    fPhiQuad = fPhiZY + TMath::PiOver2();
  } else if ( fQuadrant == 1){
    if (fInitSense < 0.0)
      fPhiZY = TMath::ATan2(fDirZ0, -1.*fDirY0);
    if (fInitSense > 0.0)
      fPhiZY = TMath::ATan2(-1.*fDirZ0, fDirY0);
    fPhiQuad = fPhiZY;
  } else if ( fQuadrant == 3){
    if (fInitSense < 0.0)
      fPhiZY = TMath::ATan2(-1.*fDirZ0, fDirY0);
    if (fInitSense > 0.0)
      fPhiZY = TMath::ATan2(fDirZ0, -1.*fDirY0);
    if (fPhiZY > 0.0)
      fPhiQuad = fPhiZY - TMath::Pi();
    else
      fPhiQuad = fPhiZY + TMath::Pi();
  } else {
    // TODO:: proper exception
    throw;
  }
  
  // Centre of the circle
  fYc = fY0 - TMath::Sin(fPhiZY)/fabs(fRho0);
  fZc = fZ0 - TMath::Cos(fPhiZY)/fabs(fRho0);
}


//*****************************************************************************
void trex::TTPCHelixPropagator::FindQuadrant(){
  if ( fFirstCluIsVert ){
    if ( (fInitSense * fRho0) > 0 ){
      fQuadrant = 2;
    } else {
      fQuadrant = 4;
    }
  } else {
    if ( (fInitSense * fRho0) > 0 ){
      fQuadrant = 3;
    } else {
      fQuadrant = 1;
    }
  }

}


//*****************************************************************************
void trex::TTPCHelixPropagator::ReloadHelixPosTanCurv(std::vector<double> Param){

  fX0    = Param[0];
  fY0    = Param[1];
  fZ0    = Param[2];
  fRho0  = Param[5];

  FindQuadrant();

  double sqrtslopes = TMath::Sqrt(1 + Param[3]*Param[3] + Param[4]*Param[4]);
  if (fQuadrant == 2 || fQuadrant == 4){
      
    fDirY0 = fInitSense * Param[4] / sqrtslopes;
    fDirZ0 = fInitSense / sqrtslopes;
  } else if (fQuadrant == 1 || fQuadrant == 3){
      
    fDirY0 = fInitSense / sqrtslopes;
    fDirZ0 = fInitSense * Param[4] / sqrtslopes;
  } else {
    // TODO:: proper exception
    throw;
  }
  fDirX0 = fInitSense * Param[3] / sqrtslopes;
    
  CalculatePhisAndCenters();

}


//*****************************************************************************
bool trex::TTPCHelixPropagator::InitHelixPosDirQoP(std::vector<double> Param, bool FirstCluIsVertical){
  fFirstCluIsVert = FirstCluIsVertical;

  if ( fFirstCluIsVert ){
    // For a vertical cluster, the Z direction cannot be zero !
    if ( ! Param[5] ){
      // TODO:: proper exception
      std::cout<<"Z dir 0 in vertical cluster"<<std::endl;
      throw;
    }
    if ( Param[5] > 0)
      fInitSense = 1;
    else
      fInitSense = -1;
  } else {
    // For a horizontal cluster, the Y direction cannot be zero !
    if ( ! Param[4] ){
      std::cout<<"Y dir 0 in horizontal cluster"<<std::endl;
      // TODO:: proper exception
      throw;
    }
    if ( Param[4] > 0)
      fInitSense = 1;
    else
      fInitSense = -1;
  }
  bool ok = ReloadHelixPosDirQoP(Param);
  return ok;
}

//*****************************************************************************
bool trex::TTPCHelixPropagator::ReloadHelixPosDirQoP(std::vector<double> Param){
  double Mag = TMath::Sqrt(Param[3]*Param[3] + Param[4]*Param[4] + Param[5]*Param[5]);
  fX0    = Param[0];
  fY0    = Param[1];
  fZ0    = Param[2];
  fDirX0 = Param[3]/Mag;
  fDirY0 = Param[4]/Mag;
  fDirZ0 = Param[5]/Mag;
  TVector3 Position(fX0, fY0, fZ0);
  TVector3 Direction(fDirX0, fDirY0, fDirZ0);
  double Momentum = fabs(1./Param[6]);
  double Charge = Param[6] / fabs(Param[6]);
  bool ok = TTPCUtils::MomentumAndCharge_to_Curvature(Position, Direction, Momentum, Charge, fRho0);

  if (!ok){
    return false;
  }
  
  // Restart from the Init point
  FindQuadrant();

  CalculatePhisAndCenters();

  return true;
}


//*****************************************************************************
int trex::TTPCHelixPropagator::GetQuadrant(){
  return fQuadrant;
}


//*****************************************************************************
int trex::TTPCHelixPropagator::GetSense(){
  return fInitSense;
}


//*****************************************************************************
bool trex::TTPCHelixPropagator::PropagateToHVCluster(trex::TTRExHVCluster& Cluster,double* length){

  /*  std::cout<<"Initial state at "<<fX0<<", "<<fY0<<", "<<fZ0<<std::endl;
  std::cout<<"Directions at "<<fDirX0<<", "<<fDirY0<<", "<<fDirZ0<<std::endl;
  std::cout<<"Propagating to cluster at"<<Cluster.X()<<", "<<Cluster.Y()<<", "<<Cluster.Z()<<std::endl;
  */
  // Start by checking the orientation of the cluster.
  // If the mode doesn't match the quadrant, change the PhiQuad accordingly.
  // If we are already in the right quadrant, don't do anything.
  unsigned int newQuadrant = fQuadrant;

  if ( Cluster.IsVertical() && !(fQuadrant == 2 || fQuadrant == 4)){
    // Which quadrant ?
    if (Cluster.Y() > fYc)
      newQuadrant = 2;
    else
      newQuadrant = 4;
    // Change the PhiQuad to match the new axis of reference.
    if ( newQuadrant == 2 && fQuadrant == 1 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    } else if ( newQuadrant == 2 && fQuadrant == 3 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 4 && fQuadrant == 1 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 4 && fQuadrant == 3 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    }
  } else if ( (!Cluster.IsVertical()) && !(fQuadrant == 1 || fQuadrant == 3)){
    // Which quadrant ?
    if (Cluster.Z() > fZc)
      newQuadrant = 1;
    else
      newQuadrant = 3;
    // Change the PhiQuad to match the new axis of reference.
    if ( newQuadrant == 1 && fQuadrant == 2 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 1 && fQuadrant == 4 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    } else if ( newQuadrant == 3 && fQuadrant == 4 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 3 && fQuadrant == 2 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    }
  }
  fQuadrant = newQuadrant;
    
  // Then calculate the DeltaZ or DeltaY depending on the cluster's orientation.
  double DeltaYorZ = 0.0;
  if (fQuadrant == 1){
    DeltaYorZ = Cluster.Y() - fY0;
  } else if (fQuadrant == 3){
    DeltaYorZ = fY0 - Cluster.Y();
  } else if (fQuadrant == 2){
    DeltaYorZ = fZ0 - Cluster.Z();
  } else if (fQuadrant == 4){
    DeltaYorZ = Cluster.Z() - fZ0;
  }

  double sinPhi = TMath::Sin(fPhiQuad) + DeltaYorZ * fabs(fRho0);
  
  if( sinPhi < -1. )
    sinPhi = -1.;
  else if( sinPhi > 1. ) 
    sinPhi = 1.;

  double newPhiQuad = TMath::ASin(sinPhi); 

  // Calculate the delta phi.
  double deltaPhi = newPhiQuad - fPhiQuad;
  if( (DeltaYorZ > 1e-3) && !deltaPhi)
    return false;
  // Then calculate the coordinate of the new point.
  double DeltaZorY = (TMath::Cos(newPhiQuad)-TMath::Cos(fPhiQuad) ) / fabs(fRho0);
  //  std::cout<<" ===> DeltaZorY = "<<DeltaZorY<<std::endl;
  // Save the new point as the new helix state.
  fPhiZY += deltaPhi;
  if ( fRho0 > 0.0 ){ // Negatively charged track
    fDirY0 = -TMath::Cos(fPhiZY);
    fDirZ0 = TMath::Sin(fPhiZY);
  } else {              // Positively charged track
    fDirY0 = TMath::Cos(fPhiZY);
    fDirZ0 = -TMath::Sin(fPhiZY);
  }
  fPhiQuad = newPhiQuad;

  double Renorm = TMath::Sqrt((1 - fDirX0*fDirX0)/(fDirY0*fDirY0 + fDirZ0*fDirZ0));
  fDirY0 *= Renorm;
  fDirZ0 *= Renorm;

  double DeltaX = -1. * (fDirX0 / TMath::Sqrt(fDirY0*fDirY0 + fDirZ0*fDirZ0)) * (deltaPhi/fRho0);

  if(length){
    *length=fabs(DeltaX)*TMath::Sqrt(1.+(fDirY0*fDirY0 + fDirZ0*fDirZ0)/fDirX0/fDirX0);
  }

  //  std::cout<<" ===> fRho0 = "<<fRho0<<std::endl;
  //std::cout<<" ===> DeltaX = "<<DeltaX<<std::endl;

  if ( fQuadrant == 1 ){
    fX0 = fX0 + DeltaX;
    fY0 = Cluster.Y();
    fZ0 = fZ0 + DeltaZorY;
  } else if ( fQuadrant == 2 ){
    fX0 = fX0 + DeltaX;
    fY0 = fY0 + DeltaZorY;
    fZ0 = Cluster.Z();
  } else if ( fQuadrant == 3 ){
    fX0 = fX0 + DeltaX;
    fY0 = Cluster.Y();
    fZ0 = fZ0 - DeltaZorY;
  } else if ( fQuadrant == 4 ){
    fX0 = fX0 + DeltaX;
    fY0 = fY0 - DeltaZorY;
    fZ0 = Cluster.Z();
  }

  //  std::cout<<"Got state at "<<fX0<<", "<<fY0<<", "<<fZ0<<std::endl;

  return true;
}

//*****************************************************************************
bool trex::TTPCHelixPropagator::FullPropagateToHVCluster(trex::TTRExHVCluster& Cluster,double* length){

  std::cout<<"FULL propagation"<<std::endl;

  // Start by checking the orientation of the cluster.
  // If the mode doesn't match the quadrant, change the PhiQuad accordingly.
  // If we are already in the right quadrant, don't do anything.
  unsigned int newQuadrant = fQuadrant;

  std::cout<<"Initial state at "<<fX0<<", "<<fY0<<", "<<fZ0<<std::endl;
  std::cout<<"Directions at "<<fDirX0<<", "<<fDirY0<<", "<<fDirZ0<<std::endl;
  std::cout<<"Propagating to cluster at"<<Cluster.X()<<", "<<Cluster.Y()<<", "<<Cluster.Z()<<std::endl;
  std::cout<<"This cluster is "<<(Cluster.IsVertical()?"Vertical":"Horizontal")<<std::endl;
  
  if ( Cluster.IsVertical() && !(fQuadrant == 2 || fQuadrant == 4)){
    // Which quadrant ?
    if (Cluster.Y() > fYc)
      newQuadrant = 2;
    else
      newQuadrant = 4;
    // Change the PhiQuad to match the new axis of reference.
    if ( newQuadrant == 2 && fQuadrant == 1 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    } else if ( newQuadrant == 2 && fQuadrant == 3 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 4 && fQuadrant == 1 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 4 && fQuadrant == 3 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    }
  } else if ( (!Cluster.IsVertical()) && !(fQuadrant == 1 || fQuadrant == 3)){
    // Which quadrant ?
    if (Cluster.Z() > fZc)
      newQuadrant = 1;
    else
      newQuadrant = 3;
    // Change the PhiQuad to match the new axis of reference.
    if ( newQuadrant == 1 && fQuadrant == 2 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 1 && fQuadrant == 4 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    } else if ( newQuadrant == 3 && fQuadrant == 4 ){
      fPhiQuad = fPhiQuad + TMath::Pi()/2.0;
    } else if ( newQuadrant == 3 && fQuadrant == 2 ){
      fPhiQuad = fPhiQuad - TMath::Pi()/2.0;
    }
  }
  fQuadrant = newQuadrant;


  //Calc new phi angle
  //-If cluster is outside of helix then put it at an extremum
  //-First use position of cluster in best-defined direction to get angle with 2-fold degeneracy
  //-Then use the other cluster coordinate (y or z) to resolve degeneracy
  double newPhi=0.;

  if(Cluster.IsVertical()){
    if(Cluster.Z()>fZc+fabs(1./fRho0)){
      newPhi=0.;
    }
    else if(Cluster.Z()<fZc-fabs(1./fRho0)){
      newPhi=TMath::Pi();
    }
    else{
      newPhi=TMath::ACos(fabs(fRho0)*(Cluster.Z()-fZc));
      if(Cluster.Y()<fYc){
	newPhi=2.*TMath::Pi()-newPhi;
      }
    }
  }
  else{    
    if(Cluster.Y()>fYc+fabs(1./fRho0)){
      newPhi=0.5*TMath::Pi();
    }
    else if(Cluster.Y()<fYc-fabs(1./fRho0)){
      newPhi=1.5*TMath::Pi();
    }
    else{
      newPhi=TMath::ASin(fabs(fRho0)*(Cluster.Y()-fYc));
      if(Cluster.Z()<fZc){
	newPhi=TMath::Pi()-newPhi;
      }
    }
  }
  double deltaPhi=newPhi-fPhiZY;
  if(deltaPhi<0.){
    deltaPhi+=2*TMath::Pi();
  }

  //Now figure out how many rotations we should go through to reach the cluster...
  
  int nRot=floor( (fX0-Cluster.X())*(TMath::Sqrt(fDirY0*fDirY0 + fDirZ0*fDirZ0)/fDirX0*fRho0/2./TMath::Pi())-deltaPhi/2./TMath::Pi()+0.5);

  double deltaX=-fDirX0/TMath::Sqrt(fDirY0*fDirY0 + fDirZ0*fDirZ0)*(deltaPhi+2.*TMath::Pi()*nRot)/fRho0;
  std::cout<<"rho0="<<fRho0<<", nRot="<<nRot<<", deltaPhi="<<deltaPhi<<", deltaX="<<deltaX<<std::endl;

  fPhiZY = newPhi;
  if ( fRho0 > 0.0 ){ // Negatively charged track
    fDirY0 = -TMath::Cos(fPhiZY);
    fDirZ0 = TMath::Sin(fPhiZY);
  } else {              // Positively charged track
    fDirY0 = TMath::Cos(fPhiZY);
    fDirZ0 = -TMath::Sin(fPhiZY);
  }

  double Renorm = TMath::Sqrt((1 - fDirX0*fDirX0)/(fDirY0*fDirY0 + fDirZ0*fDirZ0));
  fDirY0 *= Renorm;
  fDirZ0 *= Renorm;

  fZ0=fZc+cos(fPhiZY)/fabs(fRho0);
  fY0=fYc+sin(fPhiZY)/fabs(fRho0);
  fX0=fX0+deltaX;

  std::cout<<"Got state at "<<fX0<<", "<<fY0<<", "<<fZ0<<std::endl;

  //MDH TODO: Put some couts here to make sure that the phiZY and centers are not changed
  //on reinitializing, and that the final X,Y,Z are close to the cluster.

  FindQuadrant();

  CalculatePhisAndCenters();

  if(length){
    *length=deltaX*TMath::Sqrt(1.+(fDirY0*fDirY0 + fDirZ0*fDirZ0)/fDirX0/fDirX0);
  }

  return true;
}


//*****************************************************************************
void trex::TTPCHelixPropagator::GetHelixPosDirCurv(std::vector<double>& Result){

  Result.resize(7);
  Result[0] = fX0;
  Result[1] = fY0;
  Result[2] = fZ0;
  Result[3] = fDirX0;
  Result[4] = fDirY0;
  Result[5] = fDirZ0;
  Result[6] = fRho0;
}


//*****************************************************************************
void trex::TTPCHelixPropagator::GetHelixPosTanCurv(std::vector<double>& Result){

  Result.resize(7);
  Result[0] = fX0;
  Result[1] = fY0;
  Result[2] = fZ0;
  if ( fQuadrant == 1  || fQuadrant == 3){
    Result[3] = fDirX0 / fDirY0;
    Result[4] = fDirZ0 / fDirY0;
  } else {
    Result[3] = fDirX0 / fDirZ0;
    Result[4] = fDirY0 / fDirZ0;
  }
  Result[5] = fRho0;
}



//*****************************************************************************
void trex::TTPCHelixPropagator::GetHelixPosDirQoP(std::vector<double>& Result){
  Result[0] = fX0;
  Result[1] = fY0;
  Result[2] = fZ0;
  Result[3] = fDirX0;
  Result[4] = fDirY0;
  Result[5] = fDirZ0;

  TVector3 Position(Result[0], Result[1], Result[2]);
  TVector3 Direction(Result[3], Result[4], Result[5]);
  double p, q, qoverp;
  // RecPack state store Q over P so we need to convert back to the curvature.
  if (TTPCUtils::Curvature_to_MomentumAndCharge(Position, Direction, fRho0, p, q)){
    qoverp = q/p;
  } else {
    qoverp = 0.0;
  }
  Result[6] = qoverp;
}



//*****************************************************************************
void trex::TTPCHelixPropagator::PosTanCurvToPosDirQoP(std::vector<double> &ptcVect, std::vector<double> &pdqpVect){

  pdqpVect.resize(7);

  ReloadHelixPosTanCurv(ptcVect);
  pdqpVect[0] = fX0;
  pdqpVect[1] = fY0;
  pdqpVect[2] = fZ0;
  pdqpVect[3] = fDirX0;
  pdqpVect[4] = fDirY0;
  pdqpVect[5] = fDirZ0;

  TVector3 Position(pdqpVect[0], pdqpVect[1], pdqpVect[2]);
  TVector3 Direction(pdqpVect[3], pdqpVect[4], pdqpVect[5]);
  double p, q, qoverp;
  // RecPack state store Q over P so we need to convert back to the curvature.
  if (TTPCUtils::Curvature_to_MomentumAndCharge(Position, Direction, fRho0, p, q)){
    qoverp = q/p;
  } else {
    qoverp = 0.0;
  }
  pdqpVect[6] = qoverp;

  /*
  // Now the covariance matrix
  //EMatrix Jacobian = EMatrix(7,6,0);

  // Use variables for readibility
  double ux  = ptcVect[3];   // slope X
  double uyz = ptcVect[4];   // slope Y or Z depending on the first cluster's orientation
  double sqrtxyz = TMath::Sqrt( 1 + ux*ux + uyz*uyz );
  double xyz32 = TMath::Power( (1 + ux*ux + uyz*uyz), 1.5);
  double sqrtyz  = TMath::Sqrt( 1 + uyz*uyz );
  double yz32  = TMath::Power( (1 + uyz*uyz), 1.5 );

  // Little trick so I don't have to get the BField to calculate rho /(0.3 * B)
  double rho_03B = qoverp * sqrtxyz / sqrtyz;

  Jacobian[0][0] = 1.;
  Jacobian[1][1] = 1.;
  Jacobian[2][2] = 1.;
  Jacobian[3][3] = fInitSense * ( (1. / sqrtxyz) - (ux*ux / xyz32));
  Jacobian[3][4] = -1. * fInitSense *ux * uyz / xyz32;
  
  // We used slopeY !
  if ( fQuadrant == 1  || fQuadrant == 3){
    Jacobian[4][3] = -1. * fInitSense * ux / xyz32;
    Jacobian[4][4] = -1. * fInitSense * uyz / xyz32;
    Jacobian[5][3] = -1. * fInitSense * uyz * ux / xyz32;
    Jacobian[5][4] = fInitSense * ( (1. / sqrtxyz) - (uyz*uyz / xyz32));
  }
  // We used slopeZ !
  else {
    Jacobian[4][3] = -1. * fInitSense * uyz * ux / xyz32;
    Jacobian[4][4] = fInitSense * ( (1. / sqrtxyz) - (uyz*uyz / xyz32));
    Jacobian[5][3] = -1. * fInitSense * ux / xyz32;
    Jacobian[5][4] = -1. * fInitSense * uyz / xyz32;
  }
  Jacobian[6][3] = rho_03B * sqrtyz * ux / xyz32;
  Jacobian[6][4] = rho_03B * uyz * ( (yz32 / sqrtxyz) - (sqrtyz / xyz32) );
  // Another shortcut
  Jacobian[6][5] = qoverp / fRho0;
    
  pdqpCova = (Jacobian * ptcCova) * Jacobian.T(); 
  */
}

