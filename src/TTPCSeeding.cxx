#include "TTPCSeeding.hxx" 
#include <TVector.h> 

#include "TTRExHVCluster.hxx" 
#include "TTPCUtils.hxx" 
#include "TTPCHelixPropagator.hxx"


//*****************************************************************************
trex::TTPCSeeding::TTPCSeeding( ){
//*****************************************************************************
  fChi2max = 1000;

}


//*****************************************************************************
void trex::TTPCSeeding::Process(trex::TTRExPattern& Pattern){
//*****************************************************************************
  // return false when no seed was found for any path ???
  std::vector<trex::TTRExPath>& Paths = Pattern.GetPaths();

  for (auto pth = Paths.begin(); pth != Paths.end(); pth++) {
    trex::TTRExPath& path = *pth;
    // If there is already a seed result, don't seed again.
    if ( path.HasChi2Fit() || 
	 path.HasRunFit()){
      continue;
    }
    // TODO: Reset all the seed variables for safety
    FindSeed( path);
  }
}


//*****************************************************************************
void trex::TTPCSeeding::FindSeed(trex::TTRExPath& thePath){
//*****************************************************************************

  std::cout<<"Finding seed for a path..."<<std::endl;

  std::vector<TTRExHVCluster*>& HVclu = thePath.GetClusters();
  if( !HVclu.size() ) {
    std::cout<<"  No clusters!"<<std::endl;
    return;
  }
  PrepareSeeding(HVclu);
  // Do we have enough valid clusters to get a decent seed ?
  int NbValidPlanes = 0;
  for (auto Hit = HVclu.begin(); Hit != HVclu.end(); Hit++) {
    trex::TTRExHVCluster& HV = (**Hit);
    if( HV.isOkForSeed() ) NbValidPlanes++;
  }
  if (NbValidPlanes < 4 ){
    // TODO: Somehow save the fact that we tried to make a seed but failed.
  std::cout<<"  Not enough valid planes!"<<std::endl;
    thePath.SetHasRunFit(true);
    return;
  }
  std::vector<double> RiemannHelix;
  double RiemannError = Riemann(HVclu, RiemannHelix);
  std::vector<double> R2Helix;
  double R2Error = R2(HVclu, R2Helix);
  bool RiemannIsBad = ( std::isnan(RiemannError) || RiemannError > 1.e6); 
  bool R2IsBad = ( std::isnan(R2Error) || R2Error > 1.e6);
  std::cout<<"  Riemann is bad="<<RiemannIsBad<<std::endl; 
  std::cout<<"  R2 is bad="<<R2IsBad<<std::endl; 
  if( RiemannIsBad && R2IsBad ) {
    // TODO: Save "seeding failed" in path ?
    thePath.SetHasRunFit(true);
    return;
  }

  std::vector<double> frontSeedState;
  std::vector<double> backSeedState;
  if( RiemannError > R2Error || (RiemannIsBad && R2IsBad)) {
    frontSeedState = R2Helix;
  }
  else {
    frontSeedState = RiemannHelix;
  }

  if( ! IsResultValid(HVclu, frontSeedState) ) {
    // TODO: Save "seeding failed" in path ?
    std::cout<<"  Result for front seed state is bad"<<std::endl; 
    thePath.SetHasRunFit(true);
    return;
  }

  // Calculate useful quantities for the following steps like the path length and the last state.
  // Use a copy of the seed state
  std::vector<double> propagState = frontSeedState;

  // For the final propagation in particular, the seed state may not be defined
  // at the first cluster but instead at the second or third if the first clusters
  // are not selected as good enough for the seeding.
  // So just propagate to the shortest distance.

  bool firstcluvertical=(*(HVclu.begin()))->IsVertical();
  trex::TTPCHelixPropagator& hp=trex::helixPropagator();
  hp.InitHelixPosDirQoP(propagState,firstcluvertical);
  for (auto Hit = HVclu.begin(); Hit != HVclu.end(); Hit++) {
    trex::TTRExHVCluster& Cluster = **Hit;

    if (!hp.PropagateToHVCluster(Cluster)){
      continue;
    }
    
    if (*(Hit) == (*(HVclu.begin()))){
      frontSeedState = propagState;
    }
  }
  
  // TODO? inverse front and back states when the length is negative ???
  // Or should we instead change the propagation sign ?
  // Not sure what's the best option.
  backSeedState = propagState;
  
  thePath.SaveSeedStates(frontSeedState, backSeedState);
  std::cout<<"  Saving seed states!"<<RiemannIsBad<<std::endl; 
}


//*****************************************************************************
// Calculate some quantities that are used by multiple seeding algorithm
// to speed things up.
void trex::TTPCSeeding::PrepareClustersForSeeding( std::vector<trex::TTRExHVCluster*>& HVclu ){
  //*****************************************************************************
  // Select the clusters that are good enough for the sseding.
  //int NMaxPeaks = 0;
  //int NSaturation = 0;  
  int NSelVert = 0;
  int NSelHori = 0;
 
  fNbOrientChange = 0;

  bool PrevIsVert = true;
  bool FirstClu = true;
  
  for (auto tmpClu = HVclu.begin(); tmpClu != HVclu.end(); tmpClu++) {
    trex::TTRExHVCluster& Clu = **tmpClu;

    Clu.SetOkForSeed(true);  // Make sure that we start with fresh sample.

    bool ThisIsVert = Clu.IsVertical();
    if ( PrevIsVert != ThisIsVert && !FirstClu ){
      fNbOrientChange++;
      PrevIsVert = ThisIsVert;
    }

    //PD DO WE NEED THIS CRITERIA? - LEAVING IT IN FOR NOW
    // Select clusters allowing only 2 changes of orientation, not more.
    if (fNbOrientChange > 2){
      Clu.SetOkForSeed(false);
      continue;
    }

    // Cluster selected !
    if( Clu.IsVertical() ){
      NSelVert++;
    } else {
      NSelHori++;
    }

    FirstClu = false;
    PrevIsVert = Clu.IsVertical();
  }


  // We removed the clusters starting at the 3rd change or orientation.
  // So record only the first 2 changes.
  if (fNbOrientChange > 2){
    fNbOrientChange = 2;
  }
  
}


//*****************************************************************************
// Calculate some quantities that are used by multiple seeding algorithms
// to speed things up.

void trex::TTPCSeeding::PrepareSeeding( std::vector<trex::TTRExHVCluster*>& HVclu ){
//*****************************************************************************
  fZfirst = 900000.;
  fZlast = -900000.;
  fYfirst = 0.0;
  fYlast = 0.0; 
  fXfirst = 0.0;
  fXlast = 0.0; 
  fYmid = 0.0;
  fZmid = 0.0;
  double zclu,yclu,xclu;

  PrepareClustersForSeeding(HVclu);
  
  bool FoundFirst = false;
  int NbSelectedClu = 0;

  for (auto hit = HVclu.begin(); hit != HVclu.end(); hit++) {
    trex::TTRExHVCluster& clu = **hit;
    if( !clu.isOkForSeed() ) continue; 

    xclu = clu.X();
    zclu = clu.Z();
    yclu = clu.Y();

    fXlast = xclu;
    fYlast = yclu;
    fZlast = zclu;

    if( !FoundFirst){
      fXfirst = xclu;
      fYfirst = yclu;
      fZfirst = zclu;
      FoundFirst = true;
    }

    NbSelectedClu++;
  }

  int TargetClu = int(double(NbSelectedClu) / 2.);
  NbSelectedClu = 0;

  for (auto hitit = HVclu.begin(); hitit != HVclu.end(); hitit++){
    trex::TTRExHVCluster& clu = **hitit;
    if( !clu.isOkForSeed() ) continue;

    // Search for Zmid.
    if( TargetClu == NbSelectedClu )  {
      fZmid = clu.Z();
      fYmid = clu.Y();
      break;
    }
    NbSelectedClu++;
  }
}


//*****************************************************************************
// R2 (3 point seeding method)
//*****************************************************************************
double trex::TTPCSeeding::R2( std::vector<trex::TTRExHVCluster*>& HVclu, std::vector<double>& Helix){

  std::cout<<" R2 method: "<<std::endl;  
  double dy12 = fYfirst-fYmid; 
  double ay12 = fYfirst+fYmid; 
  double dy23 = fYmid-fYlast; 
  double ay23 = fYmid+fYlast; 

  double dz12 = fZfirst-fZmid; 
  double az12 = fZfirst+fZmid; 
  double dz23 = fZmid-fZlast; 
  double az23 = fZmid+fZlast; 

  double r12 = (dy12*ay12 + dz12*az12)/2.;
  double r23 = (dy23*ay23 + dz23*az23)/2.;

  double yl0 = (r12*dz23-r23*dz12)/(dy12*dz23-dy23*dz12);
  double zl0 = (r12*dy23-r23*dy12)/(dz12*dy23-dz23*dy12);

  std::cout<<"R2 yPts ("<<fYfirst<<", "<<fYmid<<", "<<fYlast<<")"<<std::endl;
  std::cout<<"R2 zPts ("<<fZfirst<<", "<<fZmid<<", "<<fZlast<<")"<<std::endl;
  std::cout<<"R2 (yl0,zl0) ("<<yl0<<", "<<zl0<<")"<<std::endl;

  double phiFirst = TMath::ATan2((fYfirst-yl0),(fZfirst-zl0)); 
  // Only for printout 
  double phiLast =  TMath::ATan2((fYlast-yl0),(fZlast-zl0)); 

  // Assuming the sense of the track is correct we have:
  // Charge = -Helicity 
  double Helicity = CalculateRhoSign(yl0, zl0);

  double r2 = (fYfirst-yl0)*(fYfirst-yl0)+(fZfirst-zl0)*(fZfirst-zl0); 

  if( std::isnan(yl0) || std::isnan(zl0) ) r2 = 1.e-20;

  double R2rho = 1./TMath::Sqrt(r2);
  
  if( R2rho < RHOMIN ) R2rho = RHOMIN;

  R2rho *= Helicity;
  double R2theta = phiFirst - (Helicity * TMath::PiOver2());

  // At this point just use the projection of the clusters
  // onto the center of the cathode to figure out the delta phi angle (less than 1pi ?)

  double R2x = fXfirst;
  double R2y = fYfirst;
  double R2z = fZfirst;

  double yzLength;

  if(fabs(phiLast-phiFirst)>1.E-6){
    yzLength=fabs((phiLast-phiFirst)/R2rho);
  }
  else{
    yzLength=std::max(fabs(fZlast-fZfirst),fabs(fYlast-fYfirst));
  }

  double xLength=fXlast-fXfirst;
  
  double renorm=1./TMath::Sqrt(yzLength*yzLength+xLength*xLength);

  Helix.push_back(R2x);
  Helix.push_back(R2y);
  Helix.push_back(R2z);
  Helix.push_back(xLength*renorm);
  Helix.push_back(TMath::Sin(R2theta)*yzLength*renorm);
  Helix.push_back(TMath::Cos(R2theta)*yzLength*renorm);

  std::cout<<"  (x,y,z)=("<<R2x<<","<<R2y<<","<<R2z<<")"<<std::endl;
  std::cout<<"  theta="<<R2theta<<", rho="<<R2rho<<std::endl;

  TVector3 Pos(Helix[0], Helix[1], Helix[2]);
  TVector3 Dir(Helix[3], Helix[4], Helix[5]);
  double p;
  double q;
  
  TTPCUtils::Curvature_to_MomentumAndCharge(Pos,Dir,R2rho,p,q);

  Helix.push_back(q/p);

  std::cout<<"Helix:"<<std::endl;
  for(int i=0;i<7;++i){
    std::cout<<"["<<i<<"]:"<<Helix[i]<<std::endl;
  }

  return FinalizeSeed(HVclu, R2rho, Helix);
}

//*****************************************************************************
// Riemann Fast Fit (Riemann sphere seeding method) 
//*****************************************************************************
double trex::TTPCSeeding::Riemann( std::vector<TTRExHVCluster*>& HVclu, std::vector<double>& Helix){
  double scale = 1000.; 

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  double rx = 0.0; 
  double ry = 0.0;
  double rz = 0.0;

  for (auto hitit = HVclu.begin(); hitit != HVclu.end(); hitit++) {
    trex::TTRExHVCluster& clu = **hitit;
    if( !clu.isOkForSeed() ) continue; 

    double yy = (clu.Y() - fYfirst)/scale;
    double zz = (clu.Z() - fZfirst)/scale;
    double r = TMath::Sqrt(zz*zz+yy*yy);
    double phi = TMath::ATan2(yy,zz);
    double xi = r/(1.+r*r)*TMath::Cos(phi);
    double yi = r/(1.+r*r)*TMath::Sin(phi);
    double zi = (r*r)/(1.+r*r);

    if( r > 0.0 ) {
      x.push_back(xi);
      y.push_back(yi);
      z.push_back(zi);
      rx += xi;
      ry += yi;
      rz += zi;
    }
  }

  rx /= (double)  x.size();
  ry /= (double)  y.size();
  rz /= (double)  z.size();

  TMatrix A(3,3);

  for(int i=0;i < 3; i++)
    for(int j=0;j < 3; j++)
      A[i][j] = 0.0;

  for(unsigned int i=0;i < x.size(); i++){
    x[i] -= rx;
    y[i] -= ry;
    z[i] -= rz;

    A[0][0] = A[0][0]+x[i]*x[i];
    A[0][1] = A[0][1]+x[i]*y[i];
    A[0][2] = A[0][2]+x[i]*z[i];

    A[1][0] = A[1][0]+y[i]*x[i];
    A[1][1] = A[1][1]+y[i]*y[i];
    A[1][2] = A[1][2]+y[i]*z[i];

    A[2][0] = A[2][0]+z[i]*x[i];
    A[2][1] = A[2][1]+z[i]*y[i];
    A[2][2] = A[2][2]+z[i]*z[i];
  }

  for(int i=0;i < 3; i++)
    for(int j=0;j < 3; j++)
      A[i][j] = A[i][j]/(double)x.size();

  TVector eigenValues;

  TMatrix M = A.EigenVectors(eigenValues);

  TVector n(3);

  n[0] = M[0][2];
  n[1] = M[1][2];
  n[2] = M[2][2];

  double c = -(n[0]*rx+n[1]*ry+n[2]*rz);

  double u0 = -n[0]/2./(c+n[2]);
  double v0 = -n[1]/2./(c+n[2]);
  double rho2 = (n[0]*n[0]+n[1]*n[1]-4.*c*(c+n[2]))/(4.*(c+n[2])*(c+n[2]));

  TVector val =  A*n;

  double Rierho = TMath::Sqrt(1./rho2)/scale;

  if( Rierho < RHOMIN ) Rierho = RHOMIN; 

  u0 *= scale;
  v0 *= scale;

  u0 += fZfirst;
  v0 += fYfirst;

  double phiFirst = TMath::ATan2((fYfirst-v0),(fZfirst-u0)); 
  // Only for printout 
  double phiLast =  TMath::ATan2((fYlast-v0),(fZlast-u0)); 

  // Assuming the sense of the track is correct we have:
  // Charge = -Helicity 
  double Helicity = CalculateRhoSign(v0, u0);

  Rierho = Rierho * Helicity;
  double Rietheta = phiFirst - (Helicity * TMath::PiOver2());

  std::cout << " Riemann:  zc  = " << u0 <<"      yc  = " << v0 << std::endl; 
  std::cout << " Riemann:  Helicity = "<<Helicity << std::endl; 
  std::cout << " Riemann:  phiFirst = "<<phiFirst <<"   phiLast = "<<phiLast << std::endl; 
  std::cout << " Riemann:  theta = "<<Rietheta<<"      rho = " << Rierho << std::endl;

  // At this point just use the projection of the clusters
  // onto the center of the cathode to figure out the delta phi angle (less than 1pi ?)

  double Riex = fXfirst;
  double Riey = fYfirst;
  double Riez = fZfirst;

  double yzLength;

  if(fabs(phiLast-phiFirst)>1.E-6){
    yzLength=fabs((phiLast-phiFirst)/Rierho);
  }
  else{
    yzLength=std::max(fabs(fZlast-fZfirst),fabs(fYlast-fYfirst));
  }

  double xLength=fXlast-fXfirst;
  
  double renorm=1./TMath::Sqrt(yzLength*yzLength+xLength*xLength);

  Helix.push_back(Riex);
  Helix.push_back(Riey);
  Helix.push_back(Riez);
  Helix.push_back(xLength*renorm);
  Helix.push_back(TMath::Sin(Rietheta)*yzLength*renorm);
  Helix.push_back(TMath::Cos(Rietheta)*yzLength*renorm);

  TVector3 Pos(Helix[0], Helix[1], Helix[2]);
  TVector3 Dir(Helix[3], Helix[4], Helix[5]);
  double p;
  double q;

  TTPCUtils::Curvature_to_MomentumAndCharge(Pos,Dir,Rierho,p,q);
  Helix.push_back(q/p);
  std::cout<<"Helix:"<<std::endl;
  for(int i=0;i<7;++i){
    std::cout<<"["<<i<<"]:"<<Helix[i]<<std::endl;
  }

  return FinalizeSeed(HVclu, Rierho, Helix);
}


//*****************************************************************************
double trex::TTPCSeeding::CalculateRhoSign(double Y0, double Z0) {
//*****************************************************************************
  // Just a short cut.
  bool largeCurv = false;
  double result;
  if (fNbOrientChange == 0 || fNbOrientChange == 1){
    if ( fabs(fZlast - fZfirst) > fabs(fYfirst - fYlast) ){
      // Is the track curves too much, using the center becomes problematic.
      if ( ( (fYlast - Y0)*(fYfirst - Y0) ) < 0.){
        largeCurv = true;
      } else {
        result = (fZlast - fZfirst) * (fYfirst - Y0);
        return ((result)/fabs(result));
      }
    } else {
      // Is the track curves too much, using the center becomes problematic.
      if ( ( (fZlast - Z0)*(fZfirst - Z0) ) < 0.){
        largeCurv = true;
      } else {
        result = (fYfirst - fYlast) * (fZfirst - Z0);
        return ((result)/fabs(result));
      }
    }
  }

  if (fNbOrientChange == 2 || largeCurv){
    if ( fabs(fZmid - fZfirst) > fabs(fYfirst - fYmid) ){
      result = (fZmid - fZfirst) * (fYfirst - Y0);
      return ((result)/fabs(result));
    } else {
      result = (fYfirst - fYmid) * (fZfirst - Z0);
      return ((result)/fabs(result));
    }
  }
  return 0.0;
}


//*****************************************************************************
double trex::TTPCSeeding::FinalizeSeed( std::vector<trex::TTRExHVCluster*>& HVclu, double rho, std::vector<double> &finalState){
//*****************************************************************************
  double rms = 0.0;
  int ntot=0;
  double ypred,zpred;
  double yclu,zclu;
  bool firstcluvertical=(*(HVclu.begin()))->IsVertical();

  if(firstcluvertical){
    if(fabs(finalState[5])<1.E-3){
      return 1.E99;
    }
  }
  else if(fabs(finalState[4])<1.E-3){
    return 1.E99;
  }

  std::vector<trex::TTRExHVCluster*>::iterator Hit = HVclu.begin();
  trex::TTRExHVCluster& Cluster = (**Hit);
  std::vector<double> propagState = finalState;
  double deltaPhi = 0.0;
  std::vector<double> prevState;
  double ClusterX1;
  double ClusterX2 = 0.0;
  ClusterX1 = Cluster.X();
  bool doDeltaPhi = true;
  double suspicious =0.0;
  trex::TTPCHelixPropagator& hp=trex::helixPropagator();
  hp.InitHelixPosDirQoP(propagState,firstcluvertical);
  for ( ; Hit != HVclu.end(); Hit++) {
    trex::TTRExHVCluster& Cluster = (**Hit);
    if( !Cluster.isOkForSeed() ) continue;
    prevState = propagState;
    double length = 0.0;
    if (!hp.PropagateToHVCluster(Cluster,&length)){
      continue;
    }
    hp.GetHelixPosDirQoP(propagState);

    ypred = propagState[1];
    zpred = propagState[2];
    yclu = Cluster.Y();
    zclu = Cluster.Z();
    double delta = TMath::Sqrt((yclu-ypred)*(yclu-ypred)+(zclu-zpred)*(zclu-zpred));
    if ( fabs(delta) > 100.){
      propagState = prevState;
      suspicious += 1.;
      continue;
    }
      
    if (doDeltaPhi)
      deltaPhi += length * rho;
    ClusterX2 = Cluster.X();
    if (deltaPhi > TMath::Pi() && doDeltaPhi){
      doDeltaPhi = false;
    }
    rms += delta*delta; 
    ntot++;
  }

  if ( (suspicious/(double)ntot) > 0.7){
    return 1.e99;
  }
  rms /= (double)ntot;

  double deltaX = ClusterX2-ClusterX1;
  finalState[3] = deltaX / TMath::Sqrt(deltaX*deltaX + (deltaPhi*deltaPhi)/(rho*rho));
  double Renorm = TMath::Sqrt((1 - finalState[3]*finalState[3])/(finalState[4]*finalState[4] + finalState[5]*finalState[5]));
  finalState[4] *= Renorm;
  finalState[5] *= Renorm;
  finalState[6] *= sqrt(1. - finalState[3]*finalState[3]);

  // Something went wrong, don't trust the result
  if ( ! (rms > 0.0) )
    rms = 1.e99;

  return rms;
}


//*****************************************************************************
bool trex::TTPCSeeding::IsResultValid( std::vector<trex::TTRExHVCluster*>& HVclu, std::vector<double>& Result){
//*****************************************************************************

  // For bad fits give tracks in X only and we can't fit those anyway.
  if (Result[4] == 0.0 && Result[5] == 0.0){
    return false;
  }

  // If the first cluster is vertical (horizontal), the z (y) direction cannot be zero.
  // TODO: maybe this constraint could be more stringent.
  trex::TTRExHVCluster& HV = **(HVclu.begin());
  if ( HV.IsVertical() && Result[5] == 0.0){
    return false;
  }
  if ( ( !HV.IsVertical()) && Result[4] == 0.0){
    return false;
  }

  return true;
}
