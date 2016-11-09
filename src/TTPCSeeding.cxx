#include "TTPCSeeding.hxx" 
#include "TTPCRecPackUtils.hxx"

#include <TVector.h> 

#include "TTPCCalibration.hxx" 
#include "TTPCHVCluster.hxx" 
#include "TTPCUtils.hxx" 
#include "TTPCDebug.hxx" 

#include <TOARuntimeParameters.hxx>
#include <TND280Event.hxx>
#include <TEventFolder.hxx>
#include <TrackingUtils.hxx>

//*****************************************************************************
ND::TTPCSeeding::TTPCSeeding( ){
//*****************************************************************************
  fChi2max = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.Seeding.MaxChi2");

  fExcludeClusterWithManyPeaks = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Seeding.ExcludeClusterWithManyPeaks");
  fExcludeSaturatedClusters = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Seeding.ExcludeSaturatedClusters");

  fUseTruthAsSeedResult = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.UseTruthAsSeedResults");
}


//*****************************************************************************
void ND::TTPCSeeding::Process(ND::THandle<ND::TTPCPattern> Pattern){
//*****************************************************************************
  // return false when no seed was found for any path ???
  std::vector< ND::THandle<ND::TTPCPath> > Paths = Pattern->GetPaths();
  
  for (std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin(); pth != Paths.end(); pth++) {
    ND::THandle<ND::TTPCPath> path = *pth;
    // If there is already a seed result, don't seed again.
    if ( path->CheckStatus(ND::TReconBase::kChi2Fit) || 
               path->CheckStatus(ND::TReconBase::kRan)){
      continue;
    }
    // TODO: Reset all the seed variables for safety
    FindSeed( path);
  }
}


//*****************************************************************************
void trex::TTPCSeeding::FindSeed(std::vector<trex::TTRExHVCluster>& HVclusters){
//*****************************************************************************
  //if ( ND::tpcDebug().Seeding(DB_INFO))
  //std::cout << " ========== Seeding on path Id: "<<thePath->GetId()<< std::endl; 

  //ND::THandle<ND::THitSelection> HVclu = thePath->GetHits();
  

  if( !HVclusters->size() ) return;

  // Little shortcut
  //fDriftVelocity = ND::tpcCalibration().GetDriftVelocity();

  PrepareSeeding(HVclusters);

    std::cout << " ====================== Seeding =====================" << std::endl; 
    //TTPCUtils::HVClustersPrintout(HVclu, ND::tpcDebug().SeededClusters(DB_VERBOSE));
    std::cout << " ----------------------------------------------------" << std::endl; 
  

  // Do we have enough valid clusters to get a decent seed ?
  int NbValidPlanes = 0;
  for (std::vector<trex::TTRExHVCluster>::iterator Hit = HVclusters->begin(); Hit != HVclusters->end(); Hit++) {
    trex::TTRExHVCluster HV = (*Hit);
    if( HV->isOkForSeed() ) NbValidPlanes++;
  }

  if (NbValidPlanes < 4 ){
    // TODO: Somehow save the fact that we tried to make a seed but failed.
    //thePath->SetStatus(ND::TReconBase::kRan);
    return;
  }
    
  std::cout << " --- Process Riemann seeding"<< std::endl; 
  std::vector<double> RiemannHelix;
  double RiemannError = Riemann(HVclusters, RiemannHelix);



  if ( ND::tpcDebug().Seeding(DB_INFO))
    std::cout << " --- Process R2 seeding"<< std::endl; 
  std::vector<double>  R2Helix;
  double R2Error = R2(HVclu, R2Helix);




  bool RiemannIsBad = ( isnan(RiemannError) || RiemannError > 1.e6); 
  bool R2IsBad = ( isnan(R2Error) || R2Error > 1.e6); 




  if( ND::tpcDebug().Seeding(DB_INFO) ){ 
    std::cout << " Riemann RMS " << RiemannError << std::endl; 
    std::cout << " 3 point RMS " << R2Error << std::endl; 
  }

  if( RiemannIsBad && R2IsBad ) {
    // TODO: Save "seeding failed" in path ?
    if( ND::tpcDebug().Seeding(DB_INFO) ){ 
      std::cout<<" ---> Both Riemann and R2 seeding algorithms FAILED ! "<<std::endl; 
    }
    thePath->SetStatus(ND::TReconBase::kRan);
    return;
  }

  State frontSeedState;
  State backSeedState;
  if( RiemannError > R2Error || (RiemannIsBad && R2IsBad)) {
    frontSeedState = R2Helix;
  }
  else {
    frontSeedState = RiemannHelix;
  }

  if( ! IsResultValid(HVclu, frontSeedState) ) {
    // TODO: Save "seeding failed" in path ?
    if( ND::tpcDebug().Seeding(DB_INFO) ){ 
      std::cout << " ---> The seed is invalid ! Don't save it."<< std::endl; 
    }
    thePath->SetStatus(ND::TReconBase::kRan);
    return;
  }

  if ( ND::tpcDebug().Seeding(DB_INFO)){
    std::cout<<" =========================================================================="<<std::endl;
    std::cout<<"  Path Id "<<thePath->GetId()<<std::endl;
    std::cout<<"  Riemann PosDirQoP: "<<RiemannHelix.vector()<<std::endl;
    std::cout<<"  R2 PosDirQoP: "<<R2Helix.vector()<<std::endl;
  }

  // We can use the truth as fake seed to check the impact of the seeding
  if (ND::tpcCalibration().IsMC() && 
    ( fUseTruthAsSeedResult || ND::tpcDebug().Seeding(DB_INFO) )) {
      State trueState = State(7);
      // This gets the truth or returns false. Includes check if input file is MC.
      bool ok = TTPCUtils::TrueStateNear3Dpoint(thePath->GetHits(), thePath->GetFirstPosition(), trueState);
      if (ok && fUseTruthAsSeedResult)
        frontSeedState = trueState;
      if (ok && ND::tpcDebug().Seeding(DB_INFO))
        std::cout<<"  Truth PosDirQoP: "<<trueState.vector()<<std::endl;
  }

  if ( ND::tpcDebug().Seeding(DB_INFO)){
    std::cout<<" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  }

  EVector vect = frontSeedState.vector();
  EMatrix cova = frontSeedState.matrix();
  // Just some values to have something.
  // TODO: These actually originate from TTPCt0 in tpcRecon and should
  // probably be reviewed
  cova[0][0] = 10000.0; //0.015;//0.15;
  cova[1][1] = 0.5; //0.015;//0.15;//0.5;
  cova[2][2] = 0.1; //0.015;//0.15;//0.5;
  cova[3][3] = 10000.0; //0.05;//0.5;
  cova[4][4] = 0.0001; //0.05;//0.5;//1.0;
  cova[5][5] = 0.0001; //0.05;//0.5;//1.0;
  cova[6][6] = 1.0e-18; //1e-10;//1.0e-9;//1.0e-7;
  frontSeedState.set_hv(HyperVector(vect,cova));

  // Calculate useful quantities for the following steps like the path length and the last state.
  // Use a copy of the seed state
  State propagState = frontSeedState;

  if( ND::tpcDebug().Seeding(DB_VVERBOSE) )
    std::cout<<" --> Get back seed state: Propagation"<<std::endl;
  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  // For the final propagation in particular, the seed state may not be defined
  // at the first cluster but instead at the second or third if the first clusters
  // are not selected as good enough for the seeding.
  // So just propagate to the shortest distance.
  ND::rpman("TREx").model_svc().model().intersector().set_length_sign(0);
  double pthLength = 0.0;
  for (ND::THitSelection::const_iterator Hit = HVclu->begin(); Hit != HVclu->end(); Hit++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = (*Hit);
    double localLength = 0.0;
    if (!TTPCRecPackUtils::PropagateToHVCluster(Cluster, propagState, localLength)){
      continue;
    }
    if( ND::tpcDebug().Seeding(DB_VVERBOSE) ){
      std::cout<<std::setw(20)<<" Clu at: "<<Cluster->GetPosition().Y()<<"   "<<Cluster->GetPosition().Z()<<std::endl;
      std::cout<<std::setw(20)<<" Propagation at: "<<propagState.vector()[1]<<"   "<<propagState.vector()[2]<<std::endl;
    }
    pthLength += localLength;

    if (Hit == HVclu->begin()){
      frontSeedState = propagState;
    }
  }

  // TODO? inverse front and back states when the length is negative ???
  // Or should we instead change the propagation sign ?
  // Not sure what's the best option.
  backSeedState = propagState;
  thePath->SetLength(pthLength);

  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);
  
  thePath->SaveSeedStates(frontSeedState, backSeedState);

}


//*****************************************************************************
// Calculate some quantities that are used by multiple seeding algorithm
// to speed things up.
void trex::TTPCSeeding::PrepareClustersForSeeding(std::vector<trex::TTRExHVCluster>& HVclu ){
//*****************************************************************************
  // Select the clusters that are good enough for the sseding.
  //int NMaxPeaks = 0;
  //int NSaturation = 0;  
  int NSelVert = 0;
  int NSelHori = 0;
 
  fNbOrientChange = 0;
  trex::TTRExHVCluster Clu;
  bool PrevIsVert = true;
  bool FirstClu = true;
  
  for (std::vector<trex::TTRExHVCluster>::iterator tmpClu = HVclu->begin(); tmpClu != HVclu->end(); tmpClu++) {
    Clu = *tmpClu;



    /*Clu->SetOkForSeed(true);  // Make sure that we start with fresh sample.

    if( Clu->GetMaxNPeaks() > 1 && fExcludeClusterWithManyPeaks )  {
      Clu->SetOkForSeed(false);
      NMaxPeaks++; 
      continue;
    }

    if( Clu->GetNSaturated() > 0  && fExcludeSaturatedClusters ) {
      Clu->SetOkForSeed(false);
      NSaturation++; 
      continue;
    }
    */



    bool ThisIsVert = Clu->IsVertical();
    if ( PrevIsVert != ThisIsVert && !FirstClu ){
      fNbOrientChange++;
      PrevIsVert = ThisIsVert;
    }


    //PD DO WE NEED THIS CRITERIA?
    // Select clusters allowing only 2 changes of orientation, not more.
    if (fNbOrientChange > 2){
      Clu->SetOkForSeed(false);
      continue;
    }



    // Cluster selected !
    if( Clu->IsVertical() ){
      NSelVert++;
    } else {
      NSelHori++;
    }

    FirstClu = false;
    PrevIsVert = Clu->IsVertical();
  }





  // We removed the clusters starting at the 3rd change or orientation.
  // So record only the first 2 changes.
  if (fNbOrientChange > 2){
      std::cout << " WARNING: The cluster orientation changed "<<fNbOrientChange<<" times !!! This is not good !" << std::endl; 
    fNbOrientChange = 2;
  }


  /*  
  if( (double)NSaturation > 0.3*(HVclu->size()) && fExcludeSaturatedClusters) {
    if( ND::tpcDebug().Seeding(DB_INFO)) {
      std::cout << " The number of saturated hits is too large. Use saturated waveforms for the seeding. " << std::endl; 
    }

    fExcludeSaturatedClusters = 0;
    PrepareClustersForSeeding(HVclu);
    fExcludeSaturatedClusters = 1;
  }
  else {
    if( ND::tpcDebug().Seeding(DB_VERBOSE)) {
  */
  


  std::cout << " ------- Hit selection prior to seeding."<<std::endl;
  std::cout << " Original number of clusters: " << HVclu->size() << std::endl; 
  std::cout << " Number of selected clusters:"<<std::endl;
  std::cout << "    Vertical:   " << NSelVert << std::endl;
  std::cout << "    Horizontal: " << NSelHori << std::endl;
  std::cout << "    Total:      " << (NSelVert + NSelHori) << std::endl;
  //std::cout << " Rejected by " << std::endl; 
  // std::cout << " restricted number of peaks       " << NMaxPeaks          << std::endl;
  //std::cout << " restricted number of saturations " << NSaturation        << std::endl;

}


//*****************************************************************************
// Calculate some quantities that are used by multiple seeding algorithms
// to speed things up.
void trex::TTPCSeeding::PrepareSeeding(std::vector<trex::TTRExCluster>& HVclu ){
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
  for (std::vector<trex::TTRExCluster>::iterator hit = HVclu->begin(); hit != HVclu->end(); hit++) {
    trex::TTRExHVCluster> clu = *hit;
    if (!clu){
      // TODO: proper exception
      throw;
    }
    if( !clu->IsOkForSeed() ) continue; 

    xclu = clu->X();
    zclu = clu->Z();
    yclu = clu->Y();

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

  for (::THitSelection::const_iterator hitit = HVclu->begin(); hitit != HVclu->end(); hitit++){
    ND::THandle<ND::TTPCHVCluster> clu = *hitit;
    if (!clu){
      // TODO: proper exception
      throw;
    }
    if( !clu->IsOkForSeed() ) continue;

    // Search for Zmid.
    if( TargetClu == NbSelectedClu )  {
      fZmid = clu->Z();
      fYmid = clu->Y();
      break;
    }
    NbSelectedClu++;
  }

    std::cout << " ======> fNbOrientChange "<< fNbOrientChange<< std::endl;
    std::cout << " First pt: " << fXfirst << ", " << fYfirst << ", " << fZfirst << std::endl;
    std::cout << " Last pt: " << fXlast << ", " << fYlast << ", " << fZlast << std::endl;
    std::cout << " Mid pt: " << fYmid << ", " << fZmid << std::endl;
}


//*****************************************************************************
// R2 (3 point seeding method)
//*****************************************************************************
double ND::TTPCSeeding::R2( ND::THandle<ND::THitSelection> HVclu, State &Helix){

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

  double phiFirst = TMath::ATan2((fYfirst-yl0),(fZfirst-zl0)); 
  // Only for printout 
  double phiLast =  TMath::ATan2((fYlast-yl0),(fZlast-zl0)); 

  // Assuming the sense of the track is correct we have:
  // Charge = -Helicity 
  double Helicity = CalculateRhoSign(yl0, zl0);

  double r2 = (fYfirst-yl0)*(fYfirst-yl0)+(fZfirst-zl0)*(fZfirst-zl0); 

  if( isnan(yl0) || isnan(zl0) ) r2 = 1.e-20;

  double R2rho = 1./TMath::Sqrt(r2);
  
  if( R2rho < RHOMIN ) R2rho = RHOMIN;

  R2rho *= Helicity;
  double R2theta = phiFirst - (Helicity * TMath::PiOver2());

  if( ND::tpcDebug().Seeding(DB_VERBOSE) ){
    std::cout << " R2:  zc  = " <<zl0 <<"      yc  = " <<yl0 << std::endl; 
    std::cout << " R2:  Helicity = "<<Helicity << std::endl; 
    std::cout << " R2:  phiFirst = "<<phiFirst <<"   phiLast = "<<phiLast << std::endl; 
    std::cout << " R2:  theta = "<<R2theta<<"      rho = " << R2rho << std::endl;
  }

  // At this point just use the projection of the clusters
  // onto the center of the cathode to figure out the delta phi angle (less than 1pi ?)

  double R2x = fXfirst;
  double R2y = fYfirst;
  double R2z = fZfirst;

  Helix = State(7);
  EVector outVect = EVector(7,0);
  EMatrix outCova = EMatrix(7,7,0);
  outVect[0] = R2x;
  outVect[1] = R2y;
  outVect[2] = R2z;
  outVect[3] = 0.0;
  outVect[4] = TMath::Sin(R2theta);
  outVect[5] = TMath::Cos(R2theta);

  TVector3 Pos(outVect[0], outVect[1], outVect[2]);
  TVector3 Dir(outVect[3], outVect[4], outVect[5]);
  double p;
  double q;
  TrackingUtils::Curvature_to_MomentumAndCharge(Pos,Dir,R2rho,p,q);
  outVect[6] = q/p;
  if( ND::tpcDebug().Seeding(DB_VERBOSE) ){
    std::cout << " R2:  dy       = " << outVect[4] << std::endl;
    std::cout << " R2:  dz       = " << outVect[5] << std::endl;
    std::cout << " R2:  momentum = " << p << std::endl;
  }

  Helix.set_hv(HyperVector(outVect,outCova));
  Helix.set_name(RP::representation,RP::pos_dir_curv);

  return FinalizeSeed(HVclu, R2rho, Helix);
}

//*****************************************************************************
// Riemann Fast Fit (Riemann sphere seeding method) 
//*****************************************************************************
double ND::TTPCSeeding::Riemann( ND::THandle<ND::THitSelection> HVclu, State &Helix){
  double scale = 1000.; 

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  double rx = 0.0; 
  double ry = 0.0;
  double rz = 0.0;

  for (ND::THitSelection::const_iterator hitit = HVclu->begin(); hitit != HVclu->end(); hitit++) {
    ND::THandle<ND::TTPCHVCluster> clu = *hitit;
    if( !clu->isOkForSeed() ) continue; 

    double yy = (clu->Y() - fYfirst)/scale;
    double zz = (clu->Z() - fZfirst)/scale;
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

  if( ND::tpcDebug().Seeding(DB_VERBOSE) ){
    std::cout << " Riemann:  zc  = " << u0 <<"      yc  = " << v0 << std::endl; 
    std::cout << " Riemann:  Helicity = "<<Helicity << std::endl; 
    std::cout << " Riemann:  phiFirst = "<<phiFirst <<"   phiLast = "<<phiLast << std::endl; 
    std::cout << " Riemann:  theta = "<<Rietheta<<"      rho = " << Rierho << std::endl;
  }

  // At this point just use the projection of the clusters
  // onto the center of the cathode to figure out the delta phi angle (less than 1pi ?)

  double Riex = fXfirst;
  double Riey = fYfirst;
  double Riez = fZfirst;

  Helix = State(7);
  EVector outVect = EVector(7,0);
  EMatrix outCova = EMatrix(7,7,0);
  outVect[0] = Riex;
  outVect[1] = Riey;
  outVect[2] = Riez;
  outVect[3] = 0.0;
  outVect[4] = TMath::Sin(Rietheta);
  outVect[5] = TMath::Cos(Rietheta);

  TVector3 Pos(outVect[0], outVect[1], outVect[2]);
  TVector3 Dir(outVect[3], outVect[4], outVect[5]);
  double p;
  double q;
  TrackingUtils::Curvature_to_MomentumAndCharge(Pos,Dir,Rierho,p,q);
  outVect[6] = q/p;
  if( ND::tpcDebug().Seeding(DB_VERBOSE) ){
    std::cout << " Riemann:  dy       = " << outVect[4] << std::endl;
    std::cout << " Riemann:  dz       = " << outVect[5] << std::endl;
    std::cout << " Riemann:  momentum = " << p << std::endl;
  }

  Helix.set_hv(HyperVector(outVect,outCova));
  Helix.set_hv(RP::sense, HyperVector(1,0));
  Helix.set_name(RP::representation,RP::pos_dir_curv);

  return FinalizeSeed(HVclu, Rierho, Helix);

}


//*****************************************************************************
double ND::TTPCSeeding::CalculateRhoSign(double Y0, double Z0) {
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
double ND::TTPCSeeding::FinalizeSeed( ND::THandle<ND::THitSelection> HVclu, double rho, State &finalState){
//*****************************************************************************
  double rms = 0.0;
  int ntot=0;
  double ypred,zpred;
  double yclu,zclu;

  EVector Vect = finalState.vector();
  EMatrix Cova = finalState.matrix();

  // First start by setting the sense
  ND::THitSelection::const_iterator Hit = HVclu->begin();
  ND::THandle<ND::TTPCHVCluster> Cluster = (*Hit);
  if (Cluster->IsVertical()){
    if ( Vect[5] > 0)
      finalState.set_hv(RP::sense, HyperVector(1,0));
    else
      finalState.set_hv(RP::sense, HyperVector(-1,0));
  } else {
    if ( Vect[4] > 0)
      finalState.set_hv(RP::sense, HyperVector(1,0));
    else
      finalState.set_hv(RP::sense, HyperVector(-1,0));
  }
    

  if( ND::tpcDebug().Seeding(DB_VVERBOSE) )
    std::cout<<" --> FinalizeSeed: Propagation"<<std::endl;
  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  State propagState = finalState;
  double deltaPhi = 0.0;
  State prevState;
  double ClusterX1;
  double ClusterX2 = 0.0;
  ClusterX1 = Cluster->CalibX();
  bool doDeltaPhi = true;
  double suspicious =0.0;
  for ( ; Hit != HVclu->end(); Hit++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = (*Hit);
    if( !Cluster->isOkForSeed() ) continue;
    if( ND::tpcDebug().Seeding(DB_VVERBOSE) )
      std::cout<<std::setw(20)<<" Cluster at: "<<Cluster->GetPosition().Y()<<"   "<<Cluster->GetPosition().Z()<<"   isVertical = "<<Cluster->IsVertical()<<std::endl;
    prevState = propagState;
    double length = 0.0;
    if (!TTPCRecPackUtils::PropagateToHVCluster(Cluster, propagState, length))
      continue;

    ypred = propagState.vector()[1];
    zpred = propagState.vector()[2];
    yclu = Cluster->Y();
    zclu = Cluster->Z();
    double delta = TMath::Sqrt((yclu-ypred)*(yclu-ypred)+(zclu-zpred)*(zclu-zpred));
    if( ND::tpcDebug().Seeding(DB_VVERBOSE) ){
      std::cout<<std::setw(20)<<" Propagation at: "<<ypred<<"   "<<zpred<<std::endl;
      std::cout<<std::setw(20)<<" => delta = "<<delta<<std::endl;
    }
    if ( fabs(delta) > 100.){
      if( ND::tpcDebug().Seeding(DB_VVERBOSE) )
        std::cout<<" --> Suspicious jump. Don't use it."<<std::endl;
      propagState = prevState;
      suspicious += 1.;
      continue;
    }
      
    if (doDeltaPhi)
      deltaPhi += length * rho;
    ClusterX2 = Cluster->CalibX();
    if (deltaPhi > TMath::Pi() && doDeltaPhi){
      doDeltaPhi = false;
    }
    rms += delta*delta; 
    ntot++;
  }

  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);

  if ( (suspicious/(double)ntot) > 0.7){
    if( ND::tpcDebug().Seeding(DB_INFO) )
      std::cout<<" ERROR: Too many suspicious jumps. Fail the seeding."<<std::endl;
    return 1.e99;
  }
  rms /= (double)ntot;

  double deltaX = ClusterX2-ClusterX1;
  Vect[3] = deltaX / TMath::Sqrt(deltaX*deltaX + (deltaPhi*deltaPhi)/(rho*rho));
  double Renorm = TMath::Sqrt((1 - Vect[3]*Vect[3])/(Vect[4]*Vect[4] + Vect[5]*Vect[5]));
  Vect[4] *= Renorm;
  Vect[5] *= Renorm;
  Vect[6] *= sqrt(1. - Vect[3]*Vect[3]);
  finalState.set_hv(HyperVector(Vect,Cova));

  if( ND::tpcDebug().Seeding(DB_VERBOSE) ){
    std::cout << " X angle calculation:  deltaX  = " << deltaX << std::endl; 
    std::cout << " X angle calculation:  deltaPhi = "<<deltaPhi << std::endl; 
    std::cout << " X angle calculation:  Renorm = "<<Renorm << std::endl; 
    std::cout << " X angle calculation:  rms = "<<rms << std::endl; 
    std::cout << " cosX = "<<Vect[3] << std::endl; 
  }

  // Something went wrong, don't trust the result
  if ( ! (rms > 0.0) )
    rms = 1.e99;

  return rms;
}


//*****************************************************************************
bool ND::TTPCSeeding::IsResultValid( ND::THandle<ND::THitSelection> HVclu, State &Result){
//*****************************************************************************
  EVector outVect = Result.vector();
  // For bad fits give tracks in X only and we can't fit those anyway.
  if (outVect[4] == 0.0 && outVect[5] == 0.0){
    if( ND::tpcDebug().Seeding(DB_VERBOSE) )
      std::cout<<" Seed invalid: The y and z directions are zero !"<<std::endl;
    return false;
  }

  // If the first cluster is vertical (horizontal), the z (y) direction cannot be zero.
  // TODO: maybe this constraint could be more stringent.
  ND::THandle<ND::TTPCHVCluster> HV = *(HVclu->begin());
  if ( HV->IsVertical() && outVect[5] == 0.0){
    if( ND::tpcDebug().Seeding(DB_VERBOSE) )
      std::cout<<" Seed invalid: The first cluster is vertical and the z direction is zero !"<<std::endl;
    return false;
  }
  if ( ( !HV->IsVertical()) && outVect[4] == 0.0){
    if( ND::tpcDebug().Seeding(DB_VERBOSE) )
      std::cout<<" Seed invalid: The first cluster is horizontal and the y direction is zero !"<<std::endl;
    return false;
  }

  return true;
}
