#include "TTPCLikFitPath.hxx"

#include <TOADatabase.hxx>
#include <HEPUnits.hxx>
#include <TOARuntimeParameters.hxx>
#include <TGeomInfo.hxx>
#include <TFieldManager.hxx>
#include <TrackingUtils.hxx>

#include "TTPCHVCluster.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCHelixPropagator.hxx"
#include "TTPCDebug.hxx"
#include "TTPCRecPackUtils.hxx"
// #include "TFieldManager.hxx"


////////////////////// Minuit global section //////////////////////////////
// We need these global variables for Minuit's fcn function.

TTPCLikFitPath *LikFitPtr;

bool gLastCall = false;



// *********************************************************************************
static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  if (iflag == 3) { // will come here only after the fit is finished.    
    gLastCall = true;
  }

  f = LikFitPtr->log_likelihood(par); 

  gLastCall = false;
}

//////////////////////////////////////////////////////////////////////////


// *********************************************************************************
// TTPCLikFitPath constructor
TTPCLikFitPath::TTPCLikFitPath(bool CalculatorMode)  {
  fCalculatorMode = CalculatorMode;
  if( ! fCalculatorMode){
    LikFitPtr = this;

    fMinuit = NULL;
  }

  // If we have no TPCs in this geometry, give up constructor before we break something.
  if(ND::TGeomInfo::Get().TPC().NumberofTPC() == 0){
    ND280Log("No TPC modules in geometry; stopping TTPCLikFitPath constructor.");
    return;
  }

  fMinuitMaxIterations = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.MinuitMaxIt");
  fMinuitPrintLevel = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.MinuitPrintLevel");

  fPadWidth = ND::TGeomInfo::Get().TPC().GetPadXPitch();
  fPadHeight = ND::TGeomInfo::Get().TPC().GetPadYPitch();

  fFitSigma    = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.FitSigma");
  fFitSigmaSeparately = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.FitSigmaSeparately");
  fFitRho      = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.FitRho");
  fFitX        = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.FitX");
  fFitDir      = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.FitDir");

  fMinimumCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MinimumCharge");
  fMaximumCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MaximumCharge");
  fMaxDeltaDrift = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MaxDeltaDrift");
  fMinPredIntCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MinPredIntCharge");
  fNoise = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.Noise");

  fErrorWeightY = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.WeightY");
  fErrorWeightX = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.WeightX");

  double transDiff = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.TransDiff");

  // Multiple scattering parameters
  fEnabledMScCorrection = (bool) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.EnabledMScCorrection");
  fErrorWeightYMSc = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.WeightYMSc");
  fErrorWeightXMSc = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.WeightXMSc");
  fRadLenHe = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.RadLenHe");
  fMuonMass = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MuonMass");
  fMomChangeForNewMSc = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.MomChangeForNewMSc");

  fInitialSigma = transDiff*transDiff;

  fFirstXYZfit = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.TryXYZFitFirst");
  fForceXZfit = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.ForceXZFit");
  fFinalXYZfit = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.FinalXYZFit");

  
  fExcludeClusterWithManyPeaks = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.ExcludeClusterWithManyPeaks");
  fExcludeSaturatedClusters = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.ExcludeSaturatedClusters");
  fExcludeClusterAtEdge = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.ExcludeClusterAtMMEdge");
  fExcludeSuspiciousPadTiming = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.ExcludeSuspicious");
  fMaxPadsPerCluster = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.MaxPadsPerCluster");
  fMinNumberOfClusters  = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.MinNbClusters");
  fMinNumberOfClustersForCalc  = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.MinNbClustersForCalculator");

  fApplyResidualDependentPenalty = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.ApplyResidualDependentPenalty");
  fPenaltyThreshold = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.PenaltyThreshold");
  fPenaltyCoefficient = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.LikFit.PenaltyCoefficient");

  // Short cut for the multiple scattering correction
  fBFieldAt0 = ND::TFieldManager::GetFieldValue(TVector3(0.0,0.0,0.0)).Mag() / unit::tesla;

  fStoreLklhdWithMScCorr = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.StoreLklhdWithMScCorr");

  Reset();
}


// *********************************************************************************
// TTPCLikFitPath destructor
TTPCLikFitPath::~TTPCLikFitPath() {
  if (fMinuit && !fCalculatorMode) delete fMinuit;
}




// *********************************************************************************
// Reset all the class variables prior to fitting
void TTPCLikFitPath::Reset(void){
  fFitYPosParam = true; // just a default.
  fMScCorr = 1.;

  fFitClu = ND::THandle<ND::THitSelection> ();

  fMeanDrift=0.0;

  fFitSteps = 0;
  fReliableFit = true;
}



// *********************************************************************************
// DefaultFixedParameters
void TTPCLikFitPath::DefaultFixedParameters(void){

  if( !fMinuit ) return; 

  if( !fFitX ){
    fMinuit->FixParameter(YPARAM);
    fMinuit->FixParameter(XPARAM);
  }

  if( !fFitDir ){
    fMinuit->FixParameter(TANXPARAM);
    fMinuit->FixParameter(TANYORZPARAM);
  }

  if( !fFitRho ){
    fMinuit->FixParameter(CURVPARAM);
  }

  if( fFitSigmaSeparately || !fFitSigma ){
    fMinuit->FixParameter(SGMPARAM);
  }

  fMinuit->FixParameter(ZPARAM);

}



// *********************************************************************************
TVector3 TanToDir(double TanX, double TanYorZ, bool WithRespToZ){
  TVector3 Dir;
  double theta = TMath::ATan(TanYorZ);
  if (WithRespToZ){
    Dir.SetY(TMath::Sin(theta));
    Dir.SetZ(TMath::Cos(theta));
    Dir.SetX(TanX * Dir.Z());
  } else {
    Dir.SetY( TMath::Cos(theta));
    Dir.SetZ( TMath::Sin(theta));
    Dir.SetX(TanX * Dir.Y());
  }
  Dir = Dir.Unit();
  return Dir;
}


// *********************************************************************************
// Setup variables and stuff for the minimization
void TTPCLikFitPath::SetupLogLklhdMinimizer(const ND::THandle<ND::TTPCPath> Path, bool UseSeedAsInit){
  if (fMinuit) delete fMinuit;

  fMinuit = new TMinuit(NPARAM);
  fMinuit->SetFCN(fcn);
  fMinuit->SetPrintLevel(fMinuitPrintLevel);

  // Fit the Z origin of the track rather than Y for high angle tracks.
  ND::THandle<ND::TTPCHVCluster> firstClu = *(Path->GetHits()->begin());
  if (firstClu->IsVertical())
    fFitYPosParam = true;
  else
    fFitYPosParam = false;

//  for (ND::THitSelection::const_iterator Hit = Path->GetHits()->begin(); Hit != Path->GetHits()->end(); Hit++) {
//    ND::THandle<ND::TTPCHVCluster> clu = *Hit;
//  std::cout<<"    "<<clu->GetPosition().X()<<"    "<<clu->GetPosition().Y()<<"    "<<clu->GetPosition().Z()<<std::endl;
//    
//  }

  ND::helixPropagator().Reset();
  State initState = Path->GetFrontSeedState();

  if (UseSeedAsInit){
    if (!Path->HasSeedState()){
      std::cout<<"TTPCLikFitPath SetupLogLklhdMinimizer: ERROR: There is no seed state available !"<<std::endl;
      //TODO: proper exception
      throw;
    }
    initState = Path->GetFrontSeedState();
  } else {
    if (!Path->HasFitState()){
      std::cout<<"TTPCLikFitPath SetupLogLklhdMinimizer: ERROR: There is no fit state available !"<<std::endl;
      //TODO: proper exception
      throw;
    }
    initState = Path->GetFrontFitState();
  }

  if ( ND::tpcDebug().LikFit(DB_VERBOSE))
    std::cout<<"TTPCLikFitPath::SetupLogLklhdMinimizer: Minimization for path Id "<<Path->GetId()<<std::endl;

  fPathLength = Path->GetLength();
  StateToInitValues(initState);

  fInitValues[SGMPARAM] = fInitialSigma;

  fStep[XPARAM]       = 0.001;
  fStep[YPARAM]       = 0.001;
  fStep[ZPARAM]       = 0.001;
  fStep[TANXPARAM]    = 0.001;
  fStep[TANYORZPARAM] = 0.001;
  fStep[CURVPARAM]    = 0.001;
  fStep[SGMPARAM]     = 0.001;

  ResetMinuitParam();

  // Just to clear things up
  fFitResults.FitState = State();
  fFitResults.Sigma = fInitialSigma;
  fFitResults.eSigma = 0.0;
  fFitResults.fitSteps = 0;

  // Probably not needed
  fLogLklhd.Total = 0.0;
  fLogLklhd.X = 0.0;
  fLogLklhd.HV = 0.0;

}


// *********************************************************************************
// Setup variables and stuff for the minimization
bool TTPCLikFitPath::SetupLogLklhdCalculator(State helixState, ND::THandle<ND::THitSelection> inputClusters, double inputLength){

  ND::helixPropagator().Reset();

  if (inputClusters->size() < fMinNumberOfClustersForCalc){
    return false;
  }

  ND::THandle<ND::TTPCHVCluster> firstClu = *(inputClusters->begin());

  // We don't fit the coordinates here but we still need to know the orientation
  // for the proper initialization of the HelixPropagator
  if (firstClu->IsVertical())
    fFitYPosParam = true;
  else
    fFitYPosParam = false;

  
  fPathLength = inputLength;
  StateToInitValues(helixState);

  fFitClu = inputClusters;

  // Probably not needed
  fLogLklhd.Total = 0.0;
  fLogLklhd.X = 0.0;
  fLogLklhd.HV = 0.0;

  return true;
}


// *********************************************************************************
// Load initial values of the likelihood based on given state
void TTPCLikFitPath::StateToInitValues( State inputState){

  if(inputState.name(RP::representation) != RP::pos_dir_curv)
    RP::rep().convert(inputState, RP::pos_dir_curv);
  // Quick and easy way to get the quadrant and the curvature.
  // TODO: should check that this Loading worked ... what to do if it didn't ?
  ND::helixPropagator().InitHelixPosDirQoP(inputState, fFitYPosParam);
  double Momentum = fabs(1./inputState.vector()[6]);
  
  double *tmpArray = new double[6];
  ND::helixPropagator().GetHelixPosTanCurv(tmpArray);
  fInitValues[XPARAM] = tmpArray[0];
  fInitValues[YPARAM] = tmpArray[1];
  fInitValues[ZPARAM] = tmpArray[2];
  fInitValues[TANXPARAM] = tmpArray[3];
  fInitValues[TANYORZPARAM] = tmpArray[4];
  fInitValues[CURVPARAM] = tmpArray[5];

  fInitQoP = inputState.vector()[6];
  
//// Old debug code that could be useful
//   if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
// //  if ( 1){
//     if (fFitYPosParam){
//       RP::rep().convert(inputState, RP::slopes_curv_z);
//     } else {
//       RP::rep().convert(inputState, RP::slopes_curv_y);
//     }
//     EVector seedVect = inputState.vector();
//     std::cout<<" ----> Compare seed tangent from RP and my code "<<std::endl;
//     std::cout<<"  "<< tmpArray[3] <<"  "<< seedVect[3] <<std::endl;
//     std::cout<<"  "<< tmpArray[4] <<"  "<< seedVect[4] <<std::endl;
// 
//     for( int i = 0; i < 5; i++)
//       tmpArray[i] = seedVect[i];
//     ND::helixPropagator().ReloadHelixPosTanCurv(tmpArray);
// 
//     double *tmpArray2 = new double[7];
//     State tmpState = Path->GetFrontSeedState();
//     ND::helixPropagator().GetHelixPosDirCurv(tmpArray2);
//     std::cout<<" ----> Compare seed position, direction, and curvature from RP and my code "<<std::endl;
//     for( int i = 0; i < 7; i++)
//       std::cout<<"  "<< tmpArray2[i] <<"  "<< tmpState.vector()[i] <<std::endl;
// 
//   }

  delete[] tmpArray;

  // Compute the multiple scattering correction
  if( fEnabledMScCorrection ){
    ComputeMScCorrection(Momentum);
  }

}


// *********************************************************************************
void TTPCLikFitPath::ComputeMScCorrection(double momentum){
  double beta = momentum/std::sqrt((momentum*momentum) + (fMuonMass*fMuonMass));
  double xx0 = fPathLength / fRadLenHe;
  fMScCorr = ( (13.6/(momentum*beta))*std::sqrt(xx0)*(1+0.038*std::log(xx0)) );
  if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
    std::cout<<" => MScCorr calculation: "<<std::endl;
    std::cout<<"    - momentum   = "<<momentum<<std::endl;
    std::cout<<"    - length     = "<<fPathLength<<std::endl;
    std::cout<<"    - correction = "<<fMScCorr<<std::endl;
  }
}


// *********************************************************************************
void TTPCLikFitPath::SelectClusters(ND::THandle<ND::THitSelection> inputClu, double XDirection, ClusterSelection &CluSel){

  // Select clusters
  CluSel.NMaxPeaks = 0; 
  CluSel.NSelVert = 0; 
  CluSel.NSelHori = 0; 
  CluSel.NHoriMMEdge = 0; 
  CluSel.NVertMMEdge = 0; 
  CluSel.NOutOfChargeWindow = 0; 
  CluSel.NOutDeltaDrift = 0; 
  CluSel.NTooManyPadsPerClu = 0; 
  CluSel.NSuspiciousPadTiming = 0; 
  CluSel.NSaturation = 0; 
  
  for (ND::THitSelection::const_iterator tmpClu = inputClu->begin(); tmpClu != inputClu->end(); tmpClu++) {
    ND::THandle<ND::TTPCHVCluster> Clu = *tmpClu;

    Clu->SetOkForFit(true);  // Make sure that we start with fresh sample.

    if( Clu->GetMaxNPeaks() > 1 && fExcludeClusterWithManyPeaks )  {
      Clu->SetOkForFit(false);
      CluSel.NMaxPeaks++; 
      continue;
    }

    if( fExcludeClusterAtEdge ) {
      if( Clu->IsAtVertEdge() > 0  && Clu->IsHorizontal()){
        Clu->SetOkForFit(false);
        CluSel.NVertMMEdge++; 
        continue;
      } else if( Clu->IsAtHoriEdge() > 0  && Clu->IsVertical()){
        Clu->SetOkForFit(false);
        CluSel.NHoriMMEdge++; 
        continue;
      } 
    }

    if( Clu->GetCharge()*fabs(XDirection) < fMinimumCharge ||  Clu->GetCharge()*fabs(XDirection) > fMaximumCharge ) {
      Clu->SetOkForFit(false);
      CluSel.NOutOfChargeWindow++;
      continue;
    }
  
    if( Clu->GetDeltaDrift()*XDirection >= fMaxDeltaDrift ) {
      Clu->SetOkForFit(false);
      CluSel.NOutDeltaDrift++;
      continue; 
    }

    if( Clu->GetHits().size() > fMaxPadsPerCluster)  {
      Clu->SetOkForFit(false);
      CluSel.NTooManyPadsPerClu++;
      continue;
    } 
    if( Clu->HasSuspiciousPadTiming() && fExcludeSuspiciousPadTiming ){
      Clu->SetOkForFit(false);
      CluSel.NSuspiciousPadTiming++;
      continue;
    }
    if( Clu->GetNSaturated() > 0  && fExcludeSaturatedClusters ) {
      Clu->SetOkForFit(false);
      CluSel.NSaturation++;
      continue;
    }

    // Cluster selected !
    if( Clu->IsVertical() ){
      CluSel.NSelVert++;
    } else {
      CluSel.NSelHori++;
    }
  }

}


// *********************************************************************************
void TTPCLikFitPath::PrepareClustersForFitting(ND::THandle<ND::THitSelection> inputClu, ND::THandle<ND::THitSelection> outputClu, double XDirection){
  ClusterSelection ClusterSelResults;

  // ==> PASS 1: with all the settings as default
  SelectClusters(inputClu, XDirection, ClusterSelResults);


  // ==> PASS 2: with saturation rejection off
  if( (double)ClusterSelResults.NSaturation > 0.3*(inputClu->size()) && fExcludeSaturatedClusters) {
    if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
      std::cout << " The number of saturated hits is too large. Use saturated waveforms for the likelihood fit. " << std::endl; 
    }

    fExcludeSaturatedClusters = 0;
    SelectClusters(inputClu, XDirection, ClusterSelResults);
    // We are going to save the clusters at the end of this pass so clear here.
    fExcludeSaturatedClusters = 1;
  }

  // ==> PASS 3: with rejection of clusters at edge of MM off
  if( (ClusterSelResults.NSelVert+ClusterSelResults.NSelHori) < fMinNumberOfClusters && fExcludeClusterAtEdge){
    if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
      std::cout << " The number of selected clusters is "<<(ClusterSelResults.NSelVert+ClusterSelResults.NSelHori)<<" below the minimum of "<<fMinNumberOfClusters<<". Use waveforms at the edge of the MMs for an unreliable likelihood fit." << std::endl; 
    }
    fExcludeClusterAtEdge = 0;
    SelectClusters(inputClu, XDirection, ClusterSelResults);
    fReliableFit = false;
    // We are going to save the clusters at the end of this pass so clear here.
    fExcludeClusterAtEdge = 1;
  }

  if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
    std::cout << " ------- Hit selection prior to fitting"<<std::endl;
    std::cout << " Original number of clusters: " << inputClu->size() << std::endl; 
    std::cout << " Number of selected clusters:"<<std::endl;
    std::cout << "    Vertical:   " << ClusterSelResults.NSelVert << std::endl;
    std::cout << "    Horizontal: " << ClusterSelResults.NSelHori << std::endl;
    std::cout << "    Total:      " << (ClusterSelResults.NSelVert + ClusterSelResults.NSelHori) << std::endl;
    std::cout << " Rejected by " << std::endl; 
    std::cout << " restricted number of peaks       " << ClusterSelResults.NMaxPeaks             << std::endl;
    std::cout << " suspicious pad timing            " << ClusterSelResults.NSuspiciousPadTiming  << std::endl;
    std::cout << " restricted number of saturations " << ClusterSelResults.NSaturation           << std::endl;
    std::cout << " charge window                    " << ClusterSelResults.NOutOfChargeWindow    << std::endl;
    std::cout << " max number of pads per row       " << ClusterSelResults.NTooManyPadsPerClu    << std::endl;
    std::cout << " Delta drift                      " << ClusterSelResults.NOutDeltaDrift        << std::endl; 
    std::cout << " MM horizontal edge               " << ClusterSelResults.NHoriMMEdge           << std::endl;
    std::cout << " MM vertical edge                 " << ClusterSelResults.NVertMMEdge           << std::endl;
  }
  
  for (ND::THitSelection::const_iterator tmpClu = inputClu->begin(); tmpClu != inputClu->end(); tmpClu++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
    if( !Cluster->isOkForFit() ) { continue;}  // Check that the plane is actually enabled.
    outputClu->push_back(Cluster);
    
  }

}



// *********************************************************************************
// Simply reset the minuit parameters to their starting values
void TTPCLikFitPath::ResetMinuitParam(){
  fMinuit->mnparm(XPARAM,      "X0",      fInitValues[XPARAM],      fStep[XPARAM],      0.,0.,fIErrorFlag);
  fMinuit->mnparm(YPARAM,      "Y0",      fInitValues[YPARAM],      fStep[YPARAM],      0.,0.,fIErrorFlag);
  fMinuit->mnparm(ZPARAM,      "Z0",      fInitValues[ZPARAM],      fStep[ZPARAM],      0.,0.,fIErrorFlag);
  fMinuit->mnparm(TANXPARAM,   "TANXO",   fInitValues[TANXPARAM],   fStep[TANXPARAM],   0.,0.,fIErrorFlag);
  fMinuit->mnparm(TANYORZPARAM,"TANYORZO",fInitValues[TANYORZPARAM],fStep[TANYORZPARAM],0.,0.,fIErrorFlag);
  fMinuit->mnparm(CURVPARAM,   "CURV0",   fInitValues[CURVPARAM],   fStep[CURVPARAM],   0.,0.,fIErrorFlag);
  fMinuit->mnparm(SGMPARAM,    "SIGMA",   fInitValues[SGMPARAM],    fStep[SGMPARAM],    0.,0.,fIErrorFlag);

  //////////////////////// RECPACK CODE. This will eventually be deleted. ////////////////////////
  // Reset the State used during the propagation of the helix prediction.
  EVector predVect = EVector(6,0);
  EMatrix predCova = EMatrix(6,6,0);
  predVect[0] = fInitValues[XPARAM];
  predVect[1] = fInitValues[YPARAM];
  predVect[2] = fInitValues[ZPARAM];
  predVect[3] = fInitValues[TANXPARAM];
  predVect[4] = fInitValues[TANYORZPARAM];
  TVector3 Position(predVect[0], predVect[1], predVect[2]);
  TVector3 Direction = TanToDir(predVect[3], predVect[4], fFitYPosParam);
  double p, q;
  // RecPack state store Q over P so we need to convert back to the curvature.
  if (TrackingUtils::Curvature_to_MomentumAndCharge(Position, Direction, x[CURVPARAM], p, q)){
    predVect[5] = q/p;
  } else {
    predVect[5] = 0.0;
  }
}

// *********************************************************************************
// Minimize the log likelihood for the given hits
int TTPCLikFitPath::LogLklhdMinimizer(ND::THandle<ND::THitSelection> inputClusters){
  GetReadyForMinimization(inputClusters);

  fFitXProj = true; 
  fFitYProj = true; 

  if( ! fFitYPosParam ){
    fMinuit->Release(ZPARAM);
    fMinuit->FixParameter(YPARAM);
  }

  
  if( ND::tpcDebug().LikFittedClusters(DB_INFO)){
    std::cout << " ================= LogLklhdMinimizer ================" << std::endl; 
    TTPCUtils::HVClustersPrintout(fFitClu, ND::tpcDebug().LikFittedClusters(DB_VERBOSE));
    std::cout << " ----------------------------------------------------" << std::endl; 
  }


  // Print results
  double arglist[4];
  // Now ready for minimization step

  arglist[0] = fMinuitMaxIterations;
  fIErrorFlag = 9999; 

  if( fFirstXYZfit )  {
    arglist[0] = fMinuitMaxIterations;
    fIErrorFlag = 9999; 
    fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);
    if( ND::tpcDebug().LikFit(DB_VERBOSE))
      std::cout << " XYZ MIGRAD ERROR " << fIErrorFlag << std::endl; 
  }
  fFitSteps += 10;
    

  if( fIErrorFlag ) {
    // ---> The first XYZ fit failed: Disconnect the X and YZ fits temporarily
    if ( ND::tpcDebug().LikFit(DB_VERBOSE) && fFirstXYZfit){
      std::cout << " >>>> XYZ fit has ierflg =" << fIErrorFlag << ", Try disconnecting X and Y/Z " << std::endl;
    }

    ResetMinuitParam();

    // ---> Fit the YZ projection only
    fFitXProj = false;
    fFitYProj = true; 

    fMinuit->FixParameter(XPARAM);
    fMinuit->FixParameter(TANXPARAM);

    arglist[0] = fMinuitMaxIterations;
    fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);
    fFitSteps += 10;
    
    bool FixMomentum = false;
    if( fIErrorFlag ) {
      if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
        std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Failing Y coordinate fit ( " << fIErrorFlag << " ) => Try with fixed momentum"<<std::endl;
        // std::cout << x0 << "  " << y0 << "  " << tx0 << "  " << phi0 << "  " << rho0 << std::endl; 
      }

      ResetMinuitParam();
      // fMinuit->mnparm(CURVPARAM,   "CURV0",   5.e-10,   fStep[CURVPARAM],   0.,0.,fIErrorFlag);
      fMinuit->FixParameter(CURVPARAM);
      fReliableFit = false;
      FixMomentum = true;
      arglist[0] = fMinuitMaxIterations;
      fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);
      if( fIErrorFlag ) {
        if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
          std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Definitely failing Y coordinate fit ( " << fIErrorFlag << " )"<<std::endl;
        }
        return fIErrorFlag; 
      }
    }

    // ---> Fit the X parameters only
    fFitXProj = true; 
    fFitYProj = false; 

    fMinuit->Release(XPARAM);
    fMinuit->Release(TANXPARAM);
    fMinuit->FixParameter(YPARAM);
    fMinuit->FixParameter(ZPARAM);
    fMinuit->FixParameter(TANYORZPARAM);
    fMinuit->FixParameter(CURVPARAM);
    fMinuit->FixParameter(SGMPARAM);

    // The TX is a delicate parameter for the seeding so, we do a scan before trying the fit. 
    int istat; 
    char namecommand[32]; 
    sprintf(namecommand,"scan %d",TANXPARAM+1); 
    fMinuit->mncomd(namecommand,istat);

    arglist[0] = fMinuitMaxIterations;
    fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);   
    fFitSteps += 10;

    if( fIErrorFlag ) {
      if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
        std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Failing X coordinate fit " << std::endl; 
      }
      return fIErrorFlag; 
    }

    // ---> Try again to fit XYZ all at once at least to get the proper error calculation.
    fFitXProj = true; 
    fFitYProj = true; 

    if (fFitYPosParam)
      fMinuit->Release(YPARAM);
    else
      fMinuit->Release(ZPARAM);

    fMinuit->Release(TANYORZPARAM);
    fMinuit->Release(CURVPARAM);
    // We couldn't fit YZ projection with variable momentum, no point trying to do it here.
    if( FixMomentum)
      fMinuit->FixParameter(CURVPARAM);
    if( fFitSigma && !fFitSigmaSeparately)
      fMinuit->Release(SGMPARAM);

    if ( ND::tpcDebug().LikFit(DB_VERBOSE))
      std::cout << " Independent YZ & XZ MIGRAD was successful " << fIErrorFlag << std::endl;
    if( !fIErrorFlag  ) {
      if( ND::tpcDebug().LikFit(DB_VERBOSE) ) std::cout << " Refit with both coordinates " << std::endl; 
      arglist[0] = fMinuitMaxIterations;
      fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);
      fFitSteps += 10;
    }
  }


  // If the Seed momentum was way off, the multiple scattering correction
  // is also way off ...
  // Let's recalculate the correction and the reminimize to get the right
  // error but only when the momentum has changed a lot to avoid rerunning
  // the minimization too often.
  if( fEnabledMScCorrection && !fIErrorFlag ){
    double *tmpArray = new double[7];
    ND::helixPropagator().GetHelixPosDirQoP(tmpArray);
    double MomChange = fabs((fabs(tmpArray[6]) - fabs(fInitQoP)) / fInitQoP);
    if( MomChange > fMomChangeForNewMSc ) {
      if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
        MinuitPrintout();
        std::cout<<" >>>> The momentum has changed by "<<(MomChange)<<" between the initial value and the likelihood fit."<<std::endl;
        std::cout<<"      Reminimize with the fit momentum for the MSc correction !"<<std::endl;
      }

      double tmpMomentum = 1./fabs(tmpArray[6]);
      ComputeMScCorrection(tmpMomentum);
      arglist[0] = fMinuitMaxIterations;
      fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);
    }
    delete[] tmpArray;
  }


  if ( ND::tpcDebug().LikFit(DB_VERBOSE)){
    MinuitPrintout();
  }


  // Store the results from the fit, in particular the covariance matrix
  // before overwriting it with the independent fit of the sigma parameter.
  StoreFittedState();

  if( fFitSigmaSeparately && !fIErrorFlag  ) {
    int IErrorForSigma = 0;
    // Fit Sigma by itself since we don't need the correlations
    // and it can make the fit go crazy when fitted with the other parameters.
    fMinuit->FixParameter(XPARAM);
    fMinuit->FixParameter(TANXPARAM);
    fMinuit->FixParameter(YPARAM);
    fMinuit->FixParameter(ZPARAM);
    fMinuit->FixParameter(TANYORZPARAM);
    fMinuit->FixParameter(CURVPARAM);
    fMinuit->Release(SGMPARAM);
    fMinuit->mnparm(SGMPARAM, "SIGMA", fInitValues[SGMPARAM], fStep[SGMPARAM], 0, 0., IErrorForSigma);
    fMinuit->mnexcm("MINI",arglist,1,IErrorForSigma);
    if ( fIErrorFlag ){
      fFitResults.Sigma = fInitialSigma;
      fFitResults.eSigma = 999.999;
    } else {
      StoreFittedSigma();
    }
  } else {
    StoreFittedSigma();
  }

  CleanUpAfterMinimization();

  return fIErrorFlag;
}




// *********************************************************************************
// Does little preparations just before minimizing
void TTPCLikFitPath::GetReadyForMinimization(ND::THandle<ND::THitSelection> inputClusters){
  fFitClu = inputClusters;

  // Find the mean drift distance for this track
  // used in the fit for sigma0, sigma1 in new diffusion fit method.
  fMeanDrift=0.0;
  int nmeas=0;
  for (ND::THitSelection::const_iterator tmpClu = fFitClu->begin() ; tmpClu != fFitClu->end(); tmpClu++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
    fMeanDrift += Cluster->GetDriftDistance();
    nmeas++;
    
  }
  if (nmeas>0) fMeanDrift /= double( nmeas );


  fMinuit->mncler();

  ResetMinuitParam();

  DefaultFixedParameters();

}


// *********************************************************************************
// Minimize the log likelihood for the given hits
int TTPCLikFitPath::SimpleLogLklhdMinimizer(bool ParamIsFree[NPARAM]){

  for ( unsigned int pc = 0; pc < NPARAM; pc++){
    if (ParamIsFree[pc])
      fMinuit->Release(pc);
    else 
      fMinuit->FixParameter(pc);
  }

  if( ! fFitYPosParam ){
    fMinuit->Release(ZPARAM);
    fMinuit->FixParameter(YPARAM);
  }

  // Print results
  double arglist[4];
  // Now ready for minimization step

  arglist[0] = fMinuitMaxIterations;
  fIErrorFlag = 9999; 
  fMinuit->mnexcm("MINI",arglist,1,fIErrorFlag);
  

  // Store the results from the fit, in particular the covariance matrix
  StoreFittedState();

  return fIErrorFlag;
}


// *********************************************************************************
// Little clean up after minimizing
void TTPCLikFitPath::CleanUpAfterMinimization(){
  fFitClu = ND::THandle<ND::THitSelection> ();
}


// *********************************************************************************
void TTPCLikFitPath::MinuitPrintout(){
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    fMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    fMinuit->mnprin(3,amin);
}

// *********************************************************************************
// *********************************************************************************
// Minimize the log likelihood for the given hits
TTPCLogLikelihood TTPCLikFitPath::LogLklhdCalculator(){

  if( ND::tpcDebug().LikFittedClusters(DB_INFO)){
    std::cout << " ================= LogLklhdCalculator ===============" << std::endl; 
    TTPCUtils::HVClustersPrintout(fFitClu, ND::tpcDebug().LikFittedClusters(DB_VERBOSE));
    std::cout << " ----------------------------------------------------" << std::endl; 
  }


  fFitXProj = true; 
  fFitYProj = true; 

  // Find the mean drift distance for this set of clusters
  // used in the fit for sigma0, sigma1 in new diffusion fit method.
  fMeanDrift=0.0;
  int nmeas=0;
  for (ND::THitSelection::const_iterator tmpClu = fFitClu->begin() ; tmpClu != fFitClu->end(); tmpClu++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
    fMeanDrift += Cluster->GetDriftDistance();
    nmeas++;
    
  }
  if (nmeas>0) fMeanDrift /= double( nmeas );

  // TODO: decide here to do both projections or may only X or only Y

  double *tmpArray = new double[7];
  tmpArray[0] = fInitValues[XPARAM];
  tmpArray[1] = fInitValues[YPARAM];
  tmpArray[2] = fInitValues[ZPARAM];
  tmpArray[3] = fInitValues[TANXPARAM];
  tmpArray[4] = fInitValues[TANYORZPARAM];
  tmpArray[5] = fInitValues[CURVPARAM];
  tmpArray[6] = fInitialSigma;

  log_likelihood(tmpArray);
  
  delete[] tmpArray;
  fFitClu = ND::THandle<ND::THitSelection> ();

  return fLogLklhd;
}


// *********************************************************************************
// The logLikelihood
// function to calculate the likelihood for all rows -> to get sigma_w
// *********************************************************************************
double TTPCLikFitPath::log_likelihood(double* x){
  double llX  = 0.0;
  double llHV = 0.0;
  double llX_NoMSc  = 0.0;
  double llHV_NoMSc = 0.0;

  if( ND::tpcDebug().LikFit(DB_VVERBOSE)) {
    std::cout<<" ________________ MINUIT TRIAL ______________________"<<std::endl;
    std::cout<<" x         = "<<x[0]<<std::endl;
    std::cout<<" y         = "<<x[1]<<std::endl;
    std::cout<<" z         = "<<x[2]<<std::endl;
    std::cout<<" tanX      = "<<x[3]<<std::endl;
    std::cout<<" tanYorZ   = "<<x[4]<<std::endl;
    std::cout<<" curvature = "<<x[5]<<std::endl;
    std::cout<<" sigma     = "<<x[6]<<std::endl;
  }
  fSigma = fabs(x[SGMPARAM]);

  if( fFitYProj ) {
    ND::helixPropagator().ReloadHelixPosTanCurv(x);
    double MScWeightYZ = fMScCorr * fMScCorr * fErrorWeightYMSc * fErrorWeightYMSc;
    double loglikehood = log_likelihoodHV();
    double Weight = 1./(fErrorWeightY * fErrorWeightY);
    llHV_NoMSc = loglikehood * Weight;
    Weight += 1./MScWeightYZ;
    llHV = loglikehood * Weight;
  }

  if( fFitXProj ) {
    ND::helixPropagator().ReloadHelixPosTanCurv(x);
    double MScWeightX = fMScCorr * fMScCorr * fErrorWeightXMSc * fErrorWeightXMSc;
    double loglikehood = log_likelihoodX();
    double Weight = 1./(fErrorWeightX * fErrorWeightX);
    llX_NoMSc = loglikehood * Weight;
    Weight += 1./MScWeightX;
    llX = loglikehood * Weight;
  }
//  std::cout<<" llHV = "<<llHV<<"    llX = "<<llX<<std::endl;

  if ( fStoreLklhdWithMScCorr){
    fLogLklhd.Total = llHV+llX;
    fLogLklhd.X = llX;
    fLogLklhd.HV = llHV;
    if( ND::tpcDebug().LikFit(DB_VVERBOSE) ) {
      std::cout<<" ________________ Likelihood: "<<llHV<<" + "<<llX<<" = "<<(llHV+llX)<<" ______________________"<<std::endl;
    }
  } else {
    fLogLklhd.Total = llHV_NoMSc+llX_NoMSc;
    fLogLklhd.X = llX_NoMSc;
    fLogLklhd.HV = llHV_NoMSc;
    if( ND::tpcDebug().LikFit(DB_VVERBOSE) ) {
      std::cout<<" ________________ Likelihood: "<<llHV_NoMSc<<" + "<<llX_NoMSc<<" = "<<(llHV_NoMSc+llX_NoMSc)<<" ______________________"<<std::endl;
    }
  }

  return llHV+llX;
}



// *********************************************************************************
double TTPCLikFitPath::log_likelihoodHV(){


  double result = 0.0;

  unsigned int iclu = 0;

  double *helixState = new double[7];

  for (ND::THitSelection::const_iterator tmpClu = fFitClu->begin() ; tmpClu != fFitClu->end(); tmpClu++, iclu++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
    double yPred;
    double zPred;
    double yDirPred;
    double zDirPred;
    // PropagateToHVCluster uses CalibZ(Y) for the vertical(horizontal) clusters
    // so the field corrections are taken into account.
    bool ok = ND::helixPropagator().PropagateToHVCluster(Cluster);
    if(fCalculatorMode && !ok){
      result = 0.0;
      break;
    }
    double cluResult = 0.0;
  
    // We want to propagate from cluster to cluster but we don't want to use
    // all clusters for the likelihood calculation.
    if (! Cluster->isUsable()){
      continue;
    }

    ND::helixPropagator().GetHelixPosDirCurv(helixState);
    yPred    = helixState[1];
    zPred    = helixState[2];
    yDirPred = helixState[4];
    zDirPred = helixState[5];

    double cluResidual;
    // Apply the field correction like in tpcRecon on the predicted position
    // out of consistency and such that we do it once per cluster.
    if (Cluster->IsVertical()){
      yPred    += Cluster->GetDeltaY();
      cluResidual = fabs(Cluster->Y() - yPred);
    } else {
      zPred    += Cluster->GetDeltaZ();
      cluResidual = fabs(Cluster->Z() - zPred);
    }

    double drift = Cluster->GetDriftDistance();
    // fInitialSigma is the tranverse diffusion squared used as constant to calculate the change
    // of transverse diffusion across the track.
    double sigmaDiff = TMath::Sqrt(TMath::Abs( fInitialSigma*(drift-fMeanDrift) + fSigma*fMeanDrift ) );

    double Longitudinal;
    double Transversal;
    double phi;
    if (Cluster->IsVertical()){
      Longitudinal = fPadWidth;
      Transversal = fPadHeight;
      phi =  TMath::ATan2(yDirPred, zDirPred);
    } else {
      Longitudinal = fPadHeight;
      Transversal = fPadWidth;
      phi =  TMath::ATan2(zDirPred, yDirPred);
    }
    if ( fabs(phi) > TMath::PiOver2() ){
      if ( phi > 0.0)
        phi -= TMath::Pi();
      else
        phi += TMath::Pi();
    }

    double qtot = q_exp(0.0, phi, sigmaDiff, Longitudinal, 1000.*Transversal);

    if( qtot <= fMinPredIntCharge ) qtot = fMinPredIntCharge;

    unsigned int ipad = 0;
    const ND::THitSelection allPads = Cluster->GetHits();
    for (ND::THitSelection::const_iterator tmpHit = allPads.begin(); tmpHit != allPads.end(); ++tmpHit, ipad++) {
      ND::THandle<ND::TTPCHitPad> hitPad = *tmpHit;
      double yPad = hitPad->Y();
      double zPad = hitPad->Z();
      double b_col;
      if (Cluster->IsVertical()){
        b_col = yPred - yPad;
      } else {
        b_col = zPred - zPad;
      }
      double qval = q_exp(b_col, phi, sigmaDiff, Longitudinal, Transversal);

      if( qval <= fMinPredIntCharge ) qval = 0.0; 

      double I =  qval/qtot;

      if( I > 1.) {
        if( ND::tpcDebug().LikFit(DB_VVERBOSE) && I > 1.0001 ){
          std::cout<<" HV Likelihood I too large "<<I<<" local pad "<<qval<< " Total "<<qtot<<" in cluster "<<iclu<<", at pad (y,z) = ("<<yPad<<", "<<zPad<<")"<<std::endl;
          std::cout<<"    b_col = "<<b_col<<"   phi = "<<phi<<"   sigmaDiff = "<<sigmaDiff<<"    Longitudinal = "<<Longitudinal<<"    Transversal = "<<Transversal<<std::endl;
          std::cout<<"    drift = "<<drift<<"   fMeanDrift = "<<fMeanDrift<<"    fSigma = "<<fSigma<<std::endl;
        }
        I = 1.;
      }
      else if( I < 0. ) {
        if( ND::tpcDebug().LikFit(DB_VVERBOSE)) 
          std::cout<<" HV Likelihood I negative "<<I<<" in cluster "<<iclu<<", at pad (y,z) = ("<<yPad<<", "<<zPad<<")"<< std::endl;
        I = 0.;
      }

      double hitCharge = 0.0;
      if (hitPad->IsFitted())
        hitCharge = hitPad->ChargeFit();
      else
        hitCharge = hitPad->GetCharge();

      if( hitCharge > 0.0 ){
        double valLH = hitCharge*TMath::Log((I+fNoise)/(1.+allPads.size()*fNoise));

        if ( !finite(valLH) ){
          if ( ND::tpcDebug().LikFit(DB_VVERBOSE)){
            std::cout<<" likelihood would be nan or infinite"<<std::endl;
            std::cout<<"  plane "<<ipad<<" qval= "<<qval<<" qtot="<<qtot<<" bcol="<<b_col<<" phi="<<phi<<" sigmaDiff="<<sigmaDiff<<std::endl;
            std::cout<<"  I+noise="<<I+fNoise<<" npads="<<allPads.size()<<" 1+npads*noise="<<1.+allPads.size()*fNoise<<std::endl;
          }
        } else {
          cluResult -= valLH;
        }
      }
    }
    
    double Penalty = 1.0;
    if (fApplyResidualDependentPenalty){
      if (cluResidual > fPenaltyThreshold)
        Penalty = ( ( cluResidual - fPenaltyThreshold)*( cluResidual - fPenaltyThreshold)  * fPenaltyCoefficient )+1;
    }
    result += Penalty * cluResult;
  }

  // If the RecPack propagation failed, then something must
  // be wrong with this set of parameters.
  if( result == 0.0 || isnan(result) ) {

    if( isnan(result) && ND::tpcDebug().LikFit(DB_VVERBOSE))
      std::cout << " Horiz Likelihood NAN value " << std::endl; 

    result = 2.01e+21;
  }

  delete[] helixState;
  return result;
}

// *********************************************************************************
double TTPCLikFitPath::log_likelihoodX(){

  double result = 0.0;

  unsigned int iclu = 0;

  double *helixState = new double[7];
//  std::cout<<" =============================================="<<std::endl;
//  std::cout<<" =                   X fit                    ="<<std::endl;
//  std::cout<<" =============================================="<<std::endl;

  for (ND::THitSelection::const_iterator tmpClu = fFitClu->begin() ; tmpClu != fFitClu->end(); tmpClu++, iclu++) {
    ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
    double xPred;
    bool ok = ND::helixPropagator().PropagateToHVCluster(Cluster);
    if(fCalculatorMode && !ok){
      result = 0.0;
      break;
    }
  
    // We want to propagate from cluster to cluster but we don't want to use
    // all clusters for the likelihood calculation.
    if (! Cluster->isUsable())
      continue;

    ND::helixPropagator().GetHelixPosDirCurv(helixState);
    xPred = helixState[0];

    double xClu = Cluster->CalibX();
    result += (xPred - xClu) * (xPred - xClu) * Cluster->GetCharge();
// if(fCalculatorMode){
//   std::cout<<" ---->>> X => "<<xPred<<" - "<<xClu<<"  -> "<< (xPred - xClu) * (xPred - xClu) * Cluster->GetCharge()<<std::endl;
//   std::cout<<"         on cluster "<<Cluster->CalibX()<<", "<<Cluster->CalibY()<<", "<<Cluster->CalibZ()<<std::endl;
// }
  }

  if( result == 0.0 || isnan(result) ) result = 2.01e+21;

  delete[] helixState;
  return result;
}



// *********************************************************************************
void TTPCLikFitPath::StoreFittedState(){
  double x,ex;
  double y,ey;
  double z,ez;
  double tanx,etanx;
  double tanyorz,etanyorz;
  double curv,ecurv;

  fMinuit->GetParameter(XPARAM,x,ex);
  fMinuit->GetParameter(YPARAM,y,ey);
  fMinuit->GetParameter(ZPARAM,z,ez);
  fMinuit->GetParameter(TANXPARAM,tanx,etanx);
  fMinuit->GetParameter(TANYORZPARAM,tanyorz,etanyorz);
  fMinuit->GetParameter(CURVPARAM,curv,ecurv);

  EVector minuitVect = EVector(6,0);
  EMatrix minuitCova = EMatrix(6,6,0);
  EVector resultVect = EVector(7,0);
  EMatrix resultCova = EMatrix(7,7,0);

  minuitVect[0] = x;
  minuitVect[1] = y;
  minuitVect[2] = z;
  minuitVect[3] = tanx;
  minuitVect[4] = tanyorz;
  minuitVect[5] = curv;

  double MinVal,Edm,errdef;
  int    nfree,ntot,istat;
  fMinuit->mnstat(MinVal,Edm,errdef,nfree,ntot,istat);

  double covar[nfree][nfree];
  bzero(covar,nfree*nfree*sizeof(double));

  fMinuit->mnemat(&covar[0][0],nfree);

  // Conversion from minuit matrix of free parameters only into the RP matrix
  // of all parameters.
  int FPar[nfree];
  unsigned NParTot = 0; 

  FPar[0] = XPARAM; 
  if( fFitYPosParam ) 
    FPar[1] = YPARAM;
  else 
    FPar[1] = ZPARAM;
  FPar[2] = TANXPARAM;
  FPar[3] = TANYORZPARAM;
  FPar[4] = CURVPARAM;
  NParTot = 5;

  if( fFitYPosParam ) 
    minuitCova[ZPARAM][ZPARAM] = 0.00001; // Non zero to initialize the ZPARAM
  else 
    minuitCova[YPARAM][YPARAM] = 0.00001; // Non zero to initialize the YPARAM

  for (unsigned int i = 0; i < NParTot ; ++i){
    for (unsigned int j = 0; j <= i; ++j) {
      int l = FPar[i]; 
      int m = FPar[j]; 
      minuitCova[l][m] = covar[i][j];
      minuitCova[m][l] = minuitCova[l][m]; 
    }
  }

  ND::helixPropagator().PosTanCurvToPosDirQoP(minuitVect, minuitCova, resultVect, resultCova);
  State RPState;
  RPState.set_hv(HyperVector(resultVect, resultCova));
  RPState.set_hv(RP::sense, HyperVector(1,0));  
  RPState.set_name(RP::representation,RP::pos_dir_curv);
  fFitResults.FitState = RPState;
  fFitResults.IsFitReliable = fReliableFit;

  fFitResults.LogLikelihood = fLogLklhd;

  // To avoid recalculating it elsewhere
  fFitResults.Curvature = curv;
  fFitResults.eCurvature = ecurv;
}  


// *********************************************************************************
void TTPCLikFitPath::StoreFittedSigma(){
  // double sigma,esigma;
  fMinuit->GetParameter(SGMPARAM,fFitResults.Sigma, fFitResults.eSigma);
}


// *********************************************************************************
void TTPCLikFitPath::SaveFitResults(ND::THandle<ND::TTPCPath> Path){
  Path->SaveFitState(fFitResults);
  Path->SaveInRealDatum("likFit_Sigma", fFitResults.Sigma);
  Path->SaveInRealDatum("likFit_SigmaUnc", fFitResults.eSigma);
  Path->SaveInRealDatum("likFit_MeanDrift", fMeanDrift);

}


// *********************************************************************************
TTPCPathFitResults TTPCLikFitPath::GetFitResults(){
  return fFitResults;
}
