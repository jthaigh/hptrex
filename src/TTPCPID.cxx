#include "TTPCPID.hxx"

#include <list>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include <HEPUnits.hxx>
#include <HEPConstants.hxx>
#include <THandle.hxx>

#include <TOADatabase.hxx>
#include <TComboHit.hxx>
#include <TMultiHit.hxx>
#include <THandle.hxx>
#include <THit.hxx>
#include <TOARuntimeParameters.hxx>
#include <TReconPID.hxx>
#include <TRecPackManager.hxx>
#include <TIntegerDatum.hxx>
#include <TRealDatum.hxx>
#include <TrackingUtils.hxx>
#include <TG4Trajectory.hxx>

#include "TTPCRecPackUtils.hxx"
#include "TTPCUtils.hxx"
#include "TTPCHitPad.hxx"

//************************************************************
ND::TTPCPID::TTPCPID(void){
  std::cout << " TPC pid package initialized " << std::endl;
  Init();
}

//************************************************************
void ND::TTPCPID::Reset(void){
  fMomentum = -99;
  fMomentumErr = -99.;
  fTotCharge = -99;
  fDedxmeas = -99;
  fDedxcorr = -99;
  SLmean = -99;
  fNTrunMM = -99;

  fNSample = -99;

  for(int i=0;i<5;i++)
  {
    fPull[i]=-99;
    fSigmaexp[i]=-99;
    fDedxexp[i]=-99;
  }
}

//************************************************************
void ND::TTPCPID::Init(void){

  fVerbose =  ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.Verbosity");
  fUseHeightBasedCharge = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.UseHeightBasedCharge"); // Use height based charge instead of integrated charge for PID?
  fHeightToAreaFactor = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.HeightToAreaFactor"); // Factor to convert height based charge to area based integrated charge

  fMIPmuon = 0;

  fSigmaMIP = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SigmaMip"); //430; ///Expected sigma for MIPs

  fSLhor = ND::TGeomInfo::Get().TPC().GetPadXPitch()/10.;
  fSLver = ND::TGeomInfo::Get().TPC().GetPadYPitch()/10.;

  fPadFraction = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.TrunMean.PadFrac");

  fExpecteddE0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ExpecteddEdx_YZstates.P0");
  fExpecteddE1 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ExpecteddEdx_YZstates.P1");
  fExpecteddE2 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ExpecteddEdx_YZstates.P2");
  fExpecteddE3 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ExpecteddEdx_YZstates.P3");
  fExpecteddE4 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ExpecteddEdx_YZstates.P4");

  fScaleFactor = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ScaleFactor");

  fSLCorrH0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorr_YZstates.H0");
  fSLCorrH1 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorr_YZstates.H1");
  fSLCorrH2 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorr_YZstates.H2");

  fRowCorrF0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorr_YZstates.F0");
  fRowCorrF1 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorr_YZstates.F1");
  fRowCorrF2 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorr_YZstates.F2");
  fRowCorrF3 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorr_YZstates.F3");
  fRowCorrF4 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorr_YZstates.F4");
  fRowCorrF5 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorr_YZstates.F5");

  // Parameters for the sigma Sample Length correction
  fSigmaSLCorrA0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorrSigma_YZstates.A0");
  fSigmaSLCorrA1 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorrSigma_YZstates.A1");
  fSigmaSLCorrA2 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorrSigma_YZstates.A2");
  fSigmaSLCorrA3 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.SampleLengthCorrSigma_YZstates.A3");

  // Parameters for the sigma N sample correction
  fSigmaRowCorrP0 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorrSigma_YZstates.P0");
  fSigmaRowCorrP1 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorrSigma_YZstates.P1");
  fSigmaRowCorrP2 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorrSigma_YZstates.P2");
  fSigmaRowCorrP3 = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.RowCorrSigma_YZstates.P3");

  fElectronCorr = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.ElectronCorr"); //430; ///Expected sigma for MIPs

  //offset in unit of charge/distance to correct for the charge dependence on the drift dist
  fChargeVsDriftParameter = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.PID.DriftParameter");

  Reset();
}

//************************************************************
ND::TTPCPID::~TTPCPID() {}

//************************************************************
void ND::TTPCPID::Process(ND::THandle<ND::TTPCPattern> Pattern){

  Reset();

  std::vector< ND::THandle<ND::TTPCPath> > Paths = Pattern->GetPaths();
  // return false when no path found
  if (Paths.size()==0) return;

  for (std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin(); pth != Paths.end(); pth++) {
    // First reset the class variables
    Reset();
    
    if (fVerbose) {
      std::cout << "PID output: ============ new Path  ============== " << std::endl;
    }

    if ( (*pth)->CheckStatus(ND::TReconBase::kLikelihoodFit) &&
         (*pth)->CheckStatus(ND::TReconBase::kSuccess)){
      TrackdEdx(*pth);
      /*double complete = 0.0;
      double clean = 0.0;
      int G4Id;
      ND::THandle<ND::TG4TrajectoryContainer> trajectories;
      ND::THandle<ND::TG4Trajectory> G4Traj = TTPCUtils::FindTrueTrajectories((*pth)->GetHits(), complete, clean);
      if (G4Traj){
        int TruPdgId = G4Traj->GetPDGEncoding();
        if (abs(TruPdgId) == 11){
          (*pth)->SetPID(ND::TReconPID::kElectron, 0.0);
        } else if (abs(TruPdgId) == 13){
          (*pth)->SetPID(ND::TReconPID::kMuon, 0.0);
        } else if (abs(TruPdgId) == 211){
          (*pth)->SetPID(ND::TReconPID::kPion, 0.0);
        } else if (abs(TruPdgId) == 321){
          (*pth)->SetPID(ND::TReconPID::kKaon, 0.0);
        } else if (abs(TruPdgId) == 2212){
          (*pth)->SetPID(ND::TReconPID::kProton, 0.0);
        }
	}*/
    }
    if (fVerbose) {
      std::cout << "PID output: Done. " << std::endl;
    }
  }
}



//************************************************************
Double_t ND::TTPCPID::ExpecteddEdx(Double_t bg){
//************************************************************

  double beta=sqrt((bg*bg)/(1.+bg*bg));

  double func=fExpecteddE1-pow(beta,fExpecteddE3)-log(fExpecteddE2+1./pow(bg,fExpecteddE4));

  func=func*fExpecteddE0/pow(beta,fExpecteddE3);

  return func;
}

//************************************************************
Double_t ND::TTPCPID::TrackMomError(Double_t momentum, Double_t momerr, Double_t Mparticle){
//************************************************************

  double bgless, bgmore, dedxless, dedxmore;
  double epsilon = 2;

  if(momentum < 6)
    bgless=1;///protect bgless from being negative
  else
    bgless = (momentum-epsilon)/Mparticle;

  bgmore = (momentum+epsilon)/Mparticle;

  dedxless = ExpecteddEdx(bgless);
  dedxmore = ExpecteddEdx(bgmore);

  return TMath::Abs((dedxless-dedxmore)/(bgless-bgmore)*momerr/Mparticle);
}



//method to compute the dEdx with Horizontal clusters enabled

//************************************************************
void ND::TTPCPID::TrackdEdx(ND::THandle<ND::TTPCPath> Path){

  double momentum;
  double errmomentum;

  //I should move it to a reasonable number.. 288?
#define NSAMPLES 500

  double clcharge[NSAMPLES];
  double MMcharge[NSAMPLES];
  bool isGoodCluster[NSAMPLES];
  double SLpoint[NSAMPLES];

  double clchargecorr[NSAMPLES];
  double MMchargecorr[NSAMPLES];
  int npoint=0;
  int MMsample=0;
  double fDedxmeascorr=0;
  double sigmade;
  double Crow;
  double Csamplelength=0;
  double Csigmarow;
  double Csigmasamplelength;

  std::string cltype[NSAMPLES];

  fMomentum=0;
  fMomentumErr=0;
  fTotCharge =0;
  fDedxmeas=0;
  fDedxcorr=0;
  SLmean=0;

  fMIPmuon=ExpecteddEdx(4.5); //different from the original method because of the new parameterization is used

  for(int i=0;i<NSAMPLES;i++)
  {
    clcharge[i]=0;
    MMcharge[i]=0;
    SLpoint[i]=0;
    MMchargecorr[i]=0;
    clchargecorr[i]=0;
    isGoodCluster[i]=false;

    cltype[i]="VER";
  }

  ND::THandle<ND::THitSelection> Clusters = Path->GetHits();

  if(Clusters->size() > NSAMPLES) {
    std::cout<<"TTPCPID: WARNING. There are "<<Clusters->size()<<" clusters in this path. More than the "<<NSAMPLES<<" samples allowed. Skip PID."<<std::endl;
    return;
  }

  // Get the list of states for all the nodes, one to one.
  fNSample=0;

  ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();

  // Compute the needed info by hand from the clusters and track state.

  TTPCRecPackUtils::SavedParam RPParam = TTPCRecPackUtils::InitForPropagationInTPC();
  ND::rpman("TREx").model_svc().model().intersector().set_length_sign(0);

  // If problem with the fitted statem then just don't do the PID.
  if( !Path->HasFitState()){
    if (fVerbose) std::cout << " problem with the fit, skip the PID " << std::endl;
    return;
  }
  State propagState = Path->GetFitState();

  propagState.set_name(RP::representation,RP::pos_dir_curv);
  EVector vect = propagState.vector();

  momentum = fabs(1./vect[6]);
  fMomentum = momentum;
  if (momentum==0 || momentum<0)
  //{
    std::cout << "ERROR: momentum <= 0! EXIT PID because not possible. Momentum = " << momentum << std::endl;
  //  return;
  //}

  // A bit brutal just to get the error.
  // Maybe this could be done by copy pasting only the important parts of the conversion.
  ND::THandle<ND::TReconState> recoState(new ND::TTrackState);
  ND::converter().State_to_TReconState(propagState,(*recoState));
  errmomentum = TrackingUtils::GetMomentumError(recoState);
  fMomentumErr = errmomentum; 

  // Loop over the clusters
  for(ND::THitSelection::const_iterator itClu = Clusters->begin(); itClu!=Clusters->end(); itClu++){
    ND::THandle<ND::TTPCHVCluster> Clu = *itClu;
    // Calculate EDeposit

    double charge = Clu->GetIntCharge();  // Get the integral of the charge (as in tpcRecon) to compute the PID
    if (fUseHeightBasedCharge){
        charge = Clu->GetCharge();      // Use the height based charge instead
        charge *= fHeightToAreaFactor;  // Multiply with fudge factor to convert from signal heigth based charge, to area based integrated charge
    }

    if (!TTPCRecPackUtils::PropagateToHVCluster(Clu, propagState))
      continue;

    //vect[0] = x; vect[1] = y; vect[2] = z; vect[3] = TVector3.X();
    //vect[4] = TVector3.Y(); vect[5] = TVector3.Z(); vect[6] = q/p;
    vect = propagState.vector();

    //compute Sample Length on each cluster and charge corrected for the angle
    if( fNSample >= NSAMPLES ) continue;

    if( Clu->IsVertical() ) {

      cltype[fNSample] = "VER";

      SLpoint[fNSample] = fSLhor/fabs(vect[5]);
      clcharge[fNSample]= charge*fabs(vect[5])/fSLhor ;

    }

    else {

      cltype[fNSample] = "HOR";

      SLpoint[fNSample] = fSLver/fabs(vect[4]);
      clcharge[fNSample]=charge*fabs(vect[4])/fSLver ;

    }

    // TODO: TEMPORARY HACK. The ClusterSelection method should change.
    ND::THandle<ND::TComboHit> combo = Clu->ConvertToOAEvent();
    //
    //select good clusters
    //
    ClusterSelection(combo, tpcGeom, cltype, fNSample, isGoodCluster);



    //
    ///apply correction for the drift distance
    //
    // TODO DANGER: what if the pattern has hits on both sides of the cathode ???
    // fDriftSense should NOT be a global variable !!!!!!!!!
    clcharge[fNSample] += fChargeVsDriftParameter*fDriftSense*(vect[0]-fDriftCenter*fDriftSense);

    

    if(SLpoint[fNSample]<100.) SLmean += SLpoint[fNSample];

    //
    // correct the sample length
    Csamplelength = fSLCorrH0+ fSLCorrH1*SLpoint[fNSample] + fSLCorrH2*SLpoint[fNSample]*SLpoint[fNSample];

    if(SLpoint[fNSample]<2.8)
      clchargecorr[fNSample]=clcharge[fNSample]/Csamplelength;
    else
      clchargecorr[fNSample]=clcharge[fNSample]/1.15; //saturate correction

    ////apply the selection on the cluster to be ready to compute the truncated mean
    if(isGoodCluster[fNSample]){

      MMsample++;
      if (MMcharge[npoint]<0)
	std::cout << "Negative charge deposited in the MM!!  "<< std::endl;
      MMcharge[npoint]=clcharge[fNSample];
      MMchargecorr[npoint]=clchargecorr[fNSample];
      npoint++;
    }

    fNSample++;

  }

  // to collect the raw charge per path before any correction and truncated mean
  for(int i=0;i <= fNSample;i++)
    fTotCharge+=clcharge[i];


  TTPCRecPackUtils::ResetAfterPropagationInTPC(RPParam);

  SLmean=SLmean/fNSample;

  ///take the 70% of the cluster with less energy
  std::sort(MMcharge,MMcharge+MMsample);
  std::sort(MMchargecorr,MMchargecorr+MMsample);

  fNTrunMM = (int)(fPadFraction * float(MMsample));

  for(int i=0;i <= fNTrunMM;i++){
    fDedxmeas+=MMcharge[i];
    fDedxmeascorr+=MMchargecorr[i];
  }


  fDedxmeas=fDedxmeas/(double)fNTrunMM;
  fDedxmeascorr=fDedxmeascorr/(double)fNTrunMM;

  ////////////compute row corrections

  //int trackpoint = fNSample;
  int trackpoint = MMsample; // this is in prod6: bug. Corrections should be applied considering the full number of clusters because are estimated this way

  if( trackpoint <= 72 ){
    Crow = fRowCorrF0+fRowCorrF1*trackpoint + fRowCorrF2*trackpoint*trackpoint;
  }
  else{
    Crow = fRowCorrF3 + (fRowCorrF4)*trackpoint+fRowCorrF5*trackpoint*trackpoint;;
  }
  fDedxcorr=fDedxmeascorr/Crow;


  ////compute and apply corrections to the sigma

  Csigmasamplelength= ( exp(fSigmaSLCorrA0+SLmean*fSigmaSLCorrA1) + exp(fSigmaSLCorrA2+SLmean*fSigmaSLCorrA3) )  ;
  Csigmarow= ( exp(fSigmaRowCorrP0+fSigmaRowCorrP1*trackpoint) + exp(fSigmaRowCorrP2+fSigmaRowCorrP3*trackpoint) );

  // check if the fScaleFactor has to be really applied?
  sigmade=fSigmaMIP*Csigmarow*Csigmasamplelength*fScaleFactor;

  if(fVerbose){
    std::cout<< "PIDoutput: --- raw charge     : " << fTotCharge << std::endl;
    std::cout<< "PIDoutput: --- dEdx meas      : " << fDedxmeas << std::endl;
    std::cout<< "PIDoutput: --- Sample Length  : " << Csamplelength<<std::endl;
    std::cout<< "PIDoutput: --- sigma Nsamples : " << Csigmasamplelength<<std::endl;
    std::cout<< "PIDoutput: --- Nb of samples  : " << fNSample <<std::endl;
    std::cout<< "PIDoutput: --- 70% of samples : " << fNTrunMM <<std::endl;
    std::cout<< "PIDoutput: --- momentum       : " << momentum  << std::endl;
    std::cout<< "PIDoutput: --- momentum error : " << errmomentum  << std::endl;
    //std::cout<< "PIDoutput: --- row corr. factor : " << Crow  << std::endl;
    //std::cout<< "PIDoutput: --- sigma row dEdx : " << Csigmarow <<std::endl;

    std::cout<< "PIDoutput: --- dEdx meascorr  : " << fDedxcorr << std::endl;

  }

  double Mparticle=0;

  for(int i=0;i<5;i++){
    switch( i ) {
    case 0: //electron
        Mparticle=0.511;
        break;
    case 1: //muon
        Mparticle=105.66;
        break;
    case 2: //pion
        Mparticle=139.57;
        break;
    case 3: //kaon
        Mparticle=493.67;
        break;
    case 4: //proton
        Mparticle=938.27;
        break;
    }

    double bg=momentum/Mparticle;
    double sigmamom;
    double sigmadedx;

    //expected energy loss and sigma for different particles
    fDedxexp[i]= ExpecteddEdx(bg);

    if (fDedxexp[i] < 0)
    {
      if (fVerbose){
        std::cout<< "NEGATIVE dedx value. Pull setted at -9999 " <<std::endl;
	std::cout << "Particle type "<< i << " dEdx_exp " << fDedxexp[i] <<" mom " << momentum << " dEdx_meas " << fDedxmeas << std::endl;
      }

      fDedxexp[i]   = -9999.;
      sigmamom     = -9999.;
      sigmadedx    = -9999.;
      fPull[i]      = -9999.;
    }
    else if (fDedxexp[i] > 1e07)
    {
      if (fVerbose){
        std::cout<< "INFINITE dedx value. Pull setted at -9999 " <<std::endl;
	std::cout << "Particle type "<< i << " dEdx_exp " << fDedxexp[i] <<" mom " << momentum << " dEdx_meas " << fDedxmeas << std::endl;
      }

      fDedxexp[i]   = -9999.;
      sigmamom     = -9999.;
      sigmadedx    = -9999.;
      fPull[i]      = -9999.;
    }
    else{
      ///small correction to the dE/dx of the electrons as in the data it's different from the expected Geant4
      if(i==0)
        fDedxexp[i]=fDedxexp[i]*fElectronCorr; ///tuned with Cosmics

      sigmamom=TrackMomError(momentum,errmomentum,Mparticle);

      sigmadedx=sigmade*TMath::Sqrt(fDedxexp[i]/fMIPmuon);
      if(fDedxexp[i]<0) sigmadedx=sigmade;

      fSigmaexp[i]=TMath::Sqrt(pow(sigmadedx,2)+pow(sigmamom,2.0));
      if( isnan(sigmamom) ) fSigmaexp[i]=sigmadedx;


      if( fSigmaexp[i] == 0. ) std::cout << " TTPCPID: Error in sigma " << std::endl;

      ///computation of the pull
      fPull[i]=(fDedxcorr-fDedxexp[i])/fSigmaexp[i];

      if (fVerbose){
	// only to debug a precise event
	//std::cout<< "PIDoutput: --- sigma de       : " << sigmade       <<std::endl;
	//std::cout<< "PIDoutput: --- sigma momentum : " << sigmamom      <<std::endl;
	//std::cout<< "PIDoutput: --- sigma dEdx     : " << sigmadedx     <<std::endl;
	//std::cout<< "PIDoutput: --- sigma expdEdx  : " << fSigmaexp[i]  <<std::endl;
      }
    }
  }

  if (fVerbose){
    std::cout<< "PIDoutput: --- Pull hypothesis: " << std::endl;
    std::cout<< "PIDoutput: ---                : - Muon     : " << fPull[1] << std::endl;
    std::cout<< "PIDoutput: ---                : - Electron : " << fPull[0] << std::endl;
    std::cout<< "PIDoutput: ---                : - Proton   : " << fPull[4] << std::endl;
    std::cout<< "PIDoutput: ---                : - Pion     : " << fPull[2] << std::endl;
    std::cout<< "PIDoutput: ---                : - Kaon     : " << fPull[3] << std::endl;
  }

  WritePiDInfo(Path);

  return;
}// end of dEdxwHorClu

//************************************************************
void ND::TTPCPID::ClusterSelection(ND::THandle<ND::TComboHit>& m_combo, ND::TTPCGeom& m_tpcGeom, std::string *m_cltype, int &m_nSample, bool *isGoodCluster){

  // int xpos[500];
  int xpos[NSAMPLES];

  int rowflag[NSAMPLES];
  int columnflag[NSAMPLES];
  int xposflag[NSAMPLES];

  //initialize some variables
  //for (int i =0; i< 500; i++){
  for (int i =0; i< NSAMPLES; i++){
    xpos[i]=0;
    xposflag[i]=0;
    rowflag[i]=0;
    columnflag[i]=0;
  }

  for(ND::THitSelection::const_iterator Hit = m_combo->GetHits().begin() ;Hit!=m_combo->GetHits().end();Hit++)
  {  ///check if the cluster is on the board of 1 MM module or on the cathod
    fDriftCenter=m_tpcGeom.GetMaxDriftDistance((*Hit)->GetGeomId())/2;
    fDriftSense= (double)m_tpcGeom.GetDriftSense((*Hit)->GetGeomId());

    // flagging vertical clusters crossing the horizontal edge
    //
    if (m_cltype[m_nSample] == "VER"){
      if (m_tpcGeom.PadIsAtHorEdge((*Hit)->GetGeomId() ) )
        rowflag[m_nSample]=1;
    }
    else if (m_cltype[m_nSample] == "HOR"){
      if (m_tpcGeom.PadIsAtVerEdge((*Hit)->GetGeomId() ) )
        columnflag[m_nSample] = 1;
    }

    TVector3 positionchannel;
    if(Hit==m_combo->GetHits().begin())
    {
      m_tpcGeom.GeomIdToGlobalXYZ((*Hit)->GetGeomId(),positionchannel);
      if(positionchannel.X()>0)
        xpos[m_nSample]=1;
      else
        xpos[m_nSample]=-1;
    }

  }//end combo->getHits

  if(fNSample>0) {
    if(xpos[fNSample]*xpos[fNSample-1]<0) {
      xposflag[fNSample-1]=1;
      xposflag[fNSample]=1;
    }
  }

  if ( xposflag[fNSample]==0 && rowflag[fNSample]==0  && columnflag[fNSample]==0)
    isGoodCluster[fNSample]= true;
} // end of ClusterSelection



//************************************************************
void ND::TTPCPID::WritePiDInfo(ND::THandle<ND::TTPCPath> path){

  /// Write the different values for the PiD
  path->SaveInRealDatum("tpcPid_Momentum", fMomentum);
  path->SaveInRealDatum("tpcPid_MomentumErr", fMomentumErr);
  path->SaveInRealDatum("tpcPid_TotCharge", fTotCharge);
  path->SaveInRealDatum("tpcPid_Craw", fDedxmeas);
  path->SaveInRealDatum("tpcPid_Ccorr", fDedxcorr);
  path->SaveInRealDatum("tpcPid_Nnodes", fNSample);
  path->SaveInRealDatum("tpcPid_NTrun", fNTrunMM);
  path->SaveInRealDatum("tpcPid_SampleLength", SLmean);
  path->SaveInRealDatum("tpcPid_dEdxexpEle",    fDedxexp[0]);
  path->SaveInRealDatum("tpcPid_SigmaEle",      fSigmaexp[0]);
  path->SaveInRealDatum("tpcPid_PullEle",       fPull[0]);
  path->SaveInRealDatum("tpcPid_dEdxexpMuon",   fDedxexp[1]);
  path->SaveInRealDatum("tpcPid_SigmaMuon",     fSigmaexp[1]);
  path->SaveInRealDatum("tpcPid_PullMuon",      fPull[1]);
  path->SaveInRealDatum("tpcPid_dEdxexpPion",   fDedxexp[2]);
  path->SaveInRealDatum("tpcPid_SigmaPion",     fSigmaexp[2]);
  path->SaveInRealDatum("tpcPid_PullPion",      fPull[2]);
  path->SaveInRealDatum("tpcPid_dEdxexpKaon",   fDedxexp[3]);
  path->SaveInRealDatum("tpcPid_SigmaKaon",     fSigmaexp[3]);
  path->SaveInRealDatum("tpcPid_PullKaon",      fPull[3]);
  path->SaveInRealDatum("tpcPid_dEdxexpProton", fDedxexp[4]);
  path->SaveInRealDatum("tpcPid_SigmaProton",   fSigmaexp[4]);
  path->SaveInRealDatum("tpcPid_PullProton",    fPull[4]);

  // list of considered particle ID hypothesis
  ND::TReconPID::ParticleId hypo[5] = {ND::TReconPID::kElectron,
    ND::TReconPID::kMuon,
    ND::TReconPID::kPion,
    ND::TReconPID::kKaon,
    ND::TReconPID::kProton};

  // compute the minimum pull
  double pull_min=1e10;
  int i_min=-3;
  for (int i=0;i<5;i++){
    if (fabs(fPull[i]) < pull_min){
      pull_min = fabs(fPull[i]);
      i_min = i;
    }
  }

  if (i_min >= 0){
    // Set the main hypothesis
    path->SetPID(hypo[i_min], GetLikelihood(fPull[i_min]));
  }
  else{
    path->SetPID(ND::TReconPID::kOther, 0.);
  }
}

//************************************************************
double ND::TTPCPID::GetLikelihood(double x){

  // assume a gaussian distribution of mean 0 and sigma 1
  double gauss = 1/sqrt(2*unit::pi)*exp(-(x*x)/2);

  return gauss;
}
