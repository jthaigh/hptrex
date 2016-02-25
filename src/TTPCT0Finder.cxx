#include "TTPCT0Finder.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCCalibration.hxx"
#include "TTPCUtils.hxx"
#include "TTPCDebug.hxx"

//#include <vector>
//#include <cmath>
#include <list>

#include <THandle.hxx>
#include <THit.hxx>
#include <TOARuntimeParameters.hxx>
#include <TRecPackManager.hxx>
#include <TGeomInfo.hxx>
#include <TVector3.h>
#include <TND280Event.hxx>
#include <TEventFolder.hxx>
//#include <recpack/KalmanFitter.h>
#include <HEPUnits.hxx>
//#include <THitP0DLastPlaneFilter.hxx>


TH1F* ND::TTPCT0Finder::hFGDt0 = NULL;
TH1F* ND::TTPCT0Finder::hECALt0 = NULL;
TH1F* ND::TTPCT0Finder::hP0Dt0 = NULL;
TH1F* ND::TTPCT0Finder::hSMRDt0 = NULL;
TH1F* ND::TTPCT0Finder::hChi2t0 = NULL;
TH1F* ND::TTPCT0Finder::hChi2t0nolog = NULL;
TH1F* ND::TTPCT0Finder::hUsedt0 = NULL;
TH1F* ND::TTPCT0Finder::hCathodet0 = NULL;
TH1F* ND::TTPCT0Finder::hUsedt0From = NULL;

TH1F* ND::TTPCT0Finder::hFGDmP0Dt0 = NULL;
TH1F* ND::TTPCT0Finder::hFGDmECALt0 = NULL;
TH1F* ND::TTPCT0Finder::hP0DmECALt0 = NULL;

TH2F* ND::TTPCT0Finder::hPositiontpc1t0 = NULL;
TH2F* ND::TTPCT0Finder::hPositiontpc2t0 = NULL;
TH2F* ND::TTPCT0Finder::hPositiontpc3t0 = NULL;

//*****************************************************************************
ND::TTPCT0Finder::TTPCT0Finder(void){
  enabledECAL = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.ECAL");
  enabledFGD = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.FGD");
  enabledP0D = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.P0D");
  enabledSMRD = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.SMRD");
  fMax_Drift = ND::TGeomInfo::Get().TPC().GetMaxDriftDistance(0,0);

  fFGD_Qmin  =  ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.FGDQmin");
  fECAL_Qmin =  ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.ECALQmin");
  fP0D_Qmin  =  ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.P0DQmin");
  fSMRD_Qmin  =  ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.SMRDQmin");

  fFGD_YedgeSize = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.FGDYedgeSize");
  fFGD_Layers    = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.FGDLayers");
  fECAL_Layers   = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.ECALLayers");
  fP0D_Layers   = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.P0DLayers");

  fChi2_Max  =  ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.Chi2Max");
  fSafetyResidualCut = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.Reco.T0.SafetyResidualCut");

  if (fEnableHistos==true && hFGDt0==NULL){
    hFGDt0  = new TH1F("hFGDt0" ,"t0 found by FGD match" ,300,-15000.0,15000.0);
    hECALt0 = new TH1F("hECALt0","t0 found by ECAL match",300,-15000.0,15000.0);
    hP0Dt0  = new TH1F("hP0Dt0" ,"t0 found by P0D match" ,300,-15000.0,15000.0);
    hSMRDt0  = new TH1F("hSMRDt0" ,"t0 found by SMRD match" ,300,-15000.0,15000.0);
    hChi2t0 = new TH1F("hChi2t0","log(Chi2) of hit used for t0",300,-50.0,10.0);
    hChi2t0nolog = new TH1F("hChi2t0nolog","Chi2 of hit used for t0",300,0.0,1.0);
    hUsedt0 = new TH1F("hUsedt0","t0 recommended by GetTimeFromSed",300,-15000.0,15000.0);
    hCathodet0  = new TH1F("hCathodet0","t0 recommended by fitT0fromCrossing",300,-15000.0,15000.0);
    hUsedt0From = new TH1F("hUsedt0From","GetTimeFromSeed t0 from FGD=1,ECAL=2,P0D=3",5,-0.5,4.5);

    hFGDmP0Dt0  = new TH1F("hFGDmP0Dt0" ,"FGD-P0D t0" ,200,-2000.0,2000.0);
    hFGDmECALt0 = new TH1F("hFGDmECALt0","FGD-ECAL t0",200,-2000.0,2000.0);
    hP0DmECALt0 = new TH1F("hP0DmECALt0","P0D-ECAL t0",200,-2000.0,2000.0);

    hPositiontpc1t0 = new TH2F("hPositiontpc1t0","Z,Y Position of hit used for TPC1 t0",300,-2500.0,3500.0,300,-2500.0,2500.0);
    hPositiontpc2t0 = new TH2F("hPositiontpc2t0","Z,Y Position of hit used for TPC2 t0",300,-2500.0,3500.0,300,-2500.0,2500.0);
    hPositiontpc3t0 = new TH2F("hPositiontpc3t0","Z,Y Position of hit used for TPC3 t0",300,-2500.0,3500.0,300,-2500.0,2500.0);

  }
}

//*****************************************************************************
ND::TTPCT0Finder::~TTPCT0Finder() {}


//*****************************************************************************
void ND::TTPCT0Finder::PrepareScintHits(ND::TND280Event& event){
  // This is event dependent so it has to be set for every event.
  fDriftVelocity = ND::tpcCalibration().GetDriftVelocity();

  if( ND::tpcDebug().T0Finder(DB_VERBOSE)   )
    std::cout << "TTPCT0Finder: Preparing hits in scintillator detectors." << std::endl;
  prepareFGDHits(event);
  prepareECALHits(event);
  prepareP0DHits(event);
  prepareSMRDHits(event);

}


//*****************************************************************************
void ND::TTPCT0Finder::CalculateDefaultT0(ND::TReconObjectContainer *allPatterns){
  // This is event dependent so it has to be set for every event.
  fDriftVelocity = ND::tpcCalibration().GetDriftVelocity();
  double minTime = +99999999.;
  double maxTime = -99999999.;
  for (ND::TReconObjectContainer::iterator pattit = allPatterns->begin(); pattit != allPatterns->end(); pattit++) {
    ND::THandle<ND::TTPCPattern> pattern = *pattit;
    double tmpMin = +99999999.;
    double tmpMax = -99999999.;
    FindTimeRange(pattern, tmpMin, tmpMax);
    minTime = std::min(minTime, tmpMin);
    maxTime = std::max(maxTime, tmpMax);
  }
  double anodet0   = minTime - ND::tpcCalibration().GetTimeOffset();
  double cathodet0 = maxTime - ND::tpcCalibration().GetTimeOffset() - ( fMax_Drift / fDriftVelocity );
  double DefaultT0 = 0.0;
  if ( cathodet0 < anodet0)
    DefaultT0 = (cathodet0 + anodet0) / 2.;
  else
    DefaultT0 = anodet0;

  ND::tpcCalibration().SetDefaultT0(DefaultT0);

}


//*****************************************************************************
void ND::TTPCT0Finder::Process(ND::THandle<ND::TTPCPattern> Pattern){
  // set default t0 offset
  TTPCT0 T0Result(ND::tpcCalibration().GetDefaultT0());
  
  // Get Paths
  std::vector< ND::THandle<ND::TTPCPath> > Paths = Pattern->GetPaths();
  if ( ND::tpcDebug().T0Finder(DB_INFO))
    std::cout<<" -- TTPCT0Finder: Search T0 for pattern with "<<Paths.size()<<" paths."<<std::endl;

  FindT0Range(Pattern, T0Result);
  double MinT0, MaxT0;
  T0Result.GetT0Range(MinT0, MaxT0);

  // Use the spread in time of the whole pattern to establish a range of T0.
  bool isT0Found = false;
  // If there are no hits available in the other detectors, don't bother.
  if( fFGD_TPC1.size() || fFGD_TPC2.size() || fFGD_TPC3.size() ||
      fECAL_TPC1.size() || fECAL_TPC2.size() || fECAL_TPC3.size() || fECAL_SIDE.size() ||
      fP0D.size() || fSMRD.size() ){

    // Order from longest to shortest the paths
    std::stable_sort(Paths.begin(), Paths.end(), TTPCUtils::SortShortestToLongest);

    // Use pop and go through the paths until one is assigned a T0
    for (std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin(); pth != Paths.end(); pth++) {
      ND::THandle<ND::TTPCPath> path = (*pth);
      if ( ND::tpcDebug().T0Finder(DB_INFO))
        std::cout<<" Search using the seed for Path Id "<<path->GetId()<<std::endl;
      // Call T0 for this path
      isT0Found = FindT0UsingSeed(path, T0Result);
      if (isT0Found)
        break;
    }
  }
  // Propagate the T0Result to all the constituents, even if we have no T0
  // just to set the guessed value depending on the position of the hits in X.
  Pattern->SetT0(T0Result);

}



//*****************************************************************************
bool ND::TTPCT0Finder::FindT0UsingSeed(ND::THandle<ND::TTPCPath> Path, TTPCT0 &T0found){
  if( !( ( fFGD_TPC1.size()+fFGD_TPC2.size()+fFGD_TPC3.size() != 0 && enabledFGD ) ||
      ( fECAL_TPC1.size()+fECAL_TPC2.size()+fECAL_TPC3.size()+fECAL_SIDE.size() != 0 && enabledECAL ) ||
      ( fP0D.size() != 0 && enabledP0D ) ||
      ( fSMRD.size() != 0 && enabledSMRD ) ) ) return false;

  bool ok = Path->HasSeedState();
  if (!ok){
    if( ND::tpcDebug().T0Finder(DB_ERROR) )
      std::cout << "  Fail to get seed state to compute T0  !!! " << std::endl;
    return false;
  }
  State seed = Path->GetFrontSeedState();

  // First which TPC are the hits in?
  int itpc=-1;
  if (Path->UsesDetector(ND::TReconBase::kTPC1)) itpc=1;
  else if (Path->UsesDetector(ND::TReconBase::kTPC2)) itpc=2;
  else if (Path->UsesDetector(ND::TReconBase::kTPC3)) itpc=3;

  TTPCT0 t0fgd(T0found);
  TTPCT0 t0ecal(T0found);
  TTPCT0 t0p0d(T0found);
  TTPCT0 t0smrd(T0found);

  if ( itpc==1 ){
    // use P0D, fg1, or barrel ecal hits
    if ( ND::tpcDebug().T0Finder(DB_VERBOSE)) std::cout<<" The path is in TPC1"<<std::endl;
    fitT0fromFGD( fFGD_TPC1, seed, t0fgd);
    if ( (!t0fgd.OneHitFound()) ){
      fitT0fromECAL( fECAL_TPC1,seed, t0ecal);
      fitT0fromP0D( fP0D, seed, t0p0d ); 
    }
  } else if (itpc==2) {
    // use fgd1, fgd2, or barrel ecal hits
    if ( ND::tpcDebug().T0Finder(DB_VERBOSE)) std::cout<<" The path is in TPC2"<<std::endl;
    fitT0fromFGD( fFGD_TPC2, seed, t0fgd ); 
    if ( (!t0fgd.OneHitFound()) ) {
    //if (1){
      fitT0fromECAL( fECAL_TPC2,seed, t0ecal); 
    //std::cout<<" ======>>> "<<t0fgd.GetChi2()<<"    "<<t0ecal.GetChi2()<<std::endl;
    }
  } else if (itpc==3){
    // use fgd2, barrel or downstream ecal hits
    if ( ND::tpcDebug().T0Finder(DB_VERBOSE)) std::cout<<" The path is in TPC3"<<std::endl;
    fitT0fromFGD( fFGD_TPC3, seed, t0fgd ); 
    if ( (!t0fgd.OneHitFound()) ){
    //if (1){
      fitT0fromECAL( fECAL_TPC3,seed, t0ecal); 
    //std::cout<<" ======>>> "<<t0fgd.GetChi2()<<"    "<<t0ecal.GetChi2()<<std::endl;
    }
  }

  // check if failed..., if so check for special cases
  if ( (!t0fgd.OneHitFound())  &&  (!t0ecal.OneHitFound())  &&  (!t0p0d.OneHitFound()) ){

    // check if could have gone around the FGD?
    double fgdymax = ND::TGeomInfo::Get().FGD().FGD1ActiveMax()[1]-50.0;
    double fgdymin = ND::TGeomInfo::Get().FGD().FGD1ActiveMin()[1]+50.0;

    if ( ND::tpcDebug().T0Finder(DB_VERBOSE))
      std::cout<<" First pass t0 not found.  y0="<<seed.vector()[1]<<" fgdymin="<<fgdymin<<" fgdymax="<<fgdymax<<std::endl;

    if ( seed.vector()[1] > fgdymax || seed.vector()[1] < fgdymin ){
      if ( ND::tpcDebug().T0Finder(DB_VERBOSE)) std::cout<<" -> Trying again"<<std::endl;

      // for tpc1, try using tpc2 ecal, tpc2 fgd, or tpc3 ecal
      if ( itpc==1 ){
        fitT0fromFGD( fFGD_TPC2, seed, t0fgd ); 
        if ( (!t0fgd.OneHitFound()) ){
          fitT0fromECAL( fECAL_TPC2,seed, t0ecal);
          if ( (!t0fgd.OneHitFound())  &&  (!t0ecal.OneHitFound()) ){
            fitT0fromECAL( fECAL_TPC3,seed, t0ecal); 
            if ( !t0ecal.OneHitFound() ){
              fitT0fromECAL( fECAL_SIDE,seed, t0ecal); 
            }
          }
        }
      }

      // for tpc2, try using tpc1 P0D or tpc1 ecal, or tpc3 ecal
      if ( itpc==2 ){
        fitT0fromECAL( fECAL_TPC3,seed, t0ecal); 
        fitT0fromP0D( fP0D, seed, t0p0d ); 
        if (  (!t0p0d.OneHitFound())  &&  (!t0ecal.OneHitFound()) ){
          fitT0fromECAL( fECAL_TPC1,seed, t0ecal);
          if ( !t0ecal.OneHitFound() ){
            fitT0fromECAL( fECAL_SIDE,seed, t0ecal); 
          }
        }
      }

      // for tpc3, try tpc2 fgd, tpc2 ecal, P0D
      if (itpc==3){
        fitT0fromFGD( fFGD_TPC2, seed, t0fgd ); 
        if ( (!t0fgd.OneHitFound()) ){
          fitT0fromECAL( fECAL_TPC2,seed, t0ecal); 
          fitT0fromP0D( fP0D, seed, t0p0d ); 
          if ( (!t0fgd.OneHitFound())  &&  (!t0ecal.OneHitFound())  &&  (!t0p0d.OneHitFound()) ){
            fitT0fromECAL( fECAL_TPC1,seed, t0ecal); 
            if ( !t0ecal.OneHitFound() ){
              fitT0fromECAL( fECAL_SIDE,seed, t0ecal); 
            }
          }
        }
      }
    }
  }


  // Only try SMRD if failed to get it from any other detector...
  if ( (!t0fgd.OneHitFound())  &&  (!t0ecal.OneHitFound())  &&  (!t0p0d.OneHitFound()) ){
    fitT0fromSMRD( fSMRD, seed, t0smrd ); 
  }

  if (fEnableHistos==true){
    if ( t0fgd.OneHitFound() && t0p0d.OneHitFound() ) hFGDmP0Dt0->Fill( t0fgd.GetT0()-t0p0d.GetT0() );
    if ( t0fgd.OneHitFound() && t0ecal.OneHitFound() ) hFGDmECALt0->Fill( t0fgd.GetT0()-t0ecal.GetT0() );
    if ( t0ecal.OneHitFound() && t0p0d.OneHitFound() ) hP0DmECALt0->Fill( t0p0d.GetT0()-t0ecal.GetT0() );
  }


  if( t0fgd.OneHitFound() ){
    T0found = t0fgd;
    FillT0Histos(1, itpc, t0fgd);
  } else if( t0ecal.OneHitFound()  && t0ecal.GetChi2() <= t0p0d.GetChi2() ) {
    T0found = t0ecal;
    FillT0Histos(2, itpc, t0ecal);
  } else if( t0p0d.OneHitFound() && t0p0d.GetChi2() <= t0ecal.GetChi2() ) {
    T0found = t0p0d;
    FillT0Histos(3, itpc, t0p0d);
  } else if( t0smrd.OneHitFound() ) {
    T0found = t0smrd;
    FillT0Histos(4, itpc, t0smrd);
  } else {
    if (fEnableHistos) hUsedt0From->Fill(0.0);
  }

  if (fEnableHistos==true) hUsedt0->Fill( T0found.GetT0() );

  if( T0found.GetSource() != kNoT0src ){
    if( ND::tpcDebug().T0Finder(DB_VERBOSE) )
      std::cout << "  Found T0 " << T0found.GetT0() << " ns ( FGD " << (t0fgd.GetT0()) << ", ECAL "  << (t0ecal.GetT0()) << ", P0D " << (t0p0d.GetT0()) << ")" << std::endl;

    return true;
  }
  return false;
}


//*****************************************************************************
void ND::TTPCT0Finder::FindTimeRange(ND::THandle<ND::TTPCPattern> Pattern, double &minTime, double &maxTime){
  // TODO: We need to come up with the rules to make this as efficient as possible
  // Find the most extreme hits
  for (ND::TReconObjectContainer::iterator constit = Pattern->GetConstituents()->begin(); constit != Pattern->GetConstituents()->end(); constit++) {
    ND::THandle<ND::THitSelection> Hits = (*constit)->GetHits();
    for (ND::THitSelection::iterator hitit = Hits->begin(); hitit != Hits->end(); ++hitit){
      ND::THandle<ND::TTPCHitPad> Hit = *hitit;
      if(Hit){
        if (Hit->IsFitted()){
          minTime = std::min(minTime, Hit->TimeFit());
          maxTime = std::max(maxTime, Hit->TimeFit());
        } else {
          minTime = std::min(minTime, Hit->GetTime());
          maxTime = std::max(maxTime, Hit->GetTime());
        }
      } else {
        ND::THandle<ND::TTPCHVCluster> Clu = *hitit;
        if (!Clu)
          continue;
        minTime = std::min(minTime, Clu->GetTime());
        maxTime = std::max(maxTime, Clu->GetTime());
      }
    }
  }
}


//*****************************************************************************
bool ND::TTPCT0Finder::FindT0Range(ND::THandle<ND::TTPCPattern> Pattern, TTPCT0 &T0found){
  // TODO: We need to come up with the rules to make this as efficient as possible
  // Find the most extreme hits
  double minTime = +999999999.;
  double maxTime = -999999999.;
  FindTimeRange(Pattern, minTime, maxTime);
  double atAnodeT0   = minTime - ND::tpcCalibration().GetTimeOffset();
  double atCathodeT0 = maxTime - ND::tpcCalibration().GetTimeOffset() - ( fMax_Drift / fDriftVelocity );

  T0found.LoadTimeRange(atCathodeT0, atAnodeT0);
  
  return true;
}


//*****************************************************************************
void ND::TTPCT0Finder::FindCathodeCrosserT0(ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPath> PathB, TTPCT0 &T0found){
  ND::THandle<ND::TTPCPath> Paths[2];
  Paths[0] = PathA;
  Paths[1] = PathB;

  // Start by calculating the cathode T0 no matter what.
  TTPCT0 CathT0;
  double maxHitTime = -999999999.;
  double maxCluTime = -999999999.;
  ND::THandle<ND::THit> maxHit;
  ND::THandle<ND::THit> maxClu;
  for (unsigned int p = 0; p < 2; p++){
    ND::THandle<ND::THitSelection> Hits = Paths[p]->GetHits();
    // Use the hit pads and the clusters because delta ray hits can mess up the mean time calculated in the clusters.
    for (ND::THitSelection::iterator cluIt = Hits->begin(); cluIt != Hits->end(); ++cluIt){
      ND::THandle<ND::TTPCHVCluster> Cluster = (*cluIt);
      
      maxCluTime = std::max(maxCluTime, Cluster->GetTime());
      maxClu = Cluster;
      for (ND::THitSelection::const_iterator hitIt = Cluster->GetHits().begin(); hitIt != Cluster->GetHits().end(); ++hitIt){
        maxHitTime = std::max(maxHitTime, (*hitIt)->GetTime());
        maxHit = (*hitIt);
      }
    }
  }
  // 150 ns corresponds to ~11mm which is taken only from the spread of X positions
  // when using the hit pads time only for cathode crosser T0 (See validation report HEAD_2015_11_17 on t2k.org).
  // When the maxHitTime is much greater than the maxCluTime, then we suspect that 
  // the mean time of the cluster is biased due to delta ray hits.
  double atCathodeT0;
  if ( (maxHitTime - maxCluTime) < 150.){
    atCathodeT0 = maxCluTime - ND::tpcCalibration().GetTimeOffset() - ( fMax_Drift / fDriftVelocity );
    CathT0.LoadCathodeHit(maxClu, atCathodeT0);
  } else {
    atCathodeT0 = maxHitTime - ND::tpcCalibration().GetTimeOffset() - ( fMax_Drift / fDriftVelocity );
    CathT0.LoadCathodeHit(maxHit, atCathodeT0);
  }

  unsigned int GoodT0 = 0;
  ND::THandle<ND::TTPCPath> GoodPath;
  for (unsigned int p = 0; p < 2; p++)
    if ( Paths[p]->GetT0Source() != kNoT0src && Paths[p]->GetT0Source() != kNegXBrECALT0src && Paths[p]->GetT0Source() != kPosXBrECALT0src){
      GoodT0++;
      GoodPath = Paths[p];
    }

  // The cathode T0 is the best thing we have here.
  if ( GoodT0 == 0){
    T0found = CathT0;
  }
  // If there is only one good T0 trust it.
  else if (GoodT0 == 1) {
    TTPCT0 ScintT0;
    ScintT0 = GoodPath->GetTTPCT0();
    // Use same value of ~11mm or 150ns used above.
    if ( fabs(ScintT0.GetT0() - CathT0.GetT0()) < 150.)
      T0found = ScintT0;
    else
      T0found = CathT0;
  }
  // If there are two, check if they are compatible.
  else {
    // If the T0s are within 40 ns, i.e. ~3mm, then pick the one from the longest segment just in case.
    if ( fabs(Paths[0]->GetT0() - Paths[1]->GetT0()) < 40){
      unsigned int NbSeedClu[2];
      for (unsigned int p = 0; p < 2; p++){
        NbSeedClu[p] = 0;
        ND::THandle<ND::THitSelection> HVclu = Paths[p]->GetHits();
        for (ND::THitSelection::const_iterator Hit = HVclu->begin(); Hit != HVclu->end(); Hit++) {
          ND::THandle<ND::TTPCHVCluster> HV = (*Hit);
          if( HV->isOkForSeed() ) NbSeedClu[p]++;
        }
      }
      if (NbSeedClu[0] > NbSeedClu[1])
        T0found = Paths[0]->GetTTPCT0();
      else
        T0found = Paths[1]->GetTTPCT0();
    }
    // If the T0s are incompatible, pick the one closest to the cathode T0
    else {
      if ( fabs(Paths[0]->GetT0() - CathT0.GetT0()) < fabs(Paths[1]->GetT0() - CathT0.GetT0()))
        T0found = Paths[0]->GetTTPCT0();
      else
        T0found = Paths[1]->GetTTPCT0();
    }
  }
}


//*****************************************************************************
void ND::TTPCT0Finder::FillT0Histos(int bin, int tpc, TTPCT0 &T0Res){
  if (!fEnableHistos) return;
  hUsedt0From->Fill(double(bin));
  hChi2t0->Fill( TMath::Log( T0Res.GetChi2() ) );
  hChi2t0nolog->Fill( T0Res.GetChi2()  );
  double yyfgd = T0Res.GetHit()->GetPosition().Y();
  double zzfgd = T0Res.GetHit()->GetPosition().Z();
  if (tpc==1) hPositiontpc1t0->Fill( zzfgd, yyfgd );
  if (tpc==2) hPositiontpc2t0->Fill( zzfgd, yyfgd );
  if (tpc==3) hPositiontpc3t0->Fill( zzfgd, yyfgd );
}



//*****************************************************************************
void ND::TTPCT0Finder::matchStateToHit(const ND::THitSelection& hits, const State &seed, TTPCT0 &T0Res){
  double minchi2 = fChi2_Max;
  double maxcharge = 0.0;
  ND::THandle<ND::THit> minHit;
  HyperVector resHV;
  HyperVector minResHV;

  // filter each of the hits. Only forward fitting is performed
  for (ND::THitSelection::const_iterator tt = hits.begin(); tt != hits.end(); tt++) {
    ND::THandle<ND::THit> hit = *tt;
    double chi2ndf = 1000000.;
    // creates a RecPack measurement
    Measurement meas;
    bool ok1 = ND::converter().THit_to_Measurement(hit, meas);
    if (!ok1) continue;

    // perform the matching 
    bool ok = ND::rpman().matching_svc().match_state_measurement(seed, meas, chi2ndf, resHV);
    if ( ok ){
      bool chi2OK = (chi2ndf < minchi2) || (TMath::Abs((chi2ndf-minchi2)/minchi2) < 1.e-5);
      bool residualsOK = true;
      for (int k=0; k < resHV.vector().num_row(); k++)
        if (fabs(resHV.vector()[k]) > fSafetyResidualCut ) residualsOK = false;
      if ( ( chi2OK && residualsOK && (hit->GetCharge() > maxcharge)) && T0Res.HitTimeWithinT0Range(hit) ){
        minchi2 = chi2ndf;
        maxcharge = hit->GetCharge();
        minHit = hit;
        minResHV = resHV;
      }
    }
  }
  if (minchi2 < fChi2_Max){
    T0Res.LoadMatchedHit(minHit, minchi2, minResHV);
  }

}


//*****************************************************************************
void ND::TTPCT0Finder::fitT0fromFGD(ND::THitSelection fgd, State& seedState, TTPCT0 &thisT0) {

  if (!enabledFGD) return;

  if (!fgd.size()) return;


  matchStateToHit(fgd, seedState, thisT0);

  // recpack track matching failed to find a hit
  // pick the best hit to use for the t0.
  if ((!thisT0.MatchFound()) && fgd.size() > 0) {
    // Make sure we have the right representation
    RP::rep().convert(seedState, RP::pos_dir_curv);
    EVector vect = seedState.vector();
    double tay = vect[1];
    double taz = vect[2];
    double tdydz = 100.;
    if (fabs(vect[5]) != 0.) tdydz = vect[4] / vect[5];

    if (ND::tpcDebug().T0Finder(DB_VERBOSE)) {
      std::cout<<" RecPack matching failed, try closest hit search"<<std::endl;
    }
    double minchi2 = fChi2_Max;
    ND::THandle<ND::THit> minHit;
    for (ND::THitSelection::const_iterator tt = fgd.begin(); tt != fgd.end(); tt++) {
      const ND::THandle<ND::THit> hit = *tt;
      if (hit->IsYHit()) {
        double hay = hit->GetPosition().Y();
        double haz = hit->GetPosition().Z();
        double adist = fabs( hay - tdydz*haz - (tay - tdydz*taz) ) / sqrt( 1.0 + tdydz*tdydz );
        double chi2 = adist*adist / 100.0;

        if (chi2 < minchi2 && thisT0.HitTimeWithinT0Range(hit) && adist < fSafetyResidualCut) {
          minchi2 = chi2;
          minHit = hit;
        }
      }
    }
    if (minchi2 < fChi2_Max){
      thisT0.LoadClosestHit(minHit, minchi2);
    }

  }

  if (ND::tpcDebug().T0Finder(DB_VERBOSE)) {
    std::cout << " FGD Best Matching " << std::endl;
    std::cout << " FGD Min Chi2  " << thisT0.GetChi2() << " time  " << thisT0.GetT0() << ", max charge " << thisT0.GetHitCharge() << std::endl;
  }



  if ( (!thisT0.OneHitFound()) && ND::tpcDebug().T0Finder(DB_VERBOSE)) {
    // Make sure we have the right representation
    RP::rep().convert(seedState, RP::pos_dir_curv);
    EVector vect = seedState.vector();
    std::cout << " Pos = ( " << vect[0] << "," << vect[1] << "," << vect[2] << ")" << std::endl;
    std::cout << " Dir = ( " << vect[3] << "," << vect[4] << "," << vect[5] << ")" << std::endl;
    std::cout << " rho = ( " << vect[6] << " ) " << std::endl;
  }

  if (fEnableHistos) hFGDt0->Fill(thisT0.GetT0());

}



//*****************************************************************************
void ND::TTPCT0Finder::fitT0fromP0D(ND::THitSelection fP0D, const State& seedState, TTPCT0 &thisT0) {

  if (!enabledP0D) return;

  if (!fP0D.size()) return;

  matchStateToHit(fP0D, seedState, thisT0);

  if( ND::tpcDebug().T0Finder(DB_VERBOSE) ){
    std::cout << " P0D Best Matching " << std::endl;
    std::cout << " P0D Min Chi2  " << thisT0.GetChi2() << " time  " << thisT0.GetT0() << ", max charge " << thisT0.GetHitCharge() << std::endl;
  }


  if (fEnableHistos==true) hP0Dt0->Fill(thisT0.GetT0());

}



//*****************************************************************************
void ND::TTPCT0Finder::fitT0fromSMRD(ND::THitSelection fSMRD, const State& seedState, TTPCT0 &thisT0) {

  if( !enabledSMRD ) return;

  if( !fSMRD.size() ) return;

  matchStateToHit(fSMRD, seedState, thisT0);

  if (ND::tpcDebug().T0Finder(DB_VERBOSE)) {
    std::cout << " SMRD Best Matching " << std::endl;
    std::cout << " SMRD Min Chi2  " << thisT0.GetChi2() << " time  " << thisT0.GetT0() << ", max charge " << thisT0.GetHitCharge() << std::endl;
  }


  if (fEnableHistos==true) hSMRDt0->Fill(thisT0.GetT0());

}



//*****************************************************************************
void ND::TTPCT0Finder::fitT0fromECAL(ND::THitSelection ecal, State& seedState, TTPCT0 &thisT0) {

  if (!enabledECAL) return;

  if (!ecal.size()) return;


  matchStateToHit(ecal, seedState, thisT0);


  // recpack track matching failed to find a hit
  // pick the best hit to use for the t0.
  double maxq = 0.0;
  if ( (!thisT0.MatchFound()) && ecal.size() > 0) {
    // Make sure we have the right representation
    RP::rep().convert(seedState, RP::pos_dir_curv);
    EVector vect = seedState.vector();
    double tay = vect[1];
    double taz = vect[2];
    double tdydz = 100.;
    if (fabs(vect[5]) > 0.) tdydz = vect[4] / vect[5];

    if (ND::tpcDebug().T0Finder(DB_VERBOSE)) {
      std::cout<<" RecPack matching failed, try closest hit search"<<std::endl;
    }

    double minchi2 = fChi2_Max;
    ND::THandle<ND::THit> minHit;
    for (ND::THitSelection::const_iterator tt = ecal.begin(); tt != ecal.end(); tt++) {
      ND::THandle<ND::THit> hit = *tt;
      if (hit->IsYHit() && hit->IsZHit()) {
        double hay = hit->GetPosition().Y();
        double haz = hit->GetPosition().Z();
        double adist = fabs( hay - tdydz*haz - (tay - tdydz*taz) ) / sqrt( 1.0 + tdydz*tdydz );
        double chi2 = adist*adist / 100.0;

        if ( (chi2 < minchi2 || ( chi2 == minchi2 && hit->GetCharge()>maxq )) && thisT0.HitTimeWithinT0Range(hit) && adist < fSafetyResidualCut) {
          maxq = hit->GetCharge();
          minchi2 = chi2;
          minHit = hit;
        }
      }
    }
    if (minchi2 < fChi2_Max){
      thisT0.LoadClosestHit(minHit, minchi2);
    }

  }


  if (ND::tpcDebug().T0Finder(DB_VERBOSE)) {
    std::cout << " ECAL Best Matching " << std::endl;
    std::cout << " ECAL Min Chi2  " << thisT0.GetChi2() << " time  " << thisT0.GetT0() << ", max charge " << thisT0.GetHitCharge() << std::endl;
  }

  if (fEnableHistos==true) hECALt0->Fill(thisT0.GetT0());

}

//*****************************************************************************
void ND::TTPCT0Finder::prepareFGDHits(ND::TND280Event& event ){

  if( ND::tpcDebug().T0Finder(DB_VERBOSE)   )
    std::cout << " Preparing FGD hits " << std::endl;
  if( !enabledFGD ) return;

  ND::THandle<ND::THitSelection> fgd_raw = event.GetHitSelection("fgd");

  if( !fgd_raw ) return;

  fFGD_TPC1.SetName("fFGD_TPC1");
  fFGD_TPC2.SetName("fFGD_TPC2");
  fFGD_TPC3.SetName("fFGD_TPC3");
  double min_charge = fFGD_Qmin;

  if( fgd_raw ) {

    if( ND::tpcDebug().T0Finder(DB_VERBOSE) )
      std::cout << " FGD RAW hits = " << fgd_raw->size() << std::endl;

    double fgdymax = ND::TGeomInfo::Get().FGD().FGD1ActiveMax()[1]-fFGD_YedgeSize;
    double fgdymin = ND::TGeomInfo::Get().FGD().FGD1ActiveMin()[1]+fFGD_YedgeSize;
    int nplanes = ND::TGeomInfo::Get().FGD().ActivePlaneCount();
    int lastFGD1plane = ND::TGeomInfo::Get().FGD().GetLastFGD1Plane();
    int firstFGD2plane = ND::TGeomInfo::Get().FGD().GetFirstFGD2Plane();

    for (ND::THitSelection::const_iterator ht = fgd_raw->begin(); ht != fgd_raw->end();++ht){
      bool isinfgd1 = ND::TGeomInfo::Get().FGD().IsInFGD1( (*ht)->GetPosition() );
      bool isinfgd2 = ND::TGeomInfo::Get().FGD().IsInFGD2( (*ht)->GetPosition() );
      int plane     = ND::TGeomInfo::Get().FGD().ActivePlane( (*ht)->GetPosition().Z() );

      if( !((*ht)->IsYHit()) ) continue;

      if( (*ht)->GetCharge() > min_charge ) {

        // relevant for TPC1: first few FGD1 planes, or FGD1 near edge
        if ( plane < fFGD_Layers ) fFGD_TPC1.push_back(*ht);
        if ( isinfgd1 == true ){
          if ( (*ht)->GetPosition().Y() > fgdymax )  fFGD_TPC1.push_back(*ht);
          if ( (*ht)->GetPosition().Y() < fgdymin )  fFGD_TPC1.push_back(*ht);
        }

        // relevant for TPC2: last few FGD1 planes, first few FGD2 planes, or y near edge
        if ( plane > lastFGD1plane-fFGD_Layers && plane < firstFGD2plane+fFGD_Layers ) fFGD_TPC2.push_back(*ht);
        if ( (*ht)->GetPosition().Y() > fgdymax )  fFGD_TPC2.push_back(*ht);
        if ( (*ht)->GetPosition().Y() < fgdymin )  fFGD_TPC2.push_back(*ht);

        // relevant for TPC3: last few planes of FGD2, or FGD2 near edge
        if ( plane >= nplanes-fFGD_Layers ) fFGD_TPC3.push_back(*ht);
        if ( isinfgd2 == true ){
          if ( (*ht)->GetPosition().Y() > fgdymax )  fFGD_TPC3.push_back(*ht);
          if ( (*ht)->GetPosition().Y() < fgdymin )  fFGD_TPC3.push_back(*ht);
        }
      }
    }
  }

  if( ND::tpcDebug().T0Finder(DB_VERBOSE) )
    std::cout << " Prepared FGD  hits for TPC1 " << fFGD_TPC1.size() << " for TPC2 "<<fFGD_TPC2.size()<<" for TPC3 "<<fFGD_TPC3.size()<<std::endl;fflush(stdout);
}


//*****************************************************************************
void ND::TTPCT0Finder::prepareECALHits(ND::TND280Event& event){

  if ( !enabledECAL ) return;

  if( ND::tpcDebug().T0Finder(DB_VERBOSE)   )
    std::cout << " Preparing FGD hits " << std::endl;

  fECAL_TPC1.SetName("fECAL_TPC1");
  fECAL_TPC2.SetName("fECAL_TPC2");
  fECAL_TPC3.SetName("fECAL_TPC3");
  fECAL_SIDE.SetName("fECAL_SIDE");

  ND::THandle<ND::THitSelection>  hits;
  std::list< ND::THandle<ND::THitSelection> > hitslist;

  double min_charge = fECAL_Qmin;

  ND::THandle< ND::THit > lasthit;
  bool lastSideECalHit = false;

  hits = event.GetHitSelection("barrel_RT");
  if( hits ) {
    hitslist.push_back( hits );
  }
  hits = event.GetHitSelection("barrel_LT");
  if( hits ) {
    hitslist.push_back( hits );
  }
  hits = event.GetHitSelection("barrel_RB");
  if( hits ) {
    hitslist.push_back( hits );
  }
  hits = event.GetHitSelection("barrel_LB");
  if( hits ) {
    hitslist.push_back( hits );
  }
  hits = event.GetHitSelection("barrel_RS");
  if( hits ) {
    hitslist.push_back( hits );
  }
  hits = event.GetHitSelection("barrel_LS");
  if( hits ) {
    hitslist.push_back( hits );
  }
  hits = event.GetHitSelection("dsecal");
  if( hits ) {
    hitslist.push_back( hits );
  }


  for (std::list< ND::THandle<ND::THitSelection> >::iterator hh = hitslist.begin(); hh != hitslist.end(); ++hh) {
    hits = *hh;
    if (hits){
      for (ND::THitSelection::const_iterator ht = hits->begin(); ht != hits->end();++ht){
        ND::THandle< ND::THit > hit = (*ht);
        // below is what "layer" was in head version, but that wasnt working for me...
        // may need to change this if geominfo changes?
        // int layer = ND::TGeomInfo::Get().ECAL().GetLayerNumber((*ht)->GetGeomId());
        ND::GeomId::Def::DetectorId det;
        ND::GeomId::Def::ECal::MagnetClams clam;
        ND::GeomId::Def::ECal::ECalModules mod;
        ND::TGeomInfo::ECAL().GetModule(hit).GetIdentity(det, clam, mod);
        int layer = ND::TGeomInfo::Get().ECAL().GetLayerNumber( hit->GetPosition());
        // Complicated but necessary to get hits that don't depend on the X position of the track
        // - For all module but the side modules, use YZ hits
        // - For the side modules, use XZ hits
        bool GoodOrientation;
        bool SideECalHit = false;
        if ((hit->GetGeomId().GetSubsystemId() == ND::GeomId::Def::kDSECal) || mod == ND::GeomId::Def::ECal::kTopModule || mod == ND::GeomId::Def::ECal::kBottomModule){
          GoodOrientation = hit->IsYHit() && hit->IsZHit();
        } else {
          GoodOrientation = hit->IsXHit() && hit->IsYHit();
          SideECalHit = true;
        }
        if ( GoodOrientation &&
            hit->GetCharge() > min_charge &&
            layer < fECAL_Layers &&
            fabs( hit->GetTime() ) < 50000.0 ){ // add a time range (avoid weird timed hits)
          bool DiffFromLast = true;
          if (lasthit){
            if (SideECalHit != lastSideECalHit){
              DiffFromLast = true;
            } else {
              if ( !SideECalHit)
                DiffFromLast = (hit->GetPosition().Y() != lasthit->GetPosition().Y() && hit->GetPosition().Z() != lasthit->GetPosition().Z());
              else
                DiffFromLast = (hit->GetPosition().Y() != lasthit->GetPosition().Y() && hit->GetPosition().X() != lasthit->GetPosition().X());
            }
          }
          if (DiffFromLast){
            if (!SideECalHit){
              // tpc1 allow ecal hits from roughly -1300 to 600mm
              if (hit->GetPosition().Z() > -1300.0 && hit->GetPosition().Z() < 600.0) {
                fECAL_TPC1.push_back(hit);
              }
              // tpc2 allow ecal hits from roughly 0 to 2000mm
              if ((hit)->GetPosition().Z() > 0.0 && (hit)->GetPosition().Z() < 2000.0) {
                fECAL_TPC2.push_back(hit);
              }
              // tpc3 allow ecal hits from roughly 1100mm onward
              if ((hit)->GetPosition().Z() > 1100.0) {
                fECAL_TPC3.push_back(hit);
              }
            } else {
              fECAL_SIDE.push_back(hit);
            }
            
            lasthit = hit;
            lastSideECalHit = SideECalHit;
          }
        }
      }
    }
  }

  if( ND::tpcDebug().T0Finder(DB_VERBOSE) )
    std::cout << " Prepared ECAL hits for TPC1 " << fECAL_TPC1.size() << " for TPC2 "<< fECAL_TPC2.size()<<" for TPC3 "<< fECAL_TPC3.size()<<" from side "<<fECAL_SIDE.size()<< std::endl;fflush(stdout);
}



//*****************************************************************************
void ND::TTPCT0Finder::prepareP0DHits(ND::TND280Event& event){

  if( !enabledP0D ) return;

  ND::THandle<ND::THitSelection> hits = event.GetHitSelection("p0d");

  if( !hits ) return;

  if( ND::tpcDebug().T0Finder(DB_VERBOSE)   )
    std::cout << " Preparing P0D hits " << std::endl;

  fP0D.SetName("p0d");

  double min_charge = fP0D_Qmin;
  if( hits ) {
    int nplanes = ND::TGeomInfo::Get().P0D().ActivePlaneCount();

    for (ND::THitSelection::const_iterator ht = hits->begin(); ht != hits->end();++ht){
      int plane = ND::TGeomInfo::Get().P0D().ActivePlane( (*ht)->GetPosition()[2] );
      if ( (*ht)->IsYHit() && (*ht)->GetCharge() > min_charge && nplanes-plane < fP0D_Layers ){
        fP0D.push_back(*ht);
      }
    }
  }

  if( ND::tpcDebug().T0Finder(DB_VERBOSE) )
    std::cout << " Prepared P0D  hits " << fP0D.size() << std::endl;fflush(stdout);

  return;
}

//*****************************************************************************
void ND::TTPCT0Finder::prepareSMRDHits(ND::TND280Event& event){

  if( !enabledSMRD ) return;

  ND::THandle<ND::THitSelection> hits = event.GetHitSelection("fSMRD");

  if( !hits ) hits = event.GetHitSelection("mrd");
  if( !hits ) hits = event.GetHitSelection("fSMRD");
  if( !hits ) hits = event.GetHitSelection("smrdReconHits");

  //  std::cout << " P0D hits " << fP0D->size() << std::endl;fflush(stdout);

  fSMRD.SetName("t0smrd");

  double min_charge = fSMRD_Qmin;
  if( hits ) {
    for (ND::THitSelection::const_iterator ht = hits->begin(); ht != hits->end();++ht){
      if ( (*ht)->IsYHit() && (*ht)->GetCharge() > min_charge ){
        fSMRD.push_back(*ht);
      }
    }
  }

  if( ND::tpcDebug().T0Finder(DB_VERBOSE) )
    std::cout << " Prepared SMRD hits " << fSMRD.size() << std::endl;fflush(stdout);

  return;
}


//*****************************************************************************
void ND::TTPCT0Finder::CleanHits(){
  if( fFGD_TPC1.size() > 0 ) fFGD_TPC1.erase(fFGD_TPC1.begin(),fFGD_TPC1.end());
  if( fFGD_TPC2.size() > 0 ) fFGD_TPC2.erase(fFGD_TPC2.begin(),fFGD_TPC2.end());
  if( fFGD_TPC3.size() > 0 ) fFGD_TPC3.erase(fFGD_TPC3.begin(),fFGD_TPC3.end());
  if( fECAL_TPC1.size() > 0 ) fECAL_TPC1.erase(fECAL_TPC1.begin(),fECAL_TPC1.end());
  if( fECAL_TPC2.size() > 0 ) fECAL_TPC2.erase(fECAL_TPC2.begin(),fECAL_TPC2.end());
  if( fECAL_TPC3.size() > 0 ) fECAL_TPC3.erase(fECAL_TPC3.begin(),fECAL_TPC3.end());
  if( fECAL_SIDE.size() > 0 ) fECAL_SIDE.erase(fECAL_SIDE.begin(),fECAL_SIDE.end());
  if( fP0D.size() > 0 ) fP0D.erase(fP0D.begin(),fP0D.end());
  if( fSMRD.size() > 0 ) fSMRD.erase(fSMRD.begin(),fSMRD.end());
}



