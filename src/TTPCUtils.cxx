#include "TTPCUtils.hxx"
#include "TTPCHVCluster.hxx"
#include "TTPCHitPad.hxx"
#include "TTPCCalibration.hxx"

#include <TG4TrajectoryPoint.hxx>
#include <TG4Trajectory.hxx>
#include <TrackTruthInfo.hxx>
#include <TND280Event.hxx>
#include <TEventFolder.hxx>
#include <TRecPackManager.hxx>
#include <TComboHit.hxx>
#include <TRealDatum.hxx>
#include <TIntegerDatum.hxx>
#include <TND280Event.hxx>
#include <TrackingUtils.hxx>
#include <ReconObjectUtils.hxx>

#include <sstream>

namespace TTPCUtils {

  //*****************************************************************************
  ND::THandle<ND::TReconVertex> CopyCreateTReconVertex(ND::TReconVertex *Original){
    ND::THandle<ND::TReconVertex> newVertex( new ND::TReconVertex() );
    
    newVertex->SetQuality(Original->GetQuality());
    newVertex->SetStatus(Original->GetStatus());
    newVertex->SetNDOF(Original->GetNDOF());
    newVertex->AddDetector(Original->GetDetectors());
    
    // The constituents in TTPCJunction are TTPCPaths that are converted into TReconTrack by TTPCPattern.
    // So TTPCPattern will take care of copying the right TReconTracks as constituents of this TReconVertex.
    // Bottom line: Don't copy the constituents here !

    // Explicitly look for TRealDatum and TIntegerDatum objects and copy them.
    for (ND::TReconBase::iterator obj = Original->begin(); 
         obj != Original->end(); ++obj) {
        ND::TDatum* datum = dynamic_cast<ND::TDatum*>(*obj);
        // Copy any TRealDatum objects.
        ND::TRealDatum *rDatum = dynamic_cast<ND::TRealDatum*>(datum);
        if (rDatum) {
            ND::TRealDatum* nDatum = new ND::TRealDatum(*rDatum);
            newVertex->AddDatum(nDatum);
            continue;
        }
        // Copy any TIntegerDatum objects.
        ND::TIntegerDatum *iDatum = dynamic_cast<ND::TIntegerDatum*>(datum);
        if (iDatum) {
            ND::TIntegerDatum* nDatum = new ND::TIntegerDatum(*iDatum);
            newVertex->AddDatum(nDatum);
            continue;
        }
    }

    ND::THandle<ND::TVertexState> oldState = Original->GetState();
    ND::THandle<ND::TVertexState> newState = newVertex->GetState();
    *newState = *oldState;

    return newVertex;
  }



  //*****************************************************************************
  double State2CluChi2(State helixState, ND::THandle<ND::TTPCHVCluster> Cluster){
    double chi2[3];
    int ndof[3];           

    return State2CluChi2(helixState, Cluster, chi2, ndof);
  }

  //*****************************************************************************
  double State2CluChi2(State helixState, ND::THandle<ND::TTPCHVCluster> Cluster, double Chi2[3], int NDOF[3]){
    double totalChi2 = 0.0;

    // Chi2 is calculated here 
   
    EVector helixStateVec(6,0);
    EMatrix helixStateCov(6,6,0);

    helixStateVec = helixState.vector();
    helixStateCov = helixState.matrix();    
   
    double xCluster = Cluster->CalibX();
    double yCluster = Cluster->CalibY();
    double zCluster = Cluster->CalibZ();

    for ( int i = 0; i < 3; i++){
      Chi2[i] = 0.0;
      NDOF[i] = 0;
    }

    // fudge factors take into account diffusion effects and is cut off at min. value 0

    double sigmaCX = 0., sigmaCY = 0., sigmaCZ = 0.;

    if(Cluster->IsVertical()){
      sigmaCX = 0.27 - TMath::Abs(helixStateVec[0]/1000.)*0.4 + TMath::Power(helixStateVec[0]/1000.,2)*0.09;
      sigmaCX = (TMath::Abs(sigmaCX)+sigmaCX)/2.;

      sigmaCY = 0.65 - TMath::Abs(helixStateVec[0]/1000.)*1.45 + TMath::Power(helixStateVec[0]/1000.,2)*2.85;
      sigmaCY = (TMath::Abs(sigmaCY)+sigmaCY)/2.;

      Chi2[0] = std::pow((helixStateVec[0]-xCluster),2)/(helixStateCov[0][0]+sigmaCX);
      NDOF[0] = 1;

      Chi2[1] = std::pow((helixStateVec[1]-yCluster),2)/(helixStateCov[1][1]+sigmaCY);
      NDOF[1] = 1;

    } else {
      sigmaCX = 0.54 - TMath::Abs(helixStateVec[0]/1000.)*0.69 + TMath::Power(helixStateVec[0]/1000.,2)*0.15;
      sigmaCX = (TMath::Abs(sigmaCX)+sigmaCX)/2.;

      sigmaCZ = 1.62 - TMath::Abs(helixStateVec[0]/1000.)*0.84 + TMath::Power(helixStateVec[0]/1000.,2)*5.68;
      sigmaCZ = (TMath::Abs(sigmaCZ)+sigmaCZ)/2.;

      Chi2[0] = std::pow((helixStateVec[0]-xCluster),2)/(helixStateCov[0][0]+sigmaCX);
      NDOF[0] = 1;

      Chi2[2] = std::pow((helixStateVec[2]-zCluster),2)/(helixStateCov[2][2]+sigmaCZ);
      NDOF[2] = 1;
    }

    totalChi2 = Chi2[0] + Chi2[1] + Chi2[2];

    return totalChi2;
  }

  //*****************************************************************************
  bool SafeSort( double first, double second ){
    // Use int because I don't trust the ordering of doubles
    // since the bug 871 was discovered.
    // Millimeter precision is enough.
    return ( int(first) < int(second) );
  }


  //*****************************************************************************
  int SenseFromTwoClusters(ND::THandle<ND::TTPCPath> Path2, ND::THandle<ND::TTPCHVCluster> ChosenClu){
  //*****************************************************************************
    int sense = 0;
    if (Path2->GetHits()->size() < 2)
      return 0;
    ND::THandle<ND::TTPCHVCluster> tmpClu;
    ND::THandle<ND::TTPCHVCluster> tmprClu;
    ND::THandle<ND::TTPCHVCluster> nextClu;
    ND::THitSelection::const_iterator Clu = Path2->GetHits()->begin();
    ND::THitSelection::const_reverse_iterator rClu = Path2->GetHits()->rbegin();
    tmpClu = *Clu;
    tmprClu = *rClu;
    if (ChosenClu == tmpClu){
      Clu++;
      nextClu = *Clu;
    } else if (ChosenClu == tmprClu){
      rClu++;
      nextClu = *rClu;
    } else {
      // PROBLEM. That will probably be a bigger problem outside of this function regardless of the sense.
      return 0;
    }

    // Find the crude track "sense" in z(y) when the first two clusters are vertical(horizontal).
    // If the first two clusters have different orientation, just use 0 and deal with it later.
    if ( ChosenClu->IsVertical() == nextClu->IsVertical()){
      if (ChosenClu->IsVertical()){
        sense = nextClu->Z() - ChosenClu->Z();
      }else{
        sense = nextClu->Y() - ChosenClu->Y();
      }
    }
    return sense;
  }


  //*****************************************************************************
  void FindClosestEnds(ND::THandle<ND::TReconBase> PathA, ND::THandle<ND::TReconBase> PathB, unsigned int &EndA, unsigned int &EndB){
    // TODO: Add check on the presence of clusters and use exception if not there.
    ND::THandle<ND::THit> FirstCluA = *(PathA->GetHits()->begin());
    ND::THandle<ND::THit> LastCluA = *(PathA->GetHits()->rbegin());
    ND::THandle<ND::THit> FirstCluB = *(PathB->GetHits()->begin());
    ND::THandle<ND::THit> LastCluB = *(PathB->GetHits()->rbegin());
    std::vector<double> Sorter;
    double FtFt = (FirstCluA->GetPosition() - FirstCluB->GetPosition()).Mag();
    double FtLt = (FirstCluA->GetPosition() - LastCluB->GetPosition()).Mag();
    double LtLt = (LastCluA->GetPosition() - LastCluB->GetPosition()).Mag();
    double LtFt = (LastCluA->GetPosition() - FirstCluB->GetPosition()).Mag();
    Sorter.push_back( FtFt );
    Sorter.push_back( FtLt );
    Sorter.push_back( LtLt );
    Sorter.push_back( LtFt );
    std::stable_sort(Sorter.begin(), Sorter.end(), TTPCUtils::SafeSort);
    double MinDist = *(Sorter.begin());
    if (MinDist == FtFt || MinDist == FtLt){
      EndA = 0;
    } else {
      EndA = 1;
    }
    if (MinDist == FtFt || MinDist == LtFt){
      EndB = 0;
    } else {
      EndB = 1;
    }
  }

  //*****************************************************************************
  ND::THandle<ND::TTPCPath> MergePaths(ND::THandle<ND::TTPCPath> PathA, ND::THandle<ND::TTPCPath> PathB){
    ND::THandle<ND::TTPCPath> newPath( new ND::TTPCPath() );

    // Create a new path putting in the clusters from the two matched paths.
    // BE CAREFUL with the ordering of the clusters !
    unsigned int UseEndA, UseEndB;
    FindClosestEnds(PathA, PathB, UseEndA, UseEndB);

    ND::THitSelection* newClusters;
    if (UseEndA == 1){
      ND::THandle<ND::THitSelection> clusters = PathA->GetHits();
      newClusters = new ND::THitSelection(*clusters);
      if (UseEndB == 0 ){
        for (ND::THitSelection::const_iterator hit = PathB->GetHits()->begin(); hit != PathB->GetHits()->end(); ++hit)
          newClusters->push_back(*hit);
      } else {
        for (ND::THitSelection::const_reverse_iterator hit = PathB->GetHits()->rbegin(); hit != PathB->GetHits()->rend(); ++hit)
          newClusters->push_back(*hit);
      }
    } else {
      if (UseEndB == 1 ){
        ND::THandle<ND::THitSelection> clusters = PathB->GetHits();
        newClusters = new ND::THitSelection(*clusters);
        for (ND::THitSelection::const_iterator hit = PathA->GetHits()->begin(); hit != PathA->GetHits()->end(); ++hit)
          newClusters->push_back(*hit);
      } else {
        newClusters = new ND::THitSelection();
        for (ND::THitSelection::const_reverse_iterator hit = PathB->GetHits()->rbegin(); hit != PathB->GetHits()->rend(); ++hit)
          newClusters->push_back(*hit);
        for (ND::THitSelection::const_iterator hit = PathA->GetHits()->begin(); hit != PathA->GetHits()->end(); ++hit)
          newClusters->push_back(*hit);
      }
    }
    newPath->AddHits(newClusters);
    newPath->SetId(ND::tpcCalibration().GetPathId());

    return newPath;
  }



  //*****************************************************************************
  ND::THandle<ND::TReconBase> MergeAndFitObjectsWithRecPack(ND::THandle<ND::TReconBase> T1, ND::THandle<ND::TReconBase> T2){
    ND::THandle<ND::TReconBase> NewTrk;
    if(!T1) return NewTrk;
    if(!T2) return NewTrk;
     
    //  if( Param->debuglevel >= INFO ){
    //    ReconPrintout::DumpReconHistory(T1);
    //    ReconPrintout::DumpReconHistory(T2);
    //  }


    ND::THandle< ND::TReconTrack > Track1 = *(T1->Get<ND::TReconObjectContainer>("constituents")->begin());
    ND::THandle< ND::TReconTrack > Track2 = *(T2->Get<ND::TReconObjectContainer>("constituents")->begin());

    bool NeedToReorder = false;
    ND::THandle<ND::TTrackState> State1 = Track1->GetState();
    ND::THandle<ND::TTrackState> State2 = Track2->GetState();
    // Find z ordering of tracks
    double z1 = State1->GetPosition().Z();
    double z2 = State2->GetPosition().Z();
    NeedToReorder = z1 > z2;

    ND::THandle<ND::TReconBase> TA;
    ND::THandle<ND::TReconBase> TB;

    if ( !NeedToReorder){
      ND::tman().ReverseObjectSenseAndCharge(*T1);
      ND::tman().ReverseObjectSenseAndCharge(*T2);
      TA = T1;
      TB = T2;
    } else {
      TA = T2;
      TB = T1;
    }

    ND::THandle<ND::TReconBase>  best;
    ND::THandle<ND::TReconBase>  worst;
    ND::tman().GetBestObject( TA, TB, best, worst);

    //so two objects were successfully matched at this point
    //then merge and refit
    
    bool ok_tpciso1 = true;
    bool ok_tpciso2 = true; 
    
    //check the covs when dealing with TPC iso objects
    ok_tpciso1 = ND::tman().CheckTpcIsoCovariance(best);
    ok_tpciso2 = ND::tman().CheckTpcIsoCovariance(worst); 

    if(!ok_tpciso1)
      return ND::THandle<ND::TReconBase>();

    //copy the main worst state for a while
    ND::THandle<ND::TReconState> worstStateClone(new ND::TPIDState()); 
    ReconObjectUtils::CopyState(*(worst->GetState()), *worstStateClone); 

    if(!ok_tpciso2)
      ND::tman().SetDefaultCovsTPCIso(worst); 


    NewTrk = ND::tman().MergeObjects("PID", TA, TB, true, false, false);
    
    if(!ok_tpciso2)
      *(worst->GetState()) = *worstStateClone;
    
    if(!NewTrk) 
      return NewTrk;

    
    //put as first constituent the best object, just in case
    // NewTrk->AddConstituent(best);
    // NewTrk->AddConstituent(worst);
    NewTrk->AddConstituent(TA);
    NewTrk->AddConstituent(TB);
    
    //store the length sign from RECPACK since may need to change it for curving back tracks
    int length_sign = ND::rpman().model_svc().model(ND::rpman().model_svc().model_name()).intersector().length_sign();
      
    //whether to allow zero length propagation
    bool allow_zero = ND::rpman().model_svc().model(ND::rpman().model_svc().model_name()).intersector().allow_zero_length();
          
    
    ND::rpman().model_svc().model().intersector().set_length_sign(0);
    
    ND::rpman().navigation_svc().navigator().set_force_allow_zero_length(true);
    
    //get two states
    ND::THandle<ND::TReconState> tstate1 = TrackingUtils::GetFirstState(*best); 
   

    //we already know which is the best object, use its first state to define the seed (tstate1)
    /// TODO: Use static parameter somehow
    double RecPackRefittingCovFactor = 100.;
    bool ok = ND::tman().FitObject(*(tstate1), RecPackRefittingCovFactor, NewTrk);
    
    // return back
    ND::rpman().model_svc().model(ND::rpman().model_svc().model_name()).intersector().set_length_sign(length_sign);
      
    ND::rpman().navigation_svc().navigator().set_force_allow_zero_length(allow_zero);
     
    if(!ok)
      return ND::THandle<ND::TReconBase>();
    
    //use info from the best object
    if (best->Get< ND::TIntegerDatum >("T0Source")){
      double t0source = TA->Get< ND::TIntegerDatum >("T0Source")->GetValue();
      ND::TIntegerDatum* t0src = new ND::TIntegerDatum("T0Source",t0source);    
      NewTrk->AddDatum(t0src);
    }

    ND::THandle<ND::TReconPID> Tlongest = best;
    if(!Tlongest)
      return ND::THandle<ND::TReconBase>();
    
    // Copy the PID info from the longest (best) segment to the new track
    TString PIDdata [19] = {"tpcPid_Craw", "tpcPid_Ccorr", "tpcPid_NTrun", "tpcPid_SampleLength", "tpcPid_dEdxexpEle", "tpcPid_SigmaEle", "tpcPid_PullEle", "tpcPid_dEdxexpMuon", "tpcPid_SigmaMuon", "tpcPid_PullMuon", "tpcPid_dEdxexpPion", "tpcPid_SigmaPion", "tpcPid_PullPion", "tpcPid_dEdxexpKaon", "tpcPid_SigmaKaon", "tpcPid_PullKaon", "tpcPid_dEdxexpProton", "tpcPid_SigmaProton", "tpcPid_PullProton"}; 
    for (int p = 0; p < 19; p++){
      if (Tlongest->Get< ND::TRealDatum >(PIDdata[p])){
        double value = Tlongest->Get< ND::TRealDatum >(PIDdata[p])->GetValue();
        ND::TRealDatum* TRD = new ND::TRealDatum(PIDdata[p],value);    
        NewTrk->AddDatum(TRD);
      }
    }
    
    ND::THandle<ND::TReconPID> NewPID = NewTrk;
    NewPID->SetParticleId(Tlongest->GetParticleId());
    NewPID->SetPIDWeight(Tlongest->GetPIDWeight());

    return NewTrk;
  }
}
