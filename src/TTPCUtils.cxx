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
  bool SortShortestToLongest( ND::THandle<ND::TTPCPath> first, ND::THandle<ND::TTPCPath> second ){

    unsigned int length1 = first->GetHits()->size();
    unsigned int length2 = second->GetHits()->size();

    bool Order = false;

    if( length1 < length2){
      Order = true;
    }

    return Order; 
  }

  //*****************************************************************************
  std::string DetectorToString( ND::THandle<ND::TReconBase> track ){
    if (track->UsesDetector(ND::TReconBase::kTPC1))
      return std::string("TPC1");
    else if (track->UsesDetector(ND::TReconBase::kTPC2))
      return std::string("TPC2");
    else if (track->UsesDetector(ND::TReconBase::kTPC3))
      return std::string("TPC3");
    // We should not arrive here !
    std::cout<<"TTPCUtils::DetectorToString: WARNING. The track you provided doesn't have TPC detector bit !"<<std::endl;
    return "UNKNOWN";
  }



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
  double GetXWithT0(double Time, double T0, ND::THandle<ND::THit> RefHit){
    ND::THandle<ND::TTPCHitPad> hitPad = RefHit;
    double RPx = hitPad->GetPosition().X();
    double Sense = hitPad->DriftSense();
    double DriftDistance = (Time - T0 - ND::tpcCalibration().GetTimeOffset()) * ND::tpcCalibration().GetDriftVelocity();
    // How far are we from the read out plane ?
    return ( RPx - (Sense * DriftDistance));
  }


  //*****************************************************************************
  State InvertDirection( State inputState ){
    if(inputState.name(RP::representation) != RP::pos_dir_curv)
      RP::rep().convert(inputState, RP::pos_dir_curv);
    State newState = inputState;
    EVector outVect = inputState.vector();
    EMatrix outCova = inputState.matrix();
    outVect[3] = -1.* outVect[3];
    outVect[4] = -1.* outVect[4];
    outVect[5] = -1.* outVect[5];
    newState.set_hv(HyperVector(outVect,outCova));

    return newState;
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



  //*****************************************************************************
  ND::THandle<ND::TG4Trajectory> FindTrueTrajectories( ND::THandle<ND::THitSelection> recoHits, double &complete, double &clean){
    ND::TND280Event* event = ND::TEventFolder::GetCurrentEvent();
    if (! event->GetContext().IsMC() ) 
      return ND::THandle<ND::TG4Trajectory> ();

    ND::THitSelection tmpHitSel;
    // Convert the clusters to TComboHits if needed cause that's what we need
    // for the truth matching functions
    for (ND::THitSelection::const_iterator tmpClu = recoHits->begin(); tmpClu != recoHits->end(); tmpClu++) {
      ND::THandle<ND::TTPCHVCluster> HVClu = *tmpClu; 
      if (HVClu){
        ND::THandle<ND::TComboHit> oaEClu = HVClu->ConvertToOAEvent();
        tmpHitSel.push_back(oaEClu);
      } else {
        tmpHitSel.push_back(*tmpClu);
      }
    }

    // First get the true track index
    int G4Id = TrackTruthInfo::GetG4TrajIDHits(tmpHitSel, complete, clean);
    ND::THandle<ND::TG4TrajectoryContainer> trajectories = event->Get<ND::TG4TrajectoryContainer>("truth/G4Trajectories");

    // Then get the trajectory
    return (trajectories->GetTrajectory(G4Id));
  }

  //*****************************************************************************
  bool TrueStateNear3Dpoint( ND::THandle<ND::THitSelection> recoHits, TVector3 pos3d, State &trueState){
    double complete = 0.0;
    double clean = 0.0;
    ND::THandle<ND::TG4TrajectoryContainer> trajectories;
    ND::THandle<ND::TG4Trajectory> G4Traj = FindTrueTrajectories(recoHits, complete, clean);
    if (!G4Traj)
      return false;
    
    // Finally find the true point closest to the first cluster
    double dist;
    ND::TG4TrajectoryPoint truPoint = TrackTruthInfo::GetCloserG4Point(G4Traj, &pos3d,dist);

    // Convert the true state into a HelixModel
    EVector vect = EVector(7,0);
    EMatrix cova = EMatrix(7,7,0);
    vect[0] = truPoint.GetPosition().X();
    vect[1] = truPoint.GetPosition().Y();
    vect[2] = truPoint.GetPosition().Z();
    TVector3 Dir = truPoint.GetMomentum()*(1/truPoint.GetMomentum().Mag());
    vect[3] = Dir.X();
    vect[4] = Dir.Y();
    vect[5] = Dir.Z();
    double Charge = 0.0;
    int PDGcode = G4Traj->GetPDGEncoding();
    unsigned int PDGcodeAbs = (unsigned int) PDGcode;
    if (PDGcodeAbs == 11 || PDGcodeAbs == 13 || PDGcodeAbs == 211){
      Charge = -1.0;
    } else if (PDGcodeAbs == 2212 || PDGcodeAbs == 211){
      Charge = 1.0;
    }
    vect[6] = Charge/truPoint.GetMomentum().Mag();
    trueState.set_hv(HyperVector(vect,cova));
    trueState.set_hv(RP::sense, HyperVector(1,0));  
    trueState.set_name(RP::representation,RP::pos_dir_curv);

    return true;
  }


  

  //**************************************************************
  void PrintTrackInfo(ND::THandle<ND::TReconBase> track) {
    TLorentzVector TrkPos;
    TVector3 TrkDir;
    double TrkMom = 0.0;
    double TrkChg = 0.0;
    ND::THandle<ND::TReconTrack> rTrack = track;
    ND::THandle<ND::TReconPID> rPID     = track;
    ND::THandle<ND::TReconBase> rPath   = track;
    if ( rTrack) {
      TrkPos = rTrack->GetPosition();
      TrkDir = rTrack->GetDirection();
      TrkMom = TrackingUtils::GetMomentum(track);
    } else if (rPID) {
      TrkPos = rPID->GetPosition();
      TrkDir = rPID->GetDirection();
      TrkMom = rPID->GetMomentum();
      TrkChg = rPID->GetCharge();
    }

    int Width = 23;
    std::cout<<std::setiosflags(std::ios::left);
    if( (rPath)->Get< ND::TIntegerDatum >("PathId") ){
      int PathId = (rPath)->Get< ND::TIntegerDatum >("PathId")->GetValue();
      std::cout<<"    -> Path ID: "<<PathId<<std::endl;
    }

    if (rPath->GetConstituents()){
      // One constituent means that we must pass the TReconTrack constituent layer.
      ND::THandle<ND::TReconBase> tmpBase;
      if (rPath->GetConstituents()->size() == 1){
        tmpBase = *(rPath->GetConstituents()->begin());
      }
      // More than one constituent should be only for cathode crossers.
      else {
        tmpBase = rPath;
      }

      if ( tmpBase->GetConstituents()){
        if (tmpBase->GetConstituents()->size() == 1){
          std::cout<<"  - 1 constituent with GasOutput ID:";
        } else {
          std::cout<<"  - "<<tmpBase->GetConstituents()->size()<<" constituents with GasOutput IDs:";
        }
        for (ND::TReconObjectContainer::const_iterator it=tmpBase->GetConstituents()->begin();it!=tmpBase->GetConstituents()->end(); it++){
          if( (*it)->Get< ND::TIntegerDatum >("PathId") ){
            int PathId = (*it)->Get< ND::TIntegerDatum >("PathId")->GetValue();
            std::cout<<" "<<PathId;
          }
        }
        std::cout<<std::endl;
      }
    }

    std::cout<<std::setw(Width)<<"  - Nb clusters        total     fitted"<<std::endl;
    if( rPath->Get< ND::TIntegerDatum >("NbVerticalClusters") && rPath->Get< ND::TIntegerDatum >("NbFittedVerticalClusters") ){
      int NbVerticalClu = rPath->Get< ND::TIntegerDatum >("NbVerticalClusters")->GetValue();
      int NbHorizontalClu = rPath->Get< ND::TIntegerDatum >("NbHorizontalClusters")->GetValue();
      int NbFittedVerticalClu = rPath->Get< ND::TIntegerDatum >("NbFittedVerticalClusters")->GetValue();
      int NbFittedHorizontalClu = rPath->Get< ND::TIntegerDatum >("NbFittedHorizontalClusters")->GetValue();
      std::cout<<std::setw(Width)<<"       Vertical:"  <<std::flush;
      std::cout<<std::setw(10)<<NbVerticalClu  <<NbFittedVerticalClu<<std::endl;
      std::cout<<std::setw(Width)<<"       Horizontal:"<<std::flush;
      std::cout<<std::setw(10)<<NbHorizontalClu<<NbFittedHorizontalClu<<std::endl;
    }

    if( rPath->Get< ND::TIntegerDatum >("T0Source") ){
      int T0Source = rPath->Get< ND::TIntegerDatum >("T0Source")->GetValue();
      std::cout<<std::setw(Width)<<"  - T0 source:"<<ConvertT0idxToName(T0Source)<<std::endl;
      if ( !T0Source){
        if( rPath->Get< ND::TRealDatum >("T0Range") ){
          double T0lowerlim = rPath->Get< ND::TRealDatum >("T0Range")->at(0);
          double T0upperlim = rPath->Get< ND::TRealDatum >("T0Range")->at(1);
          std::cout<<std::setw(Width)<<"       T0 limits:"<<"["<<T0lowerlim<<", "<<T0upperlim<<"]"<<std::endl;
        }
      }
    }
    std::cout<<std::setw(Width)<<"  - Status bits"<<std::endl;
    std::cout<<std::setw(Width)<<"       Chi2Fit:"<<rPath->CheckStatus(ND::TReconBase::kChi2Fit)<<std::endl;
    std::cout<<std::setw(Width)<<"       LikelihoodFit:"<<rPath->CheckStatus(ND::TReconBase::kLikelihoodFit)<<std::endl;
    std::cout<<std::setw(Width)<<"       KalmanFit:"<<rPath->CheckStatus(ND::TReconBase::kKalmanFit)<<std::endl;
    std::cout<<std::setw(Width)<<"       Success:"<<rPath->CheckStatus(ND::TReconBase::kSuccess)<<std::endl;
    if ( rPath->CheckStatus(ND::TReconBase::kLikelihoodFit) || rPath->CheckStatus(ND::TReconBase::kKalmanFit) || rPath->CheckStatus(ND::TReconBase::kChi2Fit)){
      ND::THandle< ND::TReconState > state = rPath->GetState();
      double ePosX = TMath::Sqrt(state->GetCovarianceValue(0,0));
      double ePosY = TMath::Sqrt(state->GetCovarianceValue(1,1));
      double ePosZ = TMath::Sqrt(state->GetCovarianceValue(2,2));
      double eDirX = TMath::Sqrt(state->GetCovarianceValue(4,4));
      double eDirY = TMath::Sqrt(state->GetCovarianceValue(5,5));
      double eDirZ = TMath::Sqrt(state->GetCovarianceValue(6,6));
      std::wstring pm = L"\u00B1";
      std::cout<<std::setw(Width)<<"  - Fit results"<<std::endl;
      std::cout<<std::setw(Width)<<"       Position:"<<"("<<TrkPos.X()<<", "<<TrkPos.Y()<<", "<<TrkPos.Z()<<", "<<TrkPos.T()<<")"<<std::endl;
      std::wcout<<std::setw(Width)<<""<<pm<<" ("<<ePosX<<", "<<ePosY<<", "<<ePosZ<<")"<<std::endl;
      //std::cout<<std::setw(Width)<<""<<std::setprecision(5)TrkPos.X()<<"|"<<TrkPos.Y()<<", "<<TrkPos.Z()<<", "<<TrkPos.T()<<")"<<std::endl;
      std::cout<<std::setw(Width)<<"       Direction:"<<"("<<TrkDir.X()<<", "<<TrkDir.Y()<<", "<<TrkDir.Z()<<")"<<std::endl;
      std::wcout<<std::setw(Width)<<""<<pm<<" ("<<eDirX<<", "<<eDirY<<", "<<eDirZ<<")"<<std::endl;
      std::cout<<std::setw(Width)<<"       Momentum:"<<TrkMom<<"  ";
      std::wcout<<pm;
      std::cout<<TrackingUtils::GetMomentumError(state)<<" MeV/c"<<std::endl;
      std::cout<<std::setw(Width)<<"       Charge:"<<TrkChg<<"  "<<std::endl;
      if( (rPath)->Get< ND::TRealDatum >("likFit_Sigma") && (rPath)->Get< ND::TRealDatum >("likFit_SigmaUnc") ){
        double Sigma = (rPath)->Get< ND::TRealDatum >("likFit_Sigma")->GetValue();
        double SigmaUnc = (rPath)->Get< ND::TRealDatum >("likFit_SigmaUnc")->GetValue();
        std::cout<<std::setw(Width)<<"       Sigma:"<<Sigma<<"  ";
        std::wcout<<pm;
        std::cout<<std::setw(Width)<<SigmaUnc<<std::endl;
      }
      if( (rPath)->Get< ND::TRealDatum >("FitLogLikelihood") ){
        std::cout<<std::setw(Width)<<"  -   Fit LogLikelihood:             "<<(rPath)->Get< ND::TRealDatum >("FitLogLikelihood")->GetValue()<<" = "<<(rPath)->Get< ND::TRealDatum >("FitLogLikelihoodX")->GetValue()<<" + "<<(rPath)->Get< ND::TRealDatum >("FitLogLikelihoodHV")->GetValue()<<std::endl;
        if( (rPath)->Get< ND::TIntegerDatum >("PathIdMatch") && (rPath)->Get< ND::TRealDatum >("PathMatchLklhd") ){
          for ( unsigned int i = 0; i < (rPath)->Get< ND::TIntegerDatum >("PathIdMatch")->size(); i++){
            std::cout<<std::setw(Width)<<"      Match across junction:"<<std::endl;
            std::cout<<std::setw(Width)<<"      - LogLikelihood: PathId "<<(rPath)->Get< ND::TIntegerDatum >("PathIdMatch")->at(i)<<" -> "<<(rPath)->Get< ND::TRealDatum >("PathMatchLklhd")->at(i)<<" = "<<(rPath)->Get< ND::TRealDatum >("PathMatchLklhdX")->at(i)<<" + "<<(rPath)->Get< ND::TRealDatum >("PathMatchLklhdHV")->at(i)<<std::endl;
          }
        }
        if( (rPath)->Get< ND::TIntegerDatum >("PatternIdMatch") && (rPath)->Get< ND::TRealDatum >("PatternMatchLklhd") ){
          for ( unsigned int i = 0; i < (rPath)->Get< ND::TIntegerDatum >("PatternIdMatch")->size(); i++){
            std::cout<<std::setw(Width)<<"      Match broken paths:"<<std::endl;
            std::cout<<std::setw(Width)<<"      - LogLikelihood: PatternId "<<(rPath)->Get< ND::TIntegerDatum >("PatternIdMatch")->at(i)<<", PathId "<<(rPath)->Get< ND::TIntegerDatum >("PatternPathIdMatch")->at(i)<<" -> "<<(rPath)->Get< ND::TRealDatum >("PatternMatchLklhd")->at(i)<<" = "<<(rPath)->Get< ND::TRealDatum >("PatternMatchLklhdX")->at(i)<<" + "<<(rPath)->Get< ND::TRealDatum >("PatternMatchLklhdHV")->at(i)<<std::endl;
          }
        }
      }
      std::cout<<std::setw(Width)<<"       chi2: "<<(rPath)->GetQuality()<<std::endl;
      std::cout<<std::setw(Width)<<"       ndof: "<<(rPath)->GetNDOF()<<std::endl;
    } else {
//\u00B1
      ND::THandle<ND::TComboHit> FirstClu = *(rPath->GetHits()->begin());
      ND::THandle<ND::TComboHit> LastClu  = *(rPath->GetHits()->rbegin());
      std::cout<<std::setw(Width)<<"  - Position"<<std::endl;
      std::cout<<std::setw(Width)<<"       First cluster:"<<"("<<FirstClu->GetPosition().X()<<", "<<FirstClu->GetPosition().Y()<<", "<<FirstClu->GetPosition().Z()<<", "<<FirstClu->GetTime()<<")"<<std::endl;
      std::cout<<std::setw(Width)<<"       Last cluster:"<<"("<<LastClu->GetPosition().X()<<", "<<LastClu->GetPosition().Y()<<", "<<LastClu->GetPosition().Z()<<", "<<LastClu->GetTime()<<")"<<std::endl;
    }

    if( rPath->Get< ND::TRealDatum >("Length") ){
      double Length = rPath->Get< ND::TRealDatum >("Length")->GetValue();
      std::cout<<std::setw(Width)<<"    Length:"<<Length<<" mm"<<std::endl;
    }

    std::cout<<std::setw(Width)<<"    Number of nodes:"<<rPath->GetNodes().size()<<std::endl;
    std::cout<<std::endl;
    if (rPID) {
      std::cout<<std::setw(Width)<<"  - PID:"<<rPID->ConvertParticleId()<<std::endl;
      if (rPID->Get< ND::TRealDatum >("tpcPid_PullMuon")){
        std::cout<<std::setw(Width)<<"       Pull Muon:" <<rPID->Get< ND::TRealDatum >("tpcPid_PullMuon")->GetValue()<<std::endl;
        std::cout<<std::setw(Width)<<"       Pull Pion:" <<rPID->Get< ND::TRealDatum >("tpcPid_PullPion")->GetValue()<<std::endl;
        std::cout<<std::setw(Width)<<"       Pull Electron:" <<rPID->Get< ND::TRealDatum >("tpcPid_PullEle")->GetValue()<<std::endl;
        std::cout<<std::setw(Width)<<"       Pull Proton:" <<rPID->Get< ND::TRealDatum >("tpcPid_PullProton")->GetValue()<<std::endl;
        std::cout<<std::setw(Width)<<"       Pull Kaon:" <<rPID->Get< ND::TRealDatum >("tpcPid_PullKaon")->GetValue()<<std::endl;
      }
      std::cout<<std::endl;
    }
    std::cout << std::resetiosflags(std::ios::left);
  }


  //**************************************************************
  void PrintConstituentMap(ND::THandle<ND::TReconBase> track, int level) {
    for (int i = 0; i < level; i++)
      std::cout<<" ";
    if( track->Get< ND::TIntegerDatum >("PathId") ){
      int PathId = track->Get< ND::TIntegerDatum >("PathId")->GetValue();
      std::cout<<"|-> "<<PathId<<std::endl;
    } else {
      std::cout<<"|-> no ID"<<std::endl;
    }

    if (track->GetConstituents()){
      for (ND::TReconObjectContainer::const_iterator it=track->GetConstituents()->begin();it!=track->GetConstituents()->end(); it++){
        PrintConstituentMap( *it, level+1);
      }
    }
  }

  //**************************************************************
  int GetIdFromOAEvent(ND::THandle<ND::TReconBase> track) {
    if( track->Get< ND::TIntegerDatum >("PathId") ){
      return track->Get< ND::TIntegerDatum >("PathId")->GetValue();
    }
    return 0;
  }

  //**************************************************************
  void ExtractPIDsFromPatterns(ND::THandle<ND::TReconObjectContainer> PatternContainer, ND::THandle<ND::TReconObjectContainer> PIDcontainer) {
    // Check that both containers are valid
    if (!PatternContainer || !PIDcontainer) return;
  
    for ( ND::TReconObjectContainer::iterator it=PatternContainer->begin(); it!=PatternContainer->end(); it++ ){
      ND::THandle<ND::TReconVertex> Pattern = (*it);
      ND::THandle<ND::TReconPID> Track = (*it);
      if (Pattern){
        // Loop over patterns
        for (ND::TReconObjectContainer::iterator ptc = Pattern->GetConstituents()->begin(); ptc != Pattern->GetConstituents()->end(); ptc++){
          ND::THandle<ND::TReconPID> Path = (*ptc);
          if (Path){
            // Use only the properly fitted tracks
            if (Path->CheckStatus(ND::TReconBase::kSuccess)){
              PIDcontainer->push_back(Path);
            }
          }
        }
      } else if (Track){
        PIDcontainer->push_back(Track);
      }
    }
  }


  // *********************************************************************************
  void HVClustersPrintout(ND::THandle<ND::THitSelection> inputClu, bool Extended){
    for (ND::THitSelection::const_iterator tmpClu = inputClu->begin() ; tmpClu != inputClu->end(); tmpClu++) {
      ND::THandle<ND::TTPCHVCluster> Cluster = *tmpClu;
      std::cout<<" * "<<Cluster->CalibX()<<", "<<Cluster->CalibY()<<", "<<Cluster->CalibZ()<<"   Vertical: "<<Cluster->IsVertical()<<std::endl;
      if( Extended){
        const ND::THitSelection allPads = Cluster->GetHits();
        for (ND::THitSelection::const_iterator tmpHit = allPads.begin(); tmpHit != allPads.end(); ++tmpHit) {
          ND::THandle<ND::TTPCHitPad> hitPad = *tmpHit;
          std::cout<<"           * T("<<hitPad->GetTime()<<", "<<hitPad->TimeFit()<<"), "<<hitPad->Y()<<", "<<hitPad->Z()<<", C("<<hitPad->ChargeFit()<<")"<<std::endl;
        }
      }
    }
  }


}
