#include "TTPCPattern.hxx"
#include "TTPCUtils.hxx"
#include "TTPCCalibration.hxx"

#include <TDatum.hxx>
#include <TIntegerDatum.hxx>
#include <TGeomInfo.hxx>

typedef std::vector<ND::THandle<ND::TReconBase> >::iterator constituent_iterator;

// *********************************************************************************
ND::TTPCPattern::TTPCPattern(){
  fUsable = true;
  fId = 0;
}


// *********************************************************************************
ND::TTPCPattern::TTPCPattern(ND::THandle<TTPCPath> path){
  AddConstituent(path);
  fUsable = true;
  fId = 0;
}


// *********************************************************************************
ND::TTPCPattern::~TTPCPattern() { }



// *********************************************************************************
void ND::TTPCPattern::SetId(unsigned int theId){
  fId = theId;
}


// *********************************************************************************
unsigned int ND::TTPCPattern::GetId(){
  return fId;
}


// *********************************************************************************
void ND::TTPCPattern::AddJunction(ND::THandle<TTPCJunction> junction){
  // TODO: Add exception cause this should not happen
  if ( GetNPaths() == 1 && GetNJunctions() == 0){
    std::cout<<"TTPCPattern: ERROR. This is a pattern with one clean track. You can't add more junctions"<<std::endl;
    return;
  }

  AddConstituent(junction);

  // Add the Paths that are constituents of the Junction to the constituent list
  // of the pattern. So make sure that the same path is not add twice.
  ND::THandle < ND::TReconObjectContainer > junctionConst = junction->GetConstituents();
  for (constituent_iterator jtc = junctionConst->begin(); jtc != junctionConst->end(); jtc++){
    bool alreadyHere = false;
    for (constituent_iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
      alreadyHere = ((*jtc) == (*ptc));
      if (alreadyHere)
        break;
    }
    if (!alreadyHere)
      AddConstituent(*jtc);
  }
    
}


// *********************************************************************************
void ND::TTPCPattern::InitialSetup(){
  if ( !fId){
    // TODO: proper exception
    std::cout<<"ERROR ND::TTPCPattern::InitialSetup called before assigning an ID to the pattern !"<<std::endl;
    throw;
  }

  ND::TReconBase::Status patDetector = GetDetectors();
  std::vector< ND::THandle<ND::TTPCPath> > Paths;
  std::vector< ND::THandle<ND::TTPCJunction> > Junctions;

  if (GetConstituents()) {
    // Calculate an average position and time.
    // This is just to have something "sensible" to put in the state.
    double meanX = 0.0;
    double meanY = 0.0;
    double meanZ = 0.0;
    double meanT = 0.0;
    double NbHits = 0.0;
    for (constituent_iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
      ND::THandle<ND::TTPCPath> path = (*ptc);
      ND::THandle<ND::TTPCJunction> junction = (*ptc);
      if(path){
        path->InitialSetup(fId);
        Paths.push_back(path);
      } else if (junction) {
        junction->InitialSetup();
        Junctions.push_back(junction);
      }
      
      ND::TReconBase::Status newDetector = (*ptc)->GetDetectors();
      if ( !patDetector){
        patDetector = newDetector;
        continue;
      } else if (newDetector != patDetector){
        // If the constituents are not all from the same TPC,
        // then throw exception.
        // TODO: proper exception
        throw;
      }
      for (ND::THitSelection::const_iterator hitit = (*ptc)->GetHits()->begin(); hitit != (*ptc)->GetHits()->end(); hitit++) {
        TVector3 Pos = (*hitit)->GetPosition();
        meanX += Pos.X();
        meanY += Pos.Y();
        meanZ += Pos.Z();
        meanT += (*hitit)->GetTime();
        NbHits += 1.;
      }
    }
    meanX /= NbHits;
    meanY /= NbHits;
    meanZ /= NbHits;
    meanT /= NbHits;
    ND::THandle<ND::TVertexState> tstate = this->GetState();
    tstate->SetPosition(TLorentzVector(meanX,meanY,meanZ,meanT));

    AddDetector(patDetector);
  } else {
    // TODO: proper exception
    throw;
  }

  // Add the junctions to the paths in the right order.
  // Loop over paths
  
  std::vector< ND::THandle<ND::TTPCPath> >::iterator pth = Paths.begin();
  for (; pth != Paths.end(); pth++){
    // Must clear the junction ids to erase previous connections when
    // recreating junctions during merging of patterns.
    (*pth)->ClearJunctionIds();

    std::vector< ND::THandle<ND::TTPCJunction> >::iterator jct = Junctions.begin();
    ND::THandle<ND::TTPCJunction> FirstJunction;
    ND::THandle<ND::TTPCJunction> LastJunction;
    // Loop over junctions
    for (; jct != Junctions.end(); jct++){
      TVector3 JuncPos = (*jct)->GetPosition().Vect();
      // Loop over the junction's paths to find a match
      for (constituent_iterator ptc = (*jct)->GetConstituents()->begin(); ptc != (*jct)->GetConstituents()->end(); ptc++){
        ND::THandle<ND::TTPCPath> tmpPath = *ptc;
        if ( !tmpPath)
          continue;
        if ( (*pth)->GetId() == tmpPath->GetId()){
          // Compare with the position of the first and last cluster and fill FirstJunction and LastJunction accordingly
          ND::THandle<ND::TTPCPath> Path = *pth;
		  ND::THandle<ND::TTPCHVCluster> tmpClu = *(Path->GetHits()->begin());
		  const TVector3 FirstPos = tmpClu->GetCalibPosition();
		  tmpClu = *(Path->GetHits()->rbegin());
		  const TVector3 LastPos = tmpClu->GetCalibPosition();
          TVector3 DistToFirst = JuncPos - FirstPos;
          TVector3 DistToLast = JuncPos - LastPos;
          if ( DistToLast.Mag() > DistToFirst.Mag()){
            FirstJunction = (*jct);
          } else {
            LastJunction = (*jct);
          }
          break;
        }
      }
    }
    if (FirstJunction){
      (*pth)->AddJunctionId(FirstJunction->GetId());
      (*pth)->SetFrontConnection(true);
    }
    if (LastJunction){
      (*pth)->AddJunctionId(LastJunction->GetId());
      (*pth)->SetBackConnection(true);
    }
  }

}

// *********************************************************************************
unsigned int ND::TTPCPattern::GetNPaths(){
  unsigned int nbPaths = 0;
  ND::THandle<ND::TReconObjectContainer> constituents = GetConstituents();
  if (constituents){
    for (constituent_iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
      ND::THandle< ND::TTPCPath > tmpRTrack= (*ptc);
      if ( tmpRTrack)
        nbPaths++;
    }
  }
  return nbPaths;
}


// *********************************************************************************
unsigned int ND::TTPCPattern::GetNJunctions(){
  unsigned int nbJunctions = 0;
  ND::THandle<ND::TReconObjectContainer> constituents = GetConstituents();
  if (constituents){
    for (constituent_iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
      ND::THandle< ND::TTPCJunction > tmpJunction= (*ptc);
      if ( tmpJunction)
        nbJunctions++;
    }
  }
  return nbJunctions;
}


// *********************************************************************************
std::vector< ND::THandle<ND::TTPCPath> > ND::TTPCPattern::GetPaths(){
  std::vector< ND::THandle<ND::TTPCPath> > Paths;
  ND::THandle<ND::TReconObjectContainer> constituents = GetConstituents();
  if (!constituents)
    return Paths;
  for (constituent_iterator ptc = constituents->begin(); ptc != constituents->end(); ptc++){
    ND::THandle< ND::TTPCPath > tmpPath= (*ptc);
    if ( tmpPath)
      Paths.push_back(tmpPath);
  }
  return Paths;
}


// *********************************************************************************
std::vector< ND::THandle<ND::TTPCJunction> > ND::TTPCPattern::GetJunctions(){
  std::vector< ND::THandle<ND::TTPCJunction> > Junctions;
  ND::THandle<ND::TReconObjectContainer> constituents = GetConstituents();
  if (!constituents)
    return Junctions;
  for (constituent_iterator ptc = constituents->begin(); ptc != constituents->end(); ptc++){
    ND::THandle< ND::TTPCJunction > tmpJunction= (*ptc);
    if ( tmpJunction)
      Junctions.push_back(tmpJunction);
  }
  return Junctions;
}


// *********************************************************************************
ND::THandle<ND::THitSelection> ND::TTPCPattern::GetConstituentHits(){
  ND::THandle<ND::THitSelection> chits(new ND::THitSelection());
  ND::THandle<ND::TReconObjectContainer> constituents = GetConstituents();
  if (constituents){
    for (constituent_iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
      for (ND::THitSelection::const_iterator hitit = (*ptc)->GetHits()->begin(); hitit != (*ptc)->GetHits()->end(); hitit++) {
        chits->push_back(*hitit);
      }
    }
  }
  return chits;
}


// *********************************************************************************
void ND::TTPCPattern::SetUsable(bool usable){
  fUsable = usable;
}


// *********************************************************************************
bool ND::TTPCPattern::IsUsable(){
  return fUsable;
}


// *********************************************************************************
void ND::TTPCPattern::SetT0(TTPCT0 &T0){
  fT0 = T0;
  ND::THandle<ND::TVertexState> tstate = this->GetState();
  TLorentzVector PosTime = tstate->GetPosition();
  ND::THandle<ND::THit> hit = *((*(GetConstituents()->begin()))->GetHits()->begin());
  double Sense = int(ND::TGeomInfo::Get().TPC().GetDriftSense( hit->GetGeomId() ));
  double RPx = hit->GetPosition().X();
  // How far are we from the read out plane ?
  double DriftDistance = (PosTime.T() - fT0.GetT0() - ND::tpcCalibration().GetTimeOffset()) * ND::tpcCalibration().GetDriftVelocity();
  PosTime.SetX( RPx - (Sense * DriftDistance));
  PosTime.SetT( T0.GetT0());
  tstate->SetPosition(PosTime); 

  for (ND::TReconObjectContainer::iterator constit = GetConstituents()->begin(); constit != GetConstituents()->end(); constit++) {
    ND::THandle<ND::TTPCPath> path = (*constit);
    ND::THandle<ND::TTPCJunction> junction = (*constit);
    if (path){
      path->SetT0(T0);
    } else if (junction){
      junction->SetT0(T0);
    }
  }
}

// *********************************************************************************
bool ND::TTPCPattern::HasT0(){
  return (fT0.GetSource() != kNoT0src);
}

// *********************************************************************************
TTRExT0Source ND::TTPCPattern::GetT0Source(){
  return fT0.GetSource();
}

// *********************************************************************************
double ND::TTPCPattern::GetT0(){
  return fT0.GetT0();
}



// *********************************************************************************
ND::THandle<ND::TReconBase> ND::TTPCPattern::ConvertToOAEvent() {
  ND::THandle<ND::TReconBase> Output;
  if ( GetNPaths() == 1 && GetNJunctions() == 0){
    ND::THandle<ND::TTPCPath> path = *(GetConstituents()->begin());
    // Could be TReconTrack or a TReconPID
    ND::THandle<TReconBase> rtrack = path->ConvertToOAEvent(); 
    Output = rtrack;
  } else {
    ND::THandle<ND::TReconVertex> tvertex = TTPCUtils::CopyCreateTReconVertex( this);
    // To differentiate with a "regular" TReconVertex, let's mark this container as composite.
    tvertex->SetComposite(true);
    // Some temporary containers for quick access
    std::vector< ND::THandle<ND::TTPCPath> > oldPaths;
    std::vector< ND::THandle<ND::TTPCJunction> > oldJunctions;
    std::vector< ND::THandle<ND::TReconBase> > newPaths;
    std::vector< ND::THandle<ND::TReconVertex> > newJunctions;
    // Convert the constituents to their oaEvent counterpart
    for (constituent_iterator ptc = GetConstituents()->begin(); ptc != GetConstituents()->end(); ptc++){
      ND::THandle<ND::TTPCPath> path = (*ptc);
      ND::THandle<ND::TTPCJunction> junction = (*ptc);
      if (path){
        THandle<TReconBase> rtrack = path->ConvertToOAEvent(); 
        oldPaths.push_back(path);
        newPaths.push_back(rtrack);
        tvertex->AddConstituent(rtrack);
      } else if (junction){
        THandle<TReconVertex> rvertex = junction->ConvertToOAEvent(); 
        oldJunctions.push_back(junction);
        newJunctions.push_back(rvertex);
        tvertex->AddConstituent(rvertex);
        ND::THandle<ND::THitSelection> vertHits = rvertex->GetHits();
      }
    }
    // The new TReconVertex (Junction) has no constituents.
    // Fill the right freshly converted TReconTrack or TReconPID.
    std::vector< ND::THandle<ND::TTPCJunction> >::iterator ojt = oldJunctions.begin();
    std::vector< ND::THandle<ND::TReconVertex> >::iterator nvt = newJunctions.begin();
    for (; ojt != oldJunctions.end(); ojt++, nvt++){
      ND::THandle<ND::TTPCJunction> oldJunc = (*ojt);
      ND::THandle<ND::TReconVertex> newJunc = (*nvt);

      for (constituent_iterator jtc = oldJunc->GetConstituents()->begin(); jtc != oldJunc->GetConstituents()->end(); jtc++){
        ND::THandle<ND::TTPCPath> path = (*jtc);
        bool foundNewPath = false;
        for (std::vector< ND::THandle<ND::TReconBase> >::iterator ptc = newPaths.begin(); ptc != newPaths.end(); ptc++){
          if( !(*ptc)->Get< ND::TIntegerDatum >("PathId")){
            // TODO: exception, this should not happen
            throw;
          }
          if ( int(path->GetId()) == (*ptc)->Get< ND::TIntegerDatum >("PathId")->GetValue()){
            newJunc->AddConstituent(*ptc);
            foundNewPath = true;
          }
        }
        if (foundNewPath) continue;
      }

    }

    Output = tvertex;
  }

  ND::TIntegerDatum* patternId = new ND::TIntegerDatum("PatternId",fId);
  Output->AddDatum(patternId);

  return Output;
}

