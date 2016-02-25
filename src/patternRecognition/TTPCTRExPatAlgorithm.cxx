// eddy
#include "TTPCTRExPatAlgorithm.hxx"

ND::TTPCTRExPatAlgorithm::TTPCTRExPatAlgorithm() {
  // has no hits by default
  fDriftVelocity = 0;
  fHasHits = false;

  // initial values
  fSubAlgorithms = std::vector<ND::TTPCTRExPatSubAlgorithm*>();
  fMasterHitMap = std::map<long, ND::TTPCUnitVolume*>();
  fMasterHitMapHigh = std::map<long, ND::TTPCUnitVolume*>();
  fMasterLayout = 0;
  fMasterVolGroupMan = 0;
  fMasterVolGroupManHigh = 0;
  fMasterAStar = 0;
  fFeatureFinder = 0;
}
ND::TTPCTRExPatAlgorithm::~TTPCTRExPatAlgorithm(){
  // clear values from last lot of processing
  CleanUp();
}

void ND::TTPCTRExPatAlgorithm::CleanUp(){
  // delete stuff from previous processing
  for(std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator it = fSubAlgorithms.begin(); it != fSubAlgorithms.end(); ++it){
    delete *it;
  };
  fSubAlgorithms.clear();

  // delete all unit volumes
  for(std::map<long, ND::TTPCUnitVolume*>::iterator vol = fMasterHitMap.begin(); vol != fMasterHitMap.end(); ++vol){
    delete vol->second;
  };
  fMasterHitMap.clear();
  fMasterHitMapHigh.clear();

  // delete stuff for whole event
  if(fMasterLayout) 
  {
    delete fMasterLayout;
    fMasterLayout = 0;
  }
  if(fMasterVolGroupMan)
  {
    delete fMasterVolGroupMan;
    fMasterVolGroupMan = 0;
  }
  if(fMasterVolGroupManHigh)
  {
    delete fMasterVolGroupManHigh;
    fMasterVolGroupManHigh = 0;
  }
  if(fMasterAStar)
  {
    delete fMasterAStar;
    fMasterAStar = 0;
  }
  if(fFeatureFinder)
  {
    delete fFeatureFinder;
    fFeatureFinder = 0;
  }

  //clean up THandles
  fSubEvents.clear();
  fDeltaHits.clear();
  fHits = ND::THandle<ND::THitSelection>();
}

void ND::TTPCTRExPatAlgorithm::PrepareHits(ND::THandle<ND::THitSelection> hits){
  // needed for pattern recognition - min and max times either side of the cathode
  double tPMin = +99999999.;
  double tPMax = -99999999.;
  double tNMin = +99999999.;
  double tNMax = -99999999.;
  for (ND::THitSelection::iterator hit = hits->begin(); hit != hits->end(); ++hit){
    ND::THandle<ND::TTPCHitPad> phit = *hit;
    if(!phit) continue;

    // add to selection
    fHits->AddHit(phit);

    // sense determines which side of the cathode we're on
    int curX = int(ND::TGeomInfo::TPC().GetDriftSense(phit->GetGeomId()));

    // needed for pattern recognition - get min and max times either side of the cathode
    std::vector<double> peakTimes = phit->GetPeakTimes();
    if(peakTimes.size()){
      for(std::vector<double>::iterator peakTimeIt = peakTimes.begin(); peakTimeIt != peakTimes.end(); ++peakTimeIt){
        if(curX < 0){
          tNMin = std::min(tNMin, *peakTimeIt);
          tNMax = std::max(tNMax, *peakTimeIt);
        }
        else{
          tPMin = std::min(tPMin, *peakTimeIt);
          tPMax = std::max(tPMax, *peakTimeIt);
        };
      };
    }
    else{
      if(curX < 0){
        tNMin = std::min(tNMin, phit->GetTime());
        tNMax = std::max(tNMax, phit->GetTime());
      }
      else{
        tPMin = std::min(tPMin, phit->GetTime());
        tPMax = std::max(tPMax, phit->GetTime());
      };
    };
  };

  // add maximum and minimum times either side of the cathode
  fMasterLayout->SetTimeRanges(tNMin, tNMax, tPMin, tPMax);

  // work out minimum and maximum cell ids in x, y and z
  int minX = +99999;
  int maxX = -99999;
  int minY = +99999;
  int maxY = -99999;
  int minZ = +99999;
  int maxZ = -99999;
  for (ND::THitSelection::iterator hitIt = hits->begin(); hitIt != hits->end(); ++hitIt){
    ND::THandle<ND::TTPCHitPad> hit = *hitIt;
    if(!hit) continue;

    // convert position to cell id in x, y and z
    ND::TTPCCellInfo3D cell = fMasterLayout->GetPadPosID(hit);

    if (cell.x < 0 || cell.y < 0 || cell.z < 0) { 
      if(ND::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << "ERROR: missing cell, returned (" << cell.x << ", " << cell.y << ", " << cell.z << ")" << std::endl;
      continue;
    };

    // find minima and maxima
    minX = std::min(minX, cell.x);
    maxX = std::max(maxX, cell.x);
    minY = std::min(minY, cell.y);
    maxY = std::max(maxY, cell.y);
    minZ = std::min(minZ, cell.z);
    maxZ = std::max(maxZ, cell.z);
  };

  // determine if any hits have been added
  fHasHits = (maxX >= 0 && maxY >= 0 && maxZ >= 0);
  if (!fHasHits) return;

  // add minimum and maximum cell ids in x, y and z to layout
  fMasterLayout->SetRanges(minX,maxX, minY,maxY, minZ,maxZ);

  // loop over positions and charges and add a new pattern recognition cell for each one
  for (ND::THitSelection::iterator hitIt = hits->begin(); hitIt != hits->end(); ++hitIt){
    ND::THandle<ND::TTPCHitPad> hit = *hitIt;
    ND::THandle<ND::TMultiHit> hitAsMulti = hit->GetMultiHit();
    if(!hit) continue;

    // get MM information
    ND::TGeometryId geomId = hit->GetGeomId();
    unsigned int tpc = ND::TGeomInfo::Get().TPC().GeomIdToTPC(geomId);
    unsigned int half = ND::TGeomInfo::Get().TPC().GeomIdToHalf(geomId);
    unsigned int mm = ND::TGeomInfo::Get().TPC().GeomIdToMM(geomId);

    // get hit fec and asic information
    ND::TChannelId chan = hitAsMulti->GetChannelId();
    ND::TTPCChannelId tpcChan(chan);
    unsigned int fec = tpcChan.GetFEC();
    unsigned int asic = tpcChan.GetAsic();

    // convert position to cell id in x, y and z
    ND::TTPCCellInfo3D cell = fMasterLayout->GetPadPosID(hit, 0);

    // get unique ASIC region for hairy tracks
    int asicRegionY = cell.y / fMasterLayout->GetASICSplittingY();
    int asicRegionZ = cell.z / fMasterLayout->GetASICSplittingZ();

    // convert cell id in x, y and z to unique id
    long id = fMasterLayout->Mash(cell.x, cell.y, cell.z);

    // ignore cells below or above minima or maxima
    if (cell.x < minX || cell.y < minY || cell.z < minZ) continue;
    if (cell.x > maxX || cell.y > maxY || cell.z > maxZ) continue;

    // see if a cell already exists at this position
    std::map<long, ND::TTPCUnitVolume*>::iterator el = fMasterHitMap.find(id);
    // if cell doesn't already exist, define a new one at this position
    if(el == fMasterHitMap.end()){
      ND::TTPCUnitVolume* curVol = new ND::TTPCUnitVolume();

      curVol->SetCell(cell.x, cell.y, cell.z, cell.edgeX, cell.edgeY, cell.edgeZ, id);
      curVol->SetAux(cell.segX, cell.segY, cell.segZ);
      curVol->SetMMLoc(tpc, half, mm);
      curVol->SetFECASIC(fec, asic);
      curVol->SetRegion(asicRegionY, asicRegionZ);
      if(cell.segX > 0){
        curVol->SetTimeOffset(fMasterLayout->GetTPMin());
      }
      else{
        curVol->SetTimeOffset(fMasterLayout->GetTNMin());
      };
      fMasterHitMap[id] = curVol;
    };
    // increment charge and average position at this cell
    fMasterHitMap[id]->AddEvent(hit);
  };
}
void ND::TTPCTRExPatAlgorithm::PopulateDeltaHits(){
  ND::THandle< ND::TTPCVolGroup > deltaHits (new ND::TTPCVolGroup(fMasterLayout) );

  for(std::map<long, ND::TTPCUnitVolume*>::iterator vol = fMasterHitMap.begin(); vol != fMasterHitMap.end(); ++vol){
    int nPeaks = vol->second->GetNPeaksSum();
    if(nPeaks != 1) deltaHits->AddHit(vol->second);
  };

  std::vector< ND::THandle<ND::TTPCVolGroup> > newDeltaHitGroups = fMasterVolGroupMan->GroupDeltaHits(deltaHits);

  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator newDeltaHitGroupIt = newDeltaHitGroups.begin(); newDeltaHitGroupIt != newDeltaHitGroups.end(); ++newDeltaHitGroupIt){
    ND::THandle<ND::TTPCVolGroup> newDeltaHitGroup = *newDeltaHitGroupIt;
    fDeltaHits.push_back(newDeltaHitGroup);
  };
}
void ND::TTPCTRExPatAlgorithm::FillInfoLowCharge(){
  double chargeCut = fMasterLayout->GetChargeCut();
  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    if(mapEl->second->GetQMax() >= chargeCut){
      mapEl->second->SetLowChargeTagged(false);
    }
    else{
      mapEl->second->SetLowChargeTagged(true);
    };
  };
}
void ND::TTPCTRExPatAlgorithm::FillInfoFullASIC(){
  std::map< std::pair<int, int>, int> subAsicOccupancies;
  std::map<long, int> asicOccupancies;
  std::set< std::pair<int, int> > foundPads;

  // increment total ASIC charge each
  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    ND::TTPCUnitVolume* vol = mapEl->second;
    long asicID = vol->GetUniqueASICID();
    std::pair<int, int> padSegment = std::make_pair(vol->GetASICRegionY(), vol->GetASICRegionZ());
    std::pair<int, int> padPos = std::make_pair(vol->GetY(), vol->GetZ());

    if(asicOccupancies.find(asicID) == asicOccupancies.end()){
      asicOccupancies[asicID] = 0;
    };
    if(subAsicOccupancies.find(padSegment) == subAsicOccupancies.end()){
      subAsicOccupancies[padSegment] = 0;
    };
    if(!foundPads.count(padPos)){
      asicOccupancies[asicID] ++;
      subAsicOccupancies[padSegment] ++;
      foundPads.insert(padPos);
    };
  };

  // find any with too high an occupancy
  std::set<long> dodgyASICs = std::set<long>();
  for(std::map<long, int>::iterator occupancyEl = asicOccupancies.begin(); occupancyEl != asicOccupancies.end(); ++occupancyEl){
    if(occupancyEl->second >= fMasterLayout->GetASICOccupancyCut()){
      dodgyASICs.insert(occupancyEl->first);
    };
  };
  std::set< std::pair<int, int> > dodgySubASICs = std::set< std::pair<int, int> >();
  for(std::map< std::pair<int, int>, int>::iterator occupancyEl = subAsicOccupancies.begin(); occupancyEl != subAsicOccupancies.end(); ++occupancyEl){
    if(occupancyEl->second >= fMasterLayout->GetASICSubOccupancyCut()){
      dodgySubASICs.insert(occupancyEl->first);
    };
  };

  // attach these dodgy hits to groups
  ND::THandle<ND::TTPCVolGroup> dodgyASICHits (new ND::TTPCVolGroup(fMasterLayout));

  // attach this information to hits
  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    ND::TTPCUnitVolume* vol = mapEl->second;
    int regY = vol->GetASICRegionY();
    int regZ = vol->GetASICRegionZ();
    // surrounding groups will be dodgy as well
    for(int i=-1; i<=1; ++i)
    for(int j=-1; j<=1; ++j){
      if(dodgySubASICs.count( std::make_pair(regY+i, regZ+j) )){
        dodgyASICHits->AddHit(vol);
      };
    };
    // also add any ASICs with too many hits
    if(dodgyASICs.count(vol->GetUniqueASICID())){
      dodgyASICHits->AddHit(vol);
    };
  };

  // expand these groups
  int asicOccExpansion = fMasterLayout->GetASICOccExpansion();

  // use A* algorithm to associate extra hits
  fMasterAStar->MergeBestHits(fMasterVolGroupMan, dodgyASICHits, false, true, asicOccExpansion);

  // add expanded hits to group
  for(std::map<long, ND::TTPCUnitVolume*>::iterator dodgyASICHitEl = dodgyASICHits->begin(); dodgyASICHitEl != dodgyASICHits->end(); ++dodgyASICHitEl){
    ND::TTPCUnitVolume* vol = dodgyASICHitEl->second;
    vol->SetFullASICTagged(true);
  };
}
void ND::TTPCTRExPatAlgorithm::FillInfoSatASIC(){
  std::map<long, int> asicSatOccupancies;

  // increment total ASIC charge each
  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    ND::TTPCUnitVolume* vol = mapEl->second;
    long asicId = vol->GetUniqueASICID();

    if(asicSatOccupancies.find(asicId) == asicSatOccupancies.end()){
      asicSatOccupancies[asicId] = 0;
    };
    asicSatOccupancies[asicId] += vol->GetNSaturated();
  };

  // find any with above threshold number of saturated hits
  std::set<long> satASICs;
  for(std::map<long, int>::iterator satOccupancyEl = asicSatOccupancies.begin(); satOccupancyEl != asicSatOccupancies.end(); ++satOccupancyEl){
    if(satOccupancyEl->second >= fMasterLayout->GetASICSaturationCut()){
      satASICs.insert(satOccupancyEl->first);
    };
  };

  // attach these dodgy hits to groups
  ND::THandle<ND::TTPCVolGroup> satASICHits (new ND::TTPCVolGroup(fMasterLayout));

  // attach this information to hits
  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    ND::TTPCUnitVolume* vol = mapEl->second;
    if(satASICs.count( vol->GetUniqueASICID() )){
      satASICHits->AddHit(vol);
    };
  };

  // expand these groups
  int asicSatExpansion = fMasterLayout->GetASICSatExpansion();

  // use A* algorithm to associate extra hits
  fMasterAStar->MergeBestHits(fMasterVolGroupMan, satASICHits, false, true, asicSatExpansion);

  // add expanded hits to group
  for(std::map<long, ND::TTPCUnitVolume*>::iterator satASICHitEl = satASICHits->begin(); satASICHitEl != satASICHits->end(); ++satASICHitEl){
    ND::TTPCUnitVolume* vol = satASICHitEl->second;
    vol->SetSatASICTagged(true);
  };
}
void ND::TTPCTRExPatAlgorithm::FillInfoNegativePeak(){
  double earlyNegativeCut = fMasterLayout->GetEarlyNegativePeakCut();
  double lateNegativeCut = fMasterLayout->GetLateNegativePeakCut();

  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    ND::TTPCUnitVolume* vol = mapEl->second;

    double earlyRatio = -(vol->GetNegativePeakEarly()/vol->GetQMax());
    double lateRatio = -(vol->GetNegativePeakLate()/vol->GetQMax());

    if(earlyRatio > earlyNegativeCut){
      vol->SetEarlyNegativeTagged(true);
    }
    else if(lateRatio > lateNegativeCut){
      vol->SetLateNegativeTagged(true);
    };
  };
}
void ND::TTPCTRExPatAlgorithm::FillHighQualityHits(){
  fMasterHitMapHigh = std::map<long, ND::TTPCUnitVolume*>();

  for(std::map<long, ND::TTPCUnitVolume*>::iterator mapEl = fMasterHitMap.begin(); mapEl != fMasterHitMap.end(); ++mapEl){
    long id = mapEl->first;
    ND::TTPCUnitVolume* vol = mapEl->second;

    if(vol->GetPathology()){
      // fill hit pads with information on status as hairy candidate
      for(std::vector< ND::THandle<ND::TTPCHitPad> >::iterator hitPadIt = vol->GetHitsBegin(); hitPadIt != vol->GetHitsEnd(); ++hitPadIt){
        ND::THandle<ND::TTPCHitPad> hitPad = *hitPadIt;
        hitPad->SetHairCandidate(true);
      };
    }
    else{
      // populate high quality hit map
      fMasterHitMapHigh[id] = vol;
    };
  };
}

void ND::TTPCTRExPatAlgorithm::FillUsedUnusedHits(ND::THitSelection* usedTREx, ND::THitSelection* used, ND::THitSelection* unused){
  ND::THitSelection* unusedTREx = new ND::THitSelection();

  // add all hits
  for(ND::THitSelection::const_iterator hitIt = fHits->begin(); hitIt != fHits->end(); ++hitIt){
    ND::THandle<ND::TTPCHitPad> hit = *hitIt;
    if(hit) unusedTREx->AddHit(hit);
  };
  // remove those contained in used
  for(ND::THitSelection::const_iterator hitIt = usedTREx->begin(); hitIt != usedTREx->end(); ++hitIt){
    ND::THandle<ND::TTPCHitPad> hit = *hitIt;
    if(hit) unusedTREx->RemoveHit(hit);
  };

  // save OA versions in containers
  for(ND::THitSelection::const_iterator hitIt = usedTREx->begin(); hitIt != usedTREx->end(); ++hitIt){
    ND::THandle<ND::TTPCHitPad> hitPad = *hitIt;
    ND::THandle<ND::TMultiHit> multiHit = hitPad->ConvertToOAEvent();
    if(multiHit) used->AddHit(multiHit);
  };
  for(ND::THitSelection::const_iterator hitIt = unusedTREx->begin(); hitIt != unusedTREx->end(); ++hitIt){
    ND::THandle<ND::TTPCHitPad> hitPad = *hitIt;
    ND::THandle<ND::TMultiHit> multiHit = hitPad->ConvertToOAEvent();
    if(multiHit) unused->AddHit(multiHit);
  };

  delete unusedTREx;
}
void ND::TTPCTRExPatAlgorithm::GetPatterns(ND::TReconObjectContainer *foundPatterns){
  // add patterns from all sub events
  for(std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator subAlgIt = fSubAlgorithms.begin(); subAlgIt != fSubAlgorithms.end(); ++subAlgIt){
    ND::TTPCTRExPatSubAlgorithm* subAlg = *subAlgIt;
    ND::THandle<ND::TTPCPattern> foundPattern = subAlg->GetPattern();

    // ensure pattern exists and contains sensible numbers of paths and junctions
    if(!foundPattern) continue;

    int nPaths = foundPattern->GetNPaths();
    int nJunctions = foundPattern->GetNJunctions();

    bool errorInPattern = false;
    if(nPaths < 1) errorInPattern = true;
    if( (nPaths < 2 && nJunctions > 0) || (nPaths > 1 && nJunctions < 1) ) errorInPattern = true;
    if(errorInPattern){
      if(ND::tpcDebug().PatternRecognition(DB_ERROR)) std::cout << " WARNING: pattern has " << nPaths << " paths and " << nJunctions << " junctions - something went wrong! " << std::endl;
    }
    else{
      foundPatterns->push_back(foundPattern);
    };
  };

  // print debug output
  VerifyHitsFull(foundPatterns);
  VerifyUsedUnusedFull(foundPatterns);
}

void ND::TTPCTRExPatAlgorithm::Process(const ND::TAlgorithmResult &in, ND::THitSelection* used, ND::THitSelection* unused){
  ND::THandle<ND::THitSelection> hits = in.GetHitSelection();

  Process(hits, used, unused);
}
void ND::TTPCTRExPatAlgorithm::Process(ND::THandle<ND::THitSelection> hits, ND::THitSelection* used, ND::THitSelection* unused){
  // setting up
  fSubAlgorithms = std::vector<ND::TTPCTRExPatSubAlgorithm*>();
  // master layout for all sub-events
  fMasterLayout = new ND::TTPCLayout();
  fDriftVelocity = fMasterLayout->GetDriftSpeed();
  // map of all hits in this event
  fMasterHitMap = std::map<long, ND::TTPCUnitVolume*>();
  // set up primary group of hits
  fHits = ND::THandle<ND::THitSelection>(new THitSelection());

  // reset group IDs
  ND::TTPCVolGroup::ResetFreeID();
  // prepare hits
  PrepareHits(hits);

  // master manager for all unit volumes
  fMasterVolGroupMan = new ND::TTPCVolGroupMan(fMasterLayout);
  fMasterVolGroupMan->AddPrimaryHits(fMasterHitMap);

  // master A* algorithm for finding paths
  fMasterAStar = new ND::TTPCAStar(fMasterLayout);
  fMasterAStar->AddHits(fMasterVolGroupMan, fMasterHitMap);

  // find all hits with pathalogical hairy behaviour
  if(fMasterLayout->GetUsePatRecPathologyCut() > 0){
    FillInfoLowCharge();
    FillInfoFullASIC();
    FillInfoSatASIC();
    FillInfoNegativePeak();
    FillHighQualityHits();
  };

  // master manager for all unit volumes passing a high charge cut
  if(fMasterLayout->GetUsePatRecPathologyCut() == 1){
    fMasterVolGroupManHigh = new ND::TTPCVolGroupMan(fMasterLayout);
    fMasterVolGroupManHigh->AddPrimaryHits(fMasterHitMapHigh);
  };

  /* REMOVED DUE TO HAIRY TRACK ISSUES
   * the following lines of code tag delta hits from their multi peak waveforms
   * they have been found to be sensitive to hairy tracks and have thus been removed
   * it may be useful to reinstate them when selecting for hairiness
  // fill fDeltaHits container with any delta hits
  PopulateDeltaHits();

  // mark delta hits from container
  for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator deltaHitGroupIt = fDeltaHits.begin(); deltaHitGroupIt != fDeltaHits.end(); ++deltaHitGroupIt){
    ND::THandle<ND::TTPCVolGroup> deltaHitGroup = *deltaHitGroupIt;
    for(std::map<long, ND::TTPCUnitVolume*>::iterator deltaHitEl = deltaHitGroup->begin(); deltaHitEl != deltaHitGroup->end(); ++deltaHitEl){
      deltaHitEl->second->SetDeltaTagged(true);
    };
  };
  */

  // split all hits up into lists of sub events, with separate group for high charge ones if needed
  std::vector< ND::THandle<ND::TTPCVolGroup> > subEvents;
  std::vector< ND::THandle<ND::TTPCVolGroup> > subEventsHigh;
  if(fMasterLayout->GetUsePatRecPathologyCut() == 1){
    subEvents = fMasterVolGroupManHigh->GetConnectedHits(ND::TTPCConnection::path, ND::TTPCHitGroupings::all, true);

    // add initial sub events to high charge selection
    for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator subEventIt = subEvents.begin(); subEventIt != subEvents.end(); ++subEventIt){
      ND::THandle<ND::TTPCVolGroup> subEventHigh (new ND::TTPCVolGroup(fMasterLayout));
      subEventHigh->AddHits(*subEventIt);
      subEventsHigh.push_back(subEventHigh);
    };

    // use A* algorithm to associate required extra hits
    fMasterAStar->MergeBestHits(fMasterVolGroupMan, subEvents);
  }
  else{
    subEvents = fMasterVolGroupMan->GetConnectedHits(ND::TTPCConnection::path);
    if(fMasterLayout->GetUsePatRecPathologyCut() > 0){
      // just add selection of high charged hits to high charged hit container
      for(std::vector< ND::THandle<ND::TTPCVolGroup> >::iterator subEventIt = subEvents.begin(); subEventIt != subEvents.end(); ++subEventIt){
        ND::THandle<ND::TTPCVolGroup> subEvent = *subEventIt;

        ND::THandle<ND::TTPCVolGroup> subEventHigh (new ND::TTPCVolGroup(fMasterLayout));
        for(std::map<long, ND::TTPCUnitVolume*>::iterator volEl = subEvent->begin(); volEl != subEvent->end(); ++volEl){
          if(!volEl->second->GetPathology()){
            subEventHigh->AddHit(volEl->second);
          };
        };
        subEventsHigh.push_back(subEventHigh);
      };
    }
    else{
      subEventsHigh = subEvents;
    };
  };

  // push all groups of decent size into sub events
  for(unsigned int i=0; i<subEvents.size(); ++i){
    ND::THandle<ND::TTPCVolGroup> subEvent = subEvents.at(i);
    ND::THandle<ND::TTPCVolGroup> subEventHigh = subEventsHigh.at(i);

    if(fMasterVolGroupMan->CheckUsability(subEventHigh)){
      // create sub-algorithm for each sub-event
      ND::TTPCTRExPatSubAlgorithm* subAlg = new ND::TTPCTRExPatSubAlgorithm(fMasterLayout);

      // set hit selection for sub-event
      subAlg->SetUpHits(subEvent->GetHitMap(), fMasterAStar);
      subAlg->SetUpHitsHigh(subEventHigh->GetHitMap(), fMasterAStar);

      // add sub-algorithm to vector
      fSubAlgorithms.push_back(subAlg);
      fSubEvents.push_back(subEvent);
    };
  };

  // first pass of processing
  int subEvent = 0;
  for(std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator algIt = fSubAlgorithms.begin(); algIt != fSubAlgorithms.end(); ++algIt){
    ND::TTPCTRExPatSubAlgorithm* alg = *algIt;
    subEvent++;
    if(ND::tpcDebug().PatternRecognition(DB_INFO)) std::cout << " ----- pattern recognition for sub-event " << subEvent << " ------ " << std::endl;
    // process each sub-algorithm
    alg->ProduceContainers();
    if(ND::tpcDebug().PatternRecognition(DB_INFO)) std::cout << " --- done pattern recognition for sub-event " << subEvent << " --- " << std::endl;
  };

  // set up container for hitpad level unused
  ND::THitSelection* usedTREx = new ND::THitSelection();

  // get patterns
  subEvent = 0;
  for(std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator algIt = fSubAlgorithms.begin(); algIt != fSubAlgorithms.end(); ++algIt){
    ND::TTPCTRExPatSubAlgorithm* alg = *algIt;
    subEvent++;
    // produce pattern for each sub-algorithm
    if(ND::tpcDebug().PatternRecognition(DB_INFO)) std::cout << " ----- pattern containers for sub-event " << subEvent << " ------ " << std::endl;
    alg->ProducePattern(usedTREx);
    if(ND::tpcDebug().PatternRecognition(DB_INFO)) std::cout << " --- done pattern containers for sub-event " << subEvent << " --- " << std::endl;
  };

  // fill unused hits
  FillUsedUnusedHits(usedTREx, used, unused);

  // clean up
  delete usedTREx;

  // print debug info on pattern constituents
  //VerifyHitsPattern();
  //VerifyUsedUnusedPattern(used, unused);
}

void ND::TTPCTRExPatAlgorithm::VerifyHitsPattern(){
  if(!ND::tpcDebug().PatternRecognition(DB_INFO)) return;

  std::cout << " ----- info on event's " << fSubAlgorithms.size() << " pattern containers ------" << std::endl;
  int subEvent=0;
  for(std::vector<ND::TTPCTRExPatSubAlgorithm*>::iterator algIt = fSubAlgorithms.begin(); algIt != fSubAlgorithms.end(); ++algIt){
    ND::TTPCTRExPatSubAlgorithm* alg = *algIt;
    subEvent++;

    ND::THandle<ND::TTPCPattern> pattern = alg->GetPattern();

    if(!pattern){
      std::cout << "\n  (Pattern) Sub event " << subEvent << " has no pattern (its " << alg->GetHitMapSize() << " hit elements may not be sufficient)" << std::endl;
    }
    else{
      std::cout << "\n  (Pattern) Pattern from sub event " << subEvent << " contains " << pattern->GetNJunctions() << " junctions and " << pattern->GetNPaths() << " paths" << std::endl;

      VerifyPatternConstituents(pattern);
    };
  };
  std::cout << "\n --- done info on event's " << fSubAlgorithms.size() << " pattern containers ---" << std::endl;
}
void ND::TTPCTRExPatAlgorithm::VerifyUsedUnusedPattern(ND::THitSelection* used, ND::THitSelection* unused){
  int nTotal = fHits->size();
  int nUsed = used->size();
  int nUnused = unused->size();

  if(ND::tpcDebug().PatternRecognition(DB_INFO)) std::cout << "\n -- total input hits -> " << nTotal << " -- " << std::endl;
  if(ND::tpcDebug().PatternRecognition(DB_INFO)) std::cout << " -- total (used + unused) hits -> (" << nUsed << " + " << nUnused << ") = " << (nUsed+nUnused) << " -- " << std::endl;

  if(nUsed+nUnused != nTotal){
    //TODO: throw proper exception
    std::cout << "FATAL ERROR:  hits missing from both used and unused" << std::endl;
    throw;
  };
}

void ND::TTPCTRExPatAlgorithm::VerifyHitsFull(ND::TReconObjectContainer* patterns){
  bool info = ND::tpcDebug().PatternRecognition(DB_INFO);
  bool fatalError = false;

  int patternID=0;
  unsigned int nPatterns = patterns->size();
  if(info) std::cout << " ----- info on event's " << nPatterns << " patterns ------" << std::endl;

  for(ND::TReconObjectContainer::iterator patternIt = patterns->begin(); patternIt != patterns->end(); ++patternIt){
    ND::THandle<ND::TTPCPattern> pattern = *patternIt;
    patternID++;

    if(!pattern){
      std::cout << "FATAL ERROR: pattern " << patternID << "/" << nPatterns << " not found when filling pattern recognition output" << std::endl;
      fatalError = true;
    }
    else{
      unsigned int nPaths = pattern->GetNPaths();
      unsigned int nJunctions = pattern->GetNJunctions();

      if(info) std::cout << "(Pattern) Pattern " << patternID << "/" << nPatterns << " contains " << nPaths << " paths and " << nJunctions << " junctions"<< std::endl;

      VerifyPatternConstituents(pattern);
    };
  };

  if(fatalError){
    //TODO: proper errors
    std::cout << "FATAL ERROR: fatal error found when verifying pattern recognition output; exiting" << std::endl;
    throw;
  };

  if(info) std::cout << "\n --- done info on event's " << fSubAlgorithms.size() << " pattern containers ---" << std::endl;
}
void ND::TTPCTRExPatAlgorithm::VerifyUsedUnusedFull(ND::TReconObjectContainer* patterns){
  bool info = ND::tpcDebug().PatternRecognition(DB_INFO);
  bool error = ND::tpcDebug().PatternRecognition(DB_ERROR);

  int nPathSharedHits = 0;
  int nPathUniqueHits = 0;
  int nJunctionHits = 0;
  int nUniqueHits = 0;

  // check input hits from paths and junctions separately, then together
  if(error){
    std::set<ND::THit*> horizontalPathHits;
    std::set<ND::THit*> verticalPathHits;
    std::set<ND::THit*> pathHits;
    std::set<ND::THit*> junctionHits;
    for(ND::TReconObjectContainer::iterator patternIt = patterns->begin(); patternIt != patterns->end(); ++patternIt){
      ND::THandle<ND::TTPCPattern> pattern = *patternIt;

      if(pattern){
        std::vector< ND::THandle<ND::TTPCPath> > paths = pattern->GetPaths();
        std::vector< ND::THandle<ND::TTPCJunction> > junctions = pattern->GetJunctions();

        for(std::vector< ND::THandle<ND::TTPCPath> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
          ND::THandle<ND::TTPCPath> path = *pathIt;
          ND::THandle<ND::THitSelection> pathHitSel = path->GetHits();

          std::set<ND::THit*> horizontalThisPathHits;
          std::set<ND::THit*> verticalThisPathHits;

          for(ND::THitSelection::iterator hvClusterIt = pathHitSel->begin(); hvClusterIt != pathHitSel->end(); ++hvClusterIt){
            ND::THandle<ND::TTPCHVCluster> hvCluster = *hvClusterIt;
            std::set<ND::THit*> thisClusterHits;

            if(hvCluster){
              const ND::THitSelection clusterHits = hvCluster->GetHits();
              for(ND::THitSelection::const_iterator hitPadIt = clusterHits.begin(); hitPadIt != clusterHits.end(); ++hitPadIt){
                ND::THit* hitPadAddr = &(**hitPadIt);

                if(thisClusterHits.count(hitPadAddr)){
                  if(error) std::cout << "ERROR:  TTPCHitPad duplicated in same HV cluster" << std::endl;
                }
                else{
                  thisClusterHits.insert(hitPadAddr);
                };
              };
            };
            for(std::set<ND::THit*>::iterator hitIt = thisClusterHits.begin(); hitIt != thisClusterHits.end(); ++hitIt){
              ND::THit* hit = *hitIt;

              if(hvCluster->IsVertical()){
                if(verticalThisPathHits.count(hit)){
                  std::cout << "ERROR:  TTPCHitPad duplicated between two vertical clusters in the same path" << std::endl;
                }
                else{
                  verticalThisPathHits.insert(hit);
                };
              }
              else if(hvCluster->IsHorizontal()){
                if(horizontalThisPathHits.count(hit)){
                  std::cout << "ERROR:  TTPCHitPad duplicated between two horizontal clusters in the same path" << std::endl;
                }
                else{
                  horizontalThisPathHits.insert(hit);
                };
              }
              else{
                std::cout << "ERROR:  Cluster is neither horizontal or vertical; can't check its hits" << std::endl;
              };
            };
          };

          // now check this path's hits against all paths
          for(std::set<ND::THit*>::iterator hitIt = horizontalThisPathHits.begin(); hitIt != horizontalThisPathHits.end(); ++hitIt){
            ND::THit* hit = *hitIt;

            if(horizontalPathHits.count(hit)){
              std::cout << "ERROR:  TTPCHitPad is duplicated between horizontal clusters in different paths" << std::endl;
            }
            else{
              horizontalPathHits.insert(hit);
              pathHits.insert(hit);
              nPathUniqueHits++;
            };
          };
          for(std::set<ND::THit*>::iterator hitIt = verticalThisPathHits.begin(); hitIt != verticalThisPathHits.end(); ++hitIt){
            ND::THit* hit = *hitIt;

            if(verticalPathHits.count(hit)){
              std::cout << "ERROR:  TTPCHitPad is duplicated between vertical clusters in different paths" << std::endl;
            }
            else{
              verticalPathHits.insert(hit);
              if(pathHits.count(hit)){
                nPathUniqueHits--;
                nPathSharedHits++;
              }
              else{
                pathHits.insert(hit);
                nPathUniqueHits++;
              };
            };
          };
        };

        for(std::vector< ND::THandle<ND::TTPCJunction> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
          ND::THandle<ND::TTPCJunction> junction = *junctionIt;
          ND::THandle<ND::THitSelection> junctionHitSel = junction->GetHits();

          std::set<ND::THit*> thisJunctionHits;
          for(ND::THitSelection::iterator junctionHitIt = junctionHitSel->begin(); junctionHitIt != junctionHitSel->end(); ++junctionHitIt){
            ND::THit* hitPadAddr = &(**junctionHitIt);

            if(thisJunctionHits.count(hitPadAddr)){
              std::cout << "ERROR:  TTPCHitPad is duplicated in the same junction" << std::endl;
            }
            else{
              thisJunctionHits.insert(hitPadAddr);
            };
          };

          for(std::set<ND::THit*>::iterator hitIt = thisJunctionHits.begin(); hitIt != thisJunctionHits.end(); ++hitIt){
            ND::THit* hit = *hitIt;

            if(junctionHits.count(hit)){
              std::cout << "ERROR:  TTPCHitPad is duplicated between different junctions" << std::endl;
            }
            else{
              junctionHits.insert(hit);
              nJunctionHits++;
            };
          };
        };
      };
    };

    nUniqueHits += junctionHits.size();
    for(std::set<ND::THit*>::iterator hitIt = pathHits.begin(); hitIt != pathHits.end(); ++hitIt){
      ND::THit* hit = *hitIt;

      if(junctionHits.count(hit)){
        std::cout << "ERROR:  TTPCHitPad is duplicated between a path and a junction" << std::endl;
      }
      else{
        nUniqueHits ++;
      };
    };

    // print summary
    if(info) std::cout << "\n - total unique hits: " << nUniqueHits << std::endl;
    if(info) std::cout <<   " |-- cluster shared hits: " << nPathSharedHits << ", cluster unique hits: " << nPathUniqueHits << std::endl;
    if(info) std::cout <<   " |-- junction hits: " << nJunctionHits << std::endl;
  };
}
void ND::TTPCTRExPatAlgorithm::VerifyPatternConstituents(ND::THandle<ND::TTPCPattern> pattern){
  bool verbose = ND::tpcDebug().PatternRecognition(DB_VERBOSE);

  if(verbose){
    std::vector< ND::THandle<ND::TTPCJunction> > junctions = pattern->GetJunctions();
    if(junctions.size() > 0) std::cout << "\n  (Pattern) Analysing " << pattern->GetNJunctions() << " junctions" << std::endl;
    else std::cout << "\n  (Pattern) Junctions empty" << std::endl;
    for(std::vector< ND::THandle<ND::TTPCJunction> >::iterator junctionIt = junctions.begin(); junctionIt != junctions.end(); ++junctionIt){
      ND::THandle<ND::TTPCJunction> junction = *junctionIt;

      std::cout << "    (Junction) Junction id " << junction->GetId() <<  "; associated with " << junction->GetNPaths() << " paths" << std::endl;

      ND::THandle<ND::THitSelection> junctionHits = junction->GetHits();
      std::cout << "      (Junction hits) Junction contains " << junctionHits->size() << " hits" << std::endl;
      bool nullJunctionHits = false;
      for(ND::THitSelection::iterator junctionHitIt = junctionHits->begin(); junctionHitIt != junctionHits->end(); ++junctionHitIt){
        ND::THandle<ND::TTPCHitPad> junctionHit = *junctionHitIt;
        if(!junctionHit){
          nullJunctionHits = true;
          break;
        };
      };
      if(nullJunctionHits) std::cout << "      (Junction hits) Junction contains NULL hit pads!" << std::endl;

      for(ND::TReconObjectContainer::iterator pathIt = junction->GetConstituents()->begin(); pathIt != junction->GetConstituents()->end(); ++pathIt){
        ND::THandle<ND::TTPCPath> path = *pathIt;
        if(path){
          std::cout << "      (Path) Path id " << path->GetId() << std::endl;
        }
        else{
          std::cout << "      (Path) NULL PATH!" << std::endl;
        };
      };
    };

    std::vector< ND::THandle<ND::TTPCPath> > paths = pattern->GetPaths();
    if(paths.size() > 0) std::cout << "\n  (Pattern) Analysing " << pattern->GetNPaths() << " paths" << std::endl;
    else std::cout << "\n  (Pattern) Paths empty" << std::endl;
    for(std::vector< ND::THandle<ND::TTPCPath> >::iterator pathIt = paths.begin(); pathIt != paths.end(); ++pathIt){
      ND::THandle<ND::TTPCPath> path = *pathIt;
      ND::THandle<ND::THitSelection> pathHits = path->GetHits();

      std::cout << "    (Path) Path id " << path->GetId() <<  "; associated with " << pathHits->size() << " HV clusters" << std::endl;
      std::cout << "      (HV Cluster) HV cluster sizes: ";
      for(ND::THitSelection::iterator hvClusterIt = pathHits->begin(); hvClusterIt != pathHits->end(); ++hvClusterIt){
        ND::THandle<ND::TTPCHVCluster> hvCluster = *hvClusterIt;
        if(hvCluster){
          const ND::THitSelection clusterHits = hvCluster->GetHits();
          bool allValid = true;

          for(ND::THitSelection::const_iterator hitPadIt = clusterHits.begin(); hitPadIt != clusterHits.end(); ++hitPadIt){
            ND::THandle<ND::TTPCHitPad> hitPad = *hitPadIt;
            if(!hitPad){
              allValid = false;
              break;
            };
          };

          if(allValid){
            std::cout << " " << clusterHits.size();
          }
          else{
            std::cout << " (NULL HIT PADS IN CLUSTER)";
          };
        }
        else{
          std::cout << " (NULL HV CLUSTER)";
        };
      };
      std::cout << "\n";
      std::cout << "      (HV Cluster) HV cluster orientations: ";
      for(ND::THitSelection::iterator hvClusterIt = pathHits->begin(); hvClusterIt != pathHits->end(); ++hvClusterIt){
        ND::THandle<ND::TTPCHVCluster> hvCluster = *hvClusterIt;
        if(hvCluster) std::cout << " " << ( (hvCluster->IsVertical()) ? "v" : "h" );
      };
      std::cout << "\n" << std::endl;
    };
  };
}
