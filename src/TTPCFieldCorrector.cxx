#include "TTPCFieldCorrector.hxx"

#include <TFieldManager.hxx>
#include <TElectronCloud.hxx>
#include <TOARuntimeParameters.hxx>
#include <TGeomInfo.hxx>
#include <HEPUnits.hxx>
#include <TGeoManager.h>

#include <iostream>

using namespace std;

TTPCFieldCorrector::TTPCFieldCorrector(bool isMC) {
  // Get a reference to the singleton instance of TPC geometry information
  ND::TTPCGeom& tpcGeom = ND::TGeomInfo::TPC();
  // Get the maximum drift x coordinate
  fMaxDriftGlobalX = tpcGeom.GetMaxDriftDistance() + tpcGeom.GetCathodeWidth() / 2.;
  fCathode =  tpcGeom.GetCathodeWidth() / 2.;

  //Before using TElectronCloud one have to initilize the magnetic field manager 
  //which is a singleton. Default it uses the magnetic field type requested in DB
  ND::TFieldManager::InitializeFieldManager();

  fDriftVelocity = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCgas.DriftSpeed");
  fUseAdaptiveStepSize = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.TPCDrift.UseAdaptiveStepSize");

  // initialize the drift manager
  double InitialDriftStepTime = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.InitialDriftStepTime");
  bool UseAdaptiveStepSize    = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.TPCDrift.UseAdaptiveStepSize");
  double OmegaTau             = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.OmegaTau");
  double Mu                   = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.ElectronMobility");
  
  if (UseAdaptiveStepSize) {
    std::cout << "USE ADAPTIVE STEP " << std::endl;
    double MaxDriftStepTime = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.MaxDriftStepTime");
    double MinDriftStepTime = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.MinDriftStepTime");
    int MaxNbrDriftSteps    = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.TPCDrift.MaxNbrDriftSteps");
    double ErrorTol         = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.ErrorTol");

    fDriftManager.InitializeFromDB(UseAdaptiveStepSize, MaxDriftStepTime, MinDriftStepTime, MaxNbrDriftSteps, ErrorTol, InitialDriftStepTime, OmegaTau, Mu, fDriftVelocity);
  } else {
    fDriftManager.InitializeFromDB(InitialDriftStepTime, OmegaTau, Mu, fDriftVelocity);
  }

  unsigned int EFDistModel = ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.TPCDrift.EFieldDistModel");

  double EFDistStrength = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.EFieldDistStrength");
  double EFDistYFactor = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.EFieldDistYFactor");
  double EFDistZFactor = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.TPCDrift.EFieldDistZFactor");
  
  //values obtained from the Alf2 fit model
  double EFDist_A[3][2];
  double EFDist_B[3][2];
  double EFDist_c[3][2];
  double EFDist_d[3][2];
  double EFDist_m[3][2];
  double EFDist_n[3][2];

  char parameter[250];  
  for ( int tpc=0; tpc<3; tpc++) {
    for ( int rp=0; rp<2; rp++) {
      sprintf(parameter,"trexRecon.TPCDrift.EFieldDistA_TPC%dRP%d",tpc+1,rp);
      EFDist_A[tpc][rp] = ND::TOARuntimeParameters::Get().GetParameterD(parameter);
      sprintf(parameter,"trexRecon.TPCDrift.EFieldDistB_TPC%dRP%d",tpc+1,rp);
      EFDist_B[tpc][rp] = ND::TOARuntimeParameters::Get().GetParameterD(parameter);
      sprintf(parameter,"trexRecon.TPCDrift.EFieldDistc_TPC%dRP%d",tpc+1,rp);
      EFDist_c[tpc][rp] = ND::TOARuntimeParameters::Get().GetParameterD(parameter);
      sprintf(parameter,"trexRecon.TPCDrift.EFieldDistd_TPC%dRP%d",tpc+1,rp);
      EFDist_d[tpc][rp] = ND::TOARuntimeParameters::Get().GetParameterD(parameter);
      sprintf(parameter,"trexRecon.TPCDrift.EFieldDistm_TPC%dRP%d",tpc+1,rp);
      EFDist_m[tpc][rp] = ND::TOARuntimeParameters::Get().GetParameterD(parameter);
      sprintf(parameter,"trexRecon.TPCDrift.EFieldDistn_TPC%dRP%d",tpc+1,rp);
      EFDist_n[tpc][rp] = ND::TOARuntimeParameters::Get().GetParameterD(parameter);
    }
  }
  

  // This is in simulate mode for MC, i.e. this will add the effect
  // of the EField distortion on a perfect MC.
  // Of course it will correct on data.

  if (isMC) {   
    EFDistStrength = -1. * EFDistStrength;    
    for (unsigned int tpc = 0; tpc < 3; tpc++){
      for (unsigned int rp = 0; rp < 2; rp++){      
        EFDist_A[tpc][rp] = -1. * EFDist_A[tpc][rp];
        EFDist_B[tpc][rp] = -1. * EFDist_B[tpc][rp]; 
      }
    }
  }

  // Use standard Scott1 field from oaUtility/TDriftManager.
  // Define the EField for each drift volume independently.
  if (EFDistModel){
    for (unsigned int tpc = 0; tpc < 3; tpc++){
      for (unsigned int rp = 0; rp < 2; rp++){
        if (EFDistModel == 1)
          fDriftManager.SetEFieldModel(tpc, rp, (new TEFieldScott1( EFDistStrength, EFDistYFactor, EFDistZFactor)));
        if (EFDistModel == 2) {
          fDriftManager.SetEFieldModel( tpc, rp, ( new TEFieldAlf2( EFDist_A[tpc][rp], EFDist_B[tpc][rp], EFDist_c[tpc][rp], EFDist_d[tpc][rp], EFDist_m[tpc][rp], EFDist_n[tpc][rp] ) ) );
        }
      }
    }
  }
}

TTPCFieldCorrector::~TTPCFieldCorrector() {
}

bool TTPCFieldCorrector::GetFieldCorrectedPoint(const TVector3& recPos, TVector3 &origPos) {
  //The x position must have been achieved using 
  //assumed drift velocity given in the data base.
  //Then we can get the drift time
  double dt = fabs(fabs(recPos.X()) - fMaxDriftGlobalX) / fDriftVelocity;

  //Since the fields were assumed to be homogenetic the z and y positions 
  //are the same (have to ignore transverse diffusion)
  //Create a cloud (# electrons are irrelevant) at the hit position
  //with oposite charge since we want to drift it back to its creation point
  double padPlaneX = fMaxDriftGlobalX; //Cathode at x = 0 globally
  if (recPos.X() < 0) padPlaneX = -padPlaneX;
  //Create the "positron" cloud at the pad plane point
  ND::TElectronCloud eCloud(1, TVector3(padPlaneX, recPos.Y(), recPos.Z()), 
      TVector3(0, 0, 0), 0, +1);
  //Now drift the cloud from the readout board back to the correct original point
  if (!fDriftManager.DriftUntilTLim(eCloud,dt)) return false;
  origPos = eCloud.GetPosition();
  return true;
}

bool TTPCFieldCorrector::GetFieldCorrectedPoint(const TVector3& padPlanePoint, double dt, TVector3 &origPos) {
  //Since the fields were assumed to be homogenetic the z and y positions 
  //are the same (have to ignore transverse diffusion)
  //Create a cloud (# electrons are irrelevant) at the pad plane point
  //with oposite charge since we want to drift it back to its creation point
  ////////
  //Just a test to make sure the input is really a pad plane point
  if (fabs(fabs(padPlanePoint.X()) - fMaxDriftGlobalX) > 0.0001) {
    cout<<"ERROR in TTPCFieldCorrector::GetFieldCorrectedPoint; input was not a pad plane point! Exit!"<<endl;
    exit(1);
  }
  //Create the "positron" cloud at the pad plane point
  ND::TElectronCloud eCloud(1, padPlanePoint, TVector3(0, 0, 0), 0, +1);
  //Now drift the cloud from the readout board back to the correct original point
  if (!fDriftManager.DriftUntilTLim(eCloud,dt)) return false;
  origPos = eCloud.GetPosition();
  return true;
}

bool TTPCFieldCorrector::GetHitOnPadPlane(const TLorentzVector& trackPoint, TLorentzVector &padPoint) {
  //Create a cloud (# electrons are irrelevant) at the track point
  ND::TElectronCloud eCloud(1, trackPoint.Vect(), TVector3(0, 0, 0), trackPoint.T());
  //Now drift the cloud from the track point to the readout board
  //Ensure that the drift is done to the correct readout plane                                                            
  double padPlaneX = fMaxDriftGlobalX; //Cathode at x = 0 globally
  if (trackPoint.X() < 0) padPlaneX = -padPlaneX;

  if (!fDriftManager.DriftUntilXLim(eCloud,padPlaneX)) return false;
  TVector3 padPos = eCloud.GetPosition();
  double arrivalTime = eCloud.GetTime();
  padPoint = TLorentzVector(padPos, arrivalTime);
  return true;
}

bool TTPCFieldCorrector::GetHitOnCathode(const TVector3& padPos, TVector3 &catPos) {

  //CORRECTION WITH REVERSE DRIFT (will take a point in pad plane and drift it back to the cathode)-> For Laser
 
  //Create a cloud (# electrons are irrelevant) at the hit position
  //with oposite charge since we want to drift it back to the cathode
 
  //Create the "positron" cloud at the pad plane point
  ND::TElectronCloud eCloud(1, TVector3(padPos.X(), padPos.Y(), padPos.Z()), 
      TVector3(0, 0, 0), 0, +1);
  
  if (!fDriftManager.DriftUntilXLim(eCloud,fCathode)) return false;
  catPos = eCloud.GetPosition();
  return true;
}
