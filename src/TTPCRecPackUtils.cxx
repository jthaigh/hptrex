#include "TTPCRecPackUtils.hxx"
#include <recpack/KalmanFitter.h>
#include <TRecPackManager.hxx>
#include <TOARuntimeParameters.hxx>


namespace TTPCRecPackUtils {
  //*****************************************************************************
  void InitRecPackManager(void) {

    // add a new RecpackManager 
    ND::tman().add_manager("TREx");
    ND::tman().select_manager("TREx");
    
    // only surfaces with measurements are tried inside a volume
    ND::rpman("TREx").navigation_svc().navigator(RP::particle_helix).set_unique_surface(true);
    
    // select the Kalman Filter as fitting algorithm
    ND::rpman("TREx").fitting_svc().select_fitter(RP::kalman);
    
    // set a big local chi2 cut for global reconstruction
    ND::rpman("TREx").fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,RP::particle_helix).set_max_local_chi2ndf(EGeo::inf());
    
    // enable energy loss correction
    ND::rpman("TREx").model_svc().enable_correction(RP::particle_helix, RP::eloss, true);
    
    // set the matching representation
    ND::rpman("TREx").matching_svc().set_matching_representation(RP::slopes_curv);
    
    //set the minimum momentum for RecPack matching in the same TPC
    double MinimumMomentum = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MinimumMomentum");
    ND::tman("TREx").SetMinimumMomentum(MinimumMomentum);
    
    // set the maximum matching residual cut for this manager
    EVector maxResidual(6,0);
    
    maxResidual[0]     = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MaximumResidualX");
    maxResidual[1]     = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MaximumResidualY");
    maxResidual[2]     = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MaximumResidualZ");
    maxResidual[3]     = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MaximumResidualTX");
    maxResidual[4]     = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MaximumResidualTY");
    maxResidual[5]     = ND::TOARuntimeParameters::Get().GetParameterD("trexRecon.RecPack.MaximumResidualQoverP");
    
    ND::tman("TREx").SetMaxResidualCut(maxResidual);
    
    
    std::string RECPACKVerbosity = ND::TOARuntimeParameters::Get().GetParameterS("trexRecon.RecPack.Verbosity");
    // Changes the RecPack verbosity (5 digits = fitting, navigation, model, matching, RayTool. Minimum is 0, maximum is 7)
    ND::tman("TREx").SetVerbosity(RECPACKVerbosity);
    
    // TODO the RecPackConverters verbosity
    //ND::converter().SetVerbosity(fConverterVerb);
    
    if(RECPACKVerbosity.compare("0000") == 0){
      std::string muted("MUTE");
      Messenger::Level mute = Messenger::str(muted);
      ND::rpman("TREx").navigation_svc().navigator(RP::particle_helix).set_verbosity(mute);
      ND::rpman("TREx").model_svc().model(RP::particle_helix).equation().set_verbosity(mute);
    }
  }


  //*****************************************************************************
  SavedParam InitForPropagationInTPC(){
    SavedParam Params;
    // Prepare RecPack for the propagation
    // Save values that we are going to modify to not mess up the rest of the reconstruction.
    Params.eloss_on       = ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).correction(RP::eloss).status("enabled");
    Params.eloss_fluct_on = ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).noiser(RP::eloss).status("enabled");
    Params.ms_on          = ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).noiser(RP::ms).status("enabled");
    Params.length_sign    = ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).intersector().length_sign();
    Params.allow_zero     = ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).intersector().allow_zero_length();

    Params.unique_surface     = ND::rpman("TREx").navigation_svc().navigator(ND::rpman("TREx").model_svc().model_name()).unique_surface();

    // Disable energy loss fluctuations and multiple scattering
    // while creating the output node states
    ND::rpman("TREx").model_svc().enable_correction(RP::particle_helix, RP::eloss, false);
    ND::rpman("TREx").model_svc().enable_noiser(RP::particle_helix, RP::eloss, false);
    ND::rpman("TREx").model_svc().enable_noiser(RP::particle_helix, RP::ms, false);
    // Important to prevent RP from misbehaving
    ND::rpman("TREx").model_svc().model().intersector().set_length_sign(1);
    ND::rpman("TREx").model_svc().model().intersector().set_allow_zero_length(true);

    ND::rpman("TREx").navigation_svc().navigator(ND::rpman("TREx").model_svc().model_name()).set_unique_surface(true);

    return Params;
  }

  //*****************************************************************************
  bool PropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState){
    double length =0.0;
    return PropagateToHVCluster(Cluster, PropagState, length);
  }

  //*****************************************************************************
  bool PropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState, double &Length){
    State newState = PropagState;

    // propagate the track to that surface
    bool ok = ND::rpman().navigation_svc().propagate_vector(Cluster->GetPropagSurf(), newState, Length);

    // This fails. Keep the input state intact. Maybe the next propagation will be more
    // successful
    if ( !ok)
      return false;

    // Always propagate the state from the previous node in order to reduce
    // risks inherent to long propagations.
    PropagState = newState;
    return true;
  }

  //*****************************************************************************
  bool FullPropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState){
    double length =0.0;
    return FullPropagateToHVCluster(Cluster, PropagState, length);
  }

  //*****************************************************************************
  bool FullPropagateToHVCluster(ND::THandle<ND::TTPCHVCluster> Cluster, State &PropagState, double &Length){
    // The propagation surface is defined by the cluster orientation.
    // Maybe for large X angle we should switch to YZ planes ...
    TVector3 projNorm(0.,0.,0.);
    if ( Cluster->IsVertical() ){
      projNorm.SetZ(1.);
    } else {
      projNorm.SetY(1.);
    }
    State newState;
    TVector3 cluPos = Cluster->GetCalibratedPosition();
    bool RPsucceeded = ND::tman().PropagateToSurface(PropagState,cluPos,projNorm,newState, Length);
    // This fails. Keep the input state intact. Maybe the next propagation will be more
    // successful
    if ( !RPsucceeded)
      return false;

    // Always propagate the state from the previous node in order to reduce
    // risks inherent to long propagations.
    PropagState = newState;
    return true;
  }

  //*****************************************************************************
  void ResetAfterPropagationInTPC(SavedParam &Params){
    // Reenable energy loss fluctuations and multiple scattering previously turned off
    ND::rpman("TREx").model_svc().enable_correction(ND::rpman("TREx").model_svc().model_name(),  RP::eloss, Params.eloss_on);
    ND::rpman("TREx").model_svc().enable_noiser(ND::rpman("TREx").model_svc().model_name(),  RP::eloss, Params.eloss_fluct_on);
    ND::rpman("TREx").model_svc().enable_noiser(ND::rpman("TREx").model_svc().model_name(),  RP::ms, Params.ms_on);

    ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).intersector().set_length_sign(Params.length_sign);
    ND::rpman("TREx").model_svc().model(ND::rpman("TREx").model_svc().model_name()).intersector().set_allow_zero_length(Params.allow_zero);

    ND::rpman("TREx").navigation_svc().navigator(ND::rpman("TREx").model_svc().model_name()).set_unique_surface(Params.unique_surface);
  }
}
