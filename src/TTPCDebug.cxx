#include "TTPCDebug.hxx"
#include <TOARuntimeParameters.hxx>


/// The static member pointer to the singleton.
ND::TTPCDebug* ND::TTPCDebug::_tpcDebug = NULL;

//*****************************************************************************
ND::TTPCDebug& ND::tpcDebug(){

  return ND::TTPCDebug::Get();

}

//*****************************************************************************
ND::TTPCDebug& ND::TTPCDebug::Get(void) {
  if (!_tpcDebug) {
    _tpcDebug = new ND::TTPCDebug();
  }

  return *_tpcDebug;
}


//*****************************************************************************
ND::TTPCDebug::TTPCDebug(){
  fGeneralSteps = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.GeneralSteps.Debug");

  fPatternRecognition = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.PatternRecognition.Debug");

  fSeeding = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.Seeding.Debug");
  fSeededClusters = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.SeededClusters.Debug");
  fMMHoriGapMerge = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.MMHoriGapMerge.Debug");
  fMMVertGapMerge = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.MMVertGapMerge.Debug");
  fT0 = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.T0.Debug");
  fLikFit = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFit.Debug");
  fLikFittedClusters = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikFittedClusters.Debug");
  fHelixPropagator = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.HelixPropagator.Debug");
  fLikelihoodMatch = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikelihoodMatch.Debug");
  fLikelihoodMerge = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.LikelihoodMerge.Debug");
  fCathCrosserMerge = (DBLEVELS) ND::TOARuntimeParameters::Get().GetParameterI("trexRecon.Reco.CathCrosserMerge.Debug");
}


//*****************************************************************************
bool ND::TTPCDebug::GeneralSteps(DBLEVELS level){
  if( fGeneralSteps >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::PatternRecognition(DBLEVELS level){
  if( fPatternRecognition >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::Seeding(DBLEVELS level){
  if( fSeeding >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::SeededClusters(DBLEVELS level){
  if( fSeededClusters >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::MatchAndMerge(DBLEVELS level){
  if( fMMVertGapMerge >= level || fMMHoriGapMerge >= level) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::MMHoriGapMerge(DBLEVELS level){
  if( fMMHoriGapMerge >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::MMVertGapMerge(DBLEVELS level){
  if( fMMVertGapMerge >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::T0Finder(DBLEVELS level){
  if( fT0 >= level ) 
    return true;
  return false;
}


//*****************************************************************************
// Turn on track printouts if LikFit or other tracking algorithm printouts are on
bool ND::TTPCDebug::Tracking(DBLEVELS level){
  if( fLikFit >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::LikFit(DBLEVELS level){
  if( fLikFit >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::LikFittedClusters(DBLEVELS level){
  if( fLikFittedClusters >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::HelixPropagator(DBLEVELS level){
  if( fHelixPropagator >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::LikelihoodMatch(DBLEVELS level){
  if( fLikelihoodMatch >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::LikelihoodMerge(DBLEVELS level){
  if( fLikelihoodMerge >= level ) 
    return true;
  return false;
}


//*****************************************************************************
bool ND::TTPCDebug::CathCrosserMerge(DBLEVELS level){
  if( fCathCrosserMerge >= level ) 
    return true;
  return false;
}


