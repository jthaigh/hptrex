#ifndef TTPCDebug_hxx_seen
#define TTPCDebug_hxx_seen

/// Definition of the debug level.
typedef enum {DB_QUIET=0,DB_ERROR,DB_INFO,DB_VERBOSE,DB_VVERBOSE} DBLEVELS; 

namespace ND {
  class TTPCDebug;
}

class ND::TTPCDebug {
public:

  virtual ~TTPCDebug() {};

  /// Get a pointer to the singleton instance of the calibration information.
  static TTPCDebug& Get(void);

  /// Return true if wanted level is higher or equal to the current
  /// debug level for general steps in TRExReco.
  bool GeneralSteps(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for pattern recognition.
  bool PatternRecognition(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for Seeding.
  bool Seeding(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for the display of the clusters used in the seeding.
  bool SeededClusters(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for MMVertGapMerge or MMHoriGapMerge.
  bool MatchAndMerge(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for MMHoriGapMerge.
  bool MMHoriGapMerge(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for MMVertGapMerge.
  bool MMVertGapMerge(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for T0Finder.
  bool T0Finder(DBLEVELS level);
  /// Return true if the debug level of the LikFit or another
  /// tracking algorithm requires it.
  bool Tracking(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for LikFit.
  bool LikFit(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for the display of the clusters used in the
  /// LogLikelihood minimizer or calculator.
  bool LikFittedClusters(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for HelixPropagator.
  bool HelixPropagator(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for LikelihoodMatch.
  bool LikelihoodMatch(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for LikelihoodMerge.
  bool LikelihoodMerge(DBLEVELS level);
  /// Return true if wanted level is higher or equal to the current
  /// debug level for CathCrosserMerge.
  bool CathCrosserMerge(DBLEVELS level);

private:
  TTPCDebug();

  static TTPCDebug* _tpcDebug;

  /// Flag for General steps
  DBLEVELS fGeneralSteps;
  /// Flag for pattern recognition
  DBLEVELS fPatternRecognition;
  /// Flag for SeedingPath
  DBLEVELS fSeeding;
  /// Flag to display the clusters used in the seeding
  DBLEVELS fSeededClusters;
  /// Flag for MMHoriGapMerge
  DBLEVELS fMMHoriGapMerge;
  /// Flag for MMVertGapMerge
  DBLEVELS fMMVertGapMerge;
  /// Flag for T0 finder
  DBLEVELS fT0;
  /// Flag for LikFitPath
  DBLEVELS fLikFit;
  /// Flag to display the clusters used in the LogLikelihood minimizer or calculator
  DBLEVELS fLikFittedClusters;
  /// Flag for HelixPropagator
  DBLEVELS fHelixPropagator;
  /// Flag for LikelihoodMatch
  DBLEVELS fLikelihoodMatch;
  /// Flag for LikelihoodMerge
  DBLEVELS fLikelihoodMerge;
  /// Flag for CathCrosserMerge
  DBLEVELS fCathCrosserMerge;

};


namespace ND{
  TTPCDebug& tpcDebug();
};

#endif

