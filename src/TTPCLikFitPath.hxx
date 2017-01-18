#ifndef TTPCLikFitPath_hxx_seen
#define TTPCLikFitPath_hxx_seen

#include <TMinuit.h>
#include "TTPCQLikelihood.hxx"
#include "TTRExPath.hxx"
#include "TTPCUtils.hxx"
#include "TTPCLayout.hxx"

#define NPARAM 9


/// We need the CLHEP package for the matrix calculations.
// using namespace CLHEP;

namespace trex {
  class TTPCLikFitPath;
}


/// The TPC likelihood Algorithm.

class trex::TTPCLikFitPath: public trex::TTPCQLikelihood{

  public:

    /// Default constructor.
    TTPCLikFitPath(bool CalculatorMode = false);

    /// Default destructor.
    virtual ~TTPCLikFitPath();

    /*!
      Select the wanted clusters for the likelihood fit.
      Save the selected ones in the output THitSelection in order to pass to the Process
      method only the clusters we actually want to fit.
      This is in order to optimize if possible the code but also give some
      flexibility to fit only a subset of the clusters by calling Process on a hand
      made selection of clusters, like for example by giving Process only one cluster
      at a time to do the Space Point Resolution study.
      Also store in the cluster using isOkForFit if we used it for fitting or not.
    */
  void PrepareClustersForFitting(std::vector<trex::TTRExHVCluster*>& inputClu, std::vector<trex::TTRExHVCluster*>& outputClu, double XDirection);

    /// Setup various parameters of the log likelihood minimization
    /// Use this to setup the minimization parameters, there initial values, the step sizes, etc...
    /// Also define the modes of the fit, which parameters to fit, etc ...
  void SetupLogLklhdMinimizer(trex::TTRExPath& Path, bool UseSeedAsInit = true);

    /// Minimize the log likelihood for the given hits.
    /// This is the place where parameters will be fixed or released but nothing else.
    /// Use SetupLogLklhdMinimizer to setup the parameters, there initial values, the step sizes, etc...

  int LogLklhdMinimizer(std::vector<trex::TTRExHVCluster*>& inputClusters);

    /// A few little calls to prepare for a minimization.
    /// Must be called explicitely before calling SimpleLogLklhdMinimizer
  void GetReadyForMinimization(std::vector<trex::TTRExHVCluster*>& inputClusters);

    /// Some clean of containers after the minimization.
    /// Must be called explicitely after calling SimpleLogLklhdMinimizer
    void CleanUpAfterMinimization();

    /// Call the printout command in Minuit for the details minimization results.
    void MinuitPrintout();

    /// Does only one minimization unlike LogLklhdMinimizer which will
    /// try to fix/release various parameters to get to the convergence.
    int SimpleLogLklhdMinimizer(bool ParamIsFree[NPARAM]);

    /// Setup various parameters of the log likelihood calculation
    /// Use this to setup the calculation of the log likelihood for a given set of clusters.
    /// Needs to know if the front of back likelihood state should be used for the propagation.
  bool SetupLogLklhdCalculator(std::vector<double> helixState, std::vector<trex::TTRExHVCluster*> inputClusters, double inputLength);

    /// Minimize the log likelihood for the hits given to SetupLogLklhdCalculator.
    /// This is the place where parameters will be fixed or released but nothing else.
    /// Use SetupLogLklhdMinimizer to setup the parameters, there initial values, the step sizes, etc...
    TTPCLogLikelihood LogLklhdCalculator();

    /// Reset all the class variables prior to the next fit.
    /// This should be called after you are done with a set of clusters
    /// to avoid memory leaks. Do not call it after PrepareClustersForFitting
    /// when you want to use these clusters.
    void Reset(void);

    /// Retrieves the value of the likelihood.
  double log_likelihood(std::vector<double>& x);

    void StoreFittedState();
    void StoreFittedSigma();

  void SaveFitResults(trex::TTRExPath& Path);

  TTPCPathFitResults GetFitResults();

  private: 
    /// Print the event information.
    void   printEvent(void);

    TMinuit* GetMinuitPtr() {return fMinuit;};

    /// Retrieve some defeult parameters. 
    void   DefaultFixedParameters(void);

    /// Reset all the minuit fitting parameters to their initial values from fInitValues.
    void ResetMinuitParam(void);

    /// Load initial values of the likelihood based on given state
  void StateToInitValues( std::vector<double> inputState);

    /// Defines the clusters or hit pads that will be matched to a path.
    struct ClusterSelection {
      public:
        int NMaxPeaks; 
        int NSelVert; 
        int NSelHori; 
        int NHoriMMEdge; 
        int NVertMMEdge; 
        int NOutOfChargeWindow; 
        int NOutDeltaDrift; 
        int NTooManyPadsPerClu; 
        int NSuspiciousPadTiming; 
        int NSaturation; 
    };

    /// Performs the selection of the clusters suitable for the fit
  void SelectClusters(std::vector<trex::TTRExHVCluster*>& inputClu, double XDirection, ClusterSelection &CluSel);

    /// Retrieves the value of the likelihood in Y along the vertical direction.
    double log_likelihoodHV();
    /// Retrieves the value of the likelihood in X along the vertical direction.
    double log_likelihoodX();

    /// Reweight factor for the multiple scattering effects on the track.
    double MScCorrection(double tx0, double curv); 

    bool fCalculatorMode;

    /// Store the initial values of the parameters.
    /// This will contain some values from the seed in particular.
    Double_t fInitValues[NPARAM];

    /// Steps used in minuit for each parameter.
    Double_t fStep[NPARAM];

    /// Fit the Y position or the Z position depending on the orientation
    /// of the first cluster.
    bool fFitYPosParam;

    /// Value to control the printlevel inside Minuit.
    int fPrintLevel;

    /// Calculated value of mean value of drift distance
    double fMeanDrift;



    /// Pointer to the TMinuit class used during the minimization.
    TMinuit *fMinuit;

    /// Error flag from Minuit.
    Int_t fIErrorFlag;

    /// Controls the print level for Minuit directly
    int fMinuitPrintLevel;
    /// Maximum number of iterations during Minuit minimization
    int fMinuitMaxIterations;

    /// BField at (0.0, 0.0, 0.0), originally used by the MSc correction
    double fBFieldAt0;
    /// Multiple scattering correction
    double fMScCorr;
    /// Switch to store the likelihood values (from minimization or calculation)
    /// with or with the multiple scattering correction applied.
    double fStoreLklhdWithMScCorr;

    /// Local copy of pad width.
    double fPadWidth;
    /// local copy of the pad height.
    double fPadHeight;


    /// Boolean to enable/disable the fitting of the transverse diffusion.
    bool fFitSigma;
    /// Boolean to fit sigma and the other parameters separately
    bool fFitSigmaSeparately;
    /// Boolean to enable/disable the fitting of the track curvature.
    bool fFitRho;
    /// Boolean to enable/disable the fitting of the track X parameters.
    bool fFitX; 
    /// Boolean to enable/disable the fitting of the track angle.
    bool fFitDir;

    /// Minimum charge allowed to select a cluster for the fit.
    double fMinimumCharge;
    /// Maximun charge allowed to select a cluster for the fit.
    double fMaximumCharge;  
    /// This is the maximum size in time*drift_speed to be use the row in the fit.  
    double fMaxDeltaDrift;  
    /// This is the minimal integral charge in a row(column) to be considered for the fit.
    double fMinPredIntCharge; 
    /// Noise value used during the likelihood fit.
    double fNoise;


    /// SigmaZ (for vertical cluster, SigmaY otherwise) value of the transfer diffusion sigmas.
    double fInitialSigma;

    /// Enable the combined XY fit. False means that the YZ and XZ fuits are done separately. 
    bool fFirstXYZfit;
    /// The X coordinate fit is always done in the XZ plane and not in YX.  
    bool fForceXZfit;
    /// Do final simultaneous fit in XYZ
    bool fFinalXYZfit;
    /// Refit with the Sigma parameter fixed when the convergence failed.
    bool fRefitWithFixedSigma;

    /// Boolean to enable/disable fitting the Y projection
    bool fFitYProj;
    /// Boolean to enable/disable fitting the X projection
    bool fFitXProj;

    /// If fReliableFit is false, the likelihood fit is performed
    /// despite some conditions (fixed momentum to converge, ...) 
    /// and the results can be stored for studies but it should not be used for physics
    bool fReliableFit;

    /// Weight factor for the likelihood to normalize the errors in X.
    double fErrorWeightX; 
    /// Weight factor for the likelihood to normalize the errors in Y.
    double fErrorWeightY; 

    /// Boolean to enable/disable multiple scattering error correction during the fit.  
    bool fEnabledMScCorrection;
    /// Weight factor for the likelihood to correct for multiple scattering in X.
    double fErrorWeightXMSc; 
    /// Weight factor for the likelihood to correct for multiple scattering in Y.
    double fErrorWeightYMSc; 
    /// Radiation length in He
    double fRadLenHe;
    /// Muon mass in MeV
    double fMuonMass;
    /// Minimum fractional change of the momentum between the seed and the fit momenta required
    /// to redo the minimization with the adjusted multiple scattering correction.
    double fMomChangeForNewMSc;
    /// Stored initial value of QoP to decide whether to adjust the MSc correction or not.
    double fInitQoP;

    /// Exclude clusters with waveforms containing many peaks from the fit.
    bool fExcludeClusterWithManyPeaks;
    /// Exclude clusters containing saturated waveforms from the fit.
    bool fExcludeSaturatedClusters;
    /// Exclude clusters at the MM edge from the fit.
    bool fExcludeClusterAtEdge;
    /// Exclude clusters with pads that have suspicious timing such as a large
    /// difference between the peak time and the fitted time.
    bool fExcludeSuspiciousPadTiming;
    /// Exclude clusters with more than fMaxPadsPerCluster waveforms.
    unsigned int fMaxPadsPerCluster;
    /// Minimum number of clusters to consider fitting a track.
    unsigned int fMinNumberOfClusters;
    /// Minimum number of clusters to consider calculating the likelihood on these clusters (for matching)
    unsigned int fMinNumberOfClustersForCalc;

    /// Switch to turn on/off the residual dependent penalty
    bool fApplyResidualDependentPenalty;
    /// Minimal value of the residual where the penalty starts
    double fPenaltyThreshold;
    /// Quadratic coefficient for the penalty term
    double fPenaltyCoefficient;


    ///
    std::vector<trex::TTRExHVCluster*> fFitClu;

    /// Value of the total likelihood
  TTPCLogLikelihood fLogLklhd;
    /// Value of the X likelihood
    double fLogLklhdX;
    /// Value of the Horizontal/Vertical cluster likelihood
    double fLogLklhdHV;

    /// Variable filled with Minuit's trial and used in the likelihood calculation.
    double fSigma;

    TTPCPathFitResults fFitResults;

    /// Records for each fit, the stages passed (XYZ fit only ? XZ and YZ fit separated ?)
    unsigned int fFitSteps;
};

// TODO: Find something better
/// Definition of the indexes used in Minuit during the likelihood fit. 
enum fitParams {
  /// X coordinate at the track origin. 
  XPARAM=0,
  /// Y coordinate at the track origin. 
  YPARAM,
  /// Z coordinate at the track origin. 
  ZPARAM,
  /// X coordinate at the track direction. 
  TANXPARAM,
  /// Tangent Y or Z depending if the first cluster is vertical or horizontal.
  TANYORZPARAM,
  /// Track curvature.
  CURVPARAM,
  /// Drift dependent term of the transverse diffusion.
  SGMPARAM };


#endif
