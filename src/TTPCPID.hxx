#ifndef TTPCPID_hxx_seen
#define TTPCPID_hxx_seen

#include <string.h>

#include <THandle.hxx>
#include <TOADatabase.hxx>
#include <TComboHit.hxx>

#include <TGeomInfo.hxx>

#include "TTPCPattern.hxx"
#include "TTPCPath.hxx"

namespace ND {
  class TTPCPID;
}

// TODO: A bit more comments would be nice.
/// Cluster hits in a single plane by charge.
class ND::TTPCPID {
  public:
    TTPCPID();

    virtual ~TTPCPID();

    /// Do the PID.

    void Process(ND::THandle<ND::TTPCPattern> Pattern);

    void TrackdEdx(ND::THandle<ND::TTPCPath> Path);

    Double_t ExpecteddEdx(double bg);
    Double_t TrackMomError(double momentum, double momerr, double Mparticle);
    void WritePiDInfo(ND::THandle<ND::TTPCPath> path);

    void ClusterSelection(ND::THandle<ND::TComboHit>& m_combo, ND::TTPCGeom& m_tpcGeom, std::string *m_cltype, int &m_nSample, bool *isGoodCluster);

    // TODO: How about using an array for those ?
    // Then use SetfExpecteddE(int index, double val);
    void SetfExpecteddE0(double val) { fExpecteddE0 = val; }
    void SetfExpecteddE1(double val) { fExpecteddE1 = val; }
    void SetfExpecteddE2(double val) { fExpecteddE2 = val; }
    void SetfExpecteddE3(double val) { fExpecteddE3 = val; }
    void SetfExpecteddE4(double val) { fExpecteddE4 = val; }

    double GetLikelihood(double x);

  private:

    void Reset(void);
    void Init(void);

    // TODO: Reduce the number of these variables. Surely some of them don't need to be
    // class variables but can be local to methods.
    // TODO: Use arrays if possible.
    // TODO: Make it clear if a variable should be applied to vertical clusters only
    // or horizontal clusters only in the name of the variable.
    bool fVerbose;

    bool fUseHeightBasedCharge;
    double fHeightToAreaFactor;

    double fMIPmuon;
    double fSigmaMIP;
    double fSLhor;
    double fSLver;
    double fPadFraction;

    double fElectronCorr;

    double fChargeVsDriftParameter;
    double fDriftCenter;
    int fDriftSense;

    double fExpecteddE0;
    double fExpecteddE1;
    double fExpecteddE2;
    double fExpecteddE3;
    double fExpecteddE4;

    double fScaleFactor;

    double fSLCorrH0;
    double fSLCorrH1;
    double fSLCorrH2; //for HorClu

    double fRowCorrF0;
    double fRowCorrF1;
    double fRowCorrF2; //for HorClu
    double fRowCorrF3; //for HorClu
    double fRowCorrF4; //for HorClu
    double fRowCorrF5; //for HorClu

    double fSigmaSLCorrA0;
    double fSigmaSLCorrA1;
    double fSigmaSLCorrA2; //for HorClu
    double fSigmaSLCorrA3; //for HorClu

    double fSigmaRowCorrP0; //for HorClu
    double fSigmaRowCorrP1;
    double fSigmaRowCorrP2; //for HorClu
    double fSigmaRowCorrP3; //for HorClu

    //double fSigmaRowCorrP2;

    double fMomentum;
    double fMomentumErr;
    double fTotCharge;
    double fDedxmeas;
    double fDedxcorr;
    double SLmean;
    int fNSample;
    int fNTrunMM;

    double fDedxexp[5];
    double fSigmaexp[5];
    double fPull[5];
};

#endif
