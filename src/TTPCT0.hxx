#ifndef TTPCT0_hxx_seen
#define TTPCT0_hxx_seen

#include <THit.hxx>
#include <TReconBase.hxx>
#include <recpack/State.h>
// #include <TRecPackManager.hxx>

#define DEFT0 1.e+20

enum TTRExT0Source {kNoT0src = 0, kFGD1T0src, kFGD2T0src, kTopBrECALT0src, kBotBrECALT0src, kNegXBrECALT0src, kPosXBrECALT0src, kDsECALT0src, kP0DT0src, kSMRDT0src, kCathodeT0src};

std::string ConvertT0idxToName(int source);

/// Handy container to store all the information relative to a T0 result
class TTPCT0 {
  public:
    TTPCT0(double Default = DEFT0);
    TTPCT0(TTPCT0 &T0Copy);
    ~TTPCT0() {};

    void Init();
    void LoadMatchedHit(ND::THandle<ND::THit> hitRes, double &chi2, HyperVector &Residuals);
    void LoadClosestHit(ND::THandle<ND::THit> hitRes, double &chi2);
    void LoadCathodeHit(ND::THandle<ND::THit> hitRes, double T0);
    void LoadTimeRange(double MiddleT0, double DeltaT0);
    ND::THandle<ND::THit> GetHit(){return fHit;};
    /// Used only to set a different default T0 in TTPCJunction for example.
    void SetDefaultT0(double t0){fT0 = t0;};
    double GetT0(){return fT0;};
    TTRExT0Source GetSource(){return fSource;};
    void GetT0Range(double &lower, double &upper);
    bool HitTimeWithinT0Range(ND::THandle<ND::THit> candHit);
    double GetChi2(){return fChi2;};
    double GetHitCharge();
    bool ClosestFound(){return fClosestHitFnd;};
    bool MatchFound(){return fMatchedHitFnd;};
    bool OneHitFound(){return (fMatchedHitFnd || fClosestHitFnd);};

    void FillTRealData(ND::THandle<ND::TReconBase> TRecB);

  private:
    bool fClosestHitFnd;
    bool fMatchedHitFnd;
    TTRExT0Source fSource;
    double fChi2;
    double fT0;
    double fT0Range[2];
    ND::THandle<ND::THit> fHit;
    HyperVector fRes;
};


#endif
