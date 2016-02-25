#ifndef TTPCQLikelihood_hxx_seen
#define TTPCQLikelihood_hxx_seen

#include "TTPCDebug.hxx"

class TTPCQLikelihood {
  public:

    TTPCQLikelihood();
    ~TTPCQLikelihood();

    double q_exp( double b, double phi, double sigma,double longi,double trans);

  private: 

    /// Minimum value of phi to be considerered during eta calculation with non zero angle.
    double fMinimumPhi_Eta; 

    double eta( double b, double phi, double sigma,double longi,double trans);

};

#endif 
