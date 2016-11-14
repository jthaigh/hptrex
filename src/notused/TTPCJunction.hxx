#ifndef TTPCJunction_hxx_seen
#define TTPCJunction_hxx_seen

#include <THandle.hxx>
#include <TReconVertex.hxx>
#include <TVertexState.hxx>

#include "TTPCT0.hxx"

namespace ND {
  class TTPCJunction;
}


class ND::TTPCJunction : public TReconVertex {
public:
  TTPCJunction();
  TTPCJunction(const TLorentzVector &Position);
  virtual ~TTPCJunction();

  void SetId(unsigned int theId);
  unsigned int GetId();

  /// Setup things like the detector bit for the path.
  /// This should be called right after AddHits.
  void InitialSetup();

  void SetT0(TTPCT0 &T0);
  bool HasT0();
  TTRExT0Source GetT0Source();
  double GetT0();

  /// Calculate the X position of the given hit using the junction's T0
  double GetCalibX(ND::THandle<ND::THit> Hit);

  /// Get number of paths associated with this junction
  unsigned int GetNPaths(){ return (unsigned int)(GetConstituents()->size()); }

  /// Check if a path with this Id is connected to this junction
  bool IsPathConnected(unsigned int WantedPathId);

  ND::THandle<ND::TReconVertex> ConvertToOAEvent();

private:
  unsigned int fId;

  /// T0 information that can be applied to all the hits.
  TTPCT0 fT0;

  ClassDef(TTPCJunction,1);
};


#endif
