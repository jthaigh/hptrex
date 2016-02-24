#ifndef TTPCPattern_hxx_seen
#define TTPCPattern_hxx_seen

#include <THandle.hxx>
#include <TReconVertex.hxx>

#include <TTPCPath.hxx>
#include <TTPCJunction.hxx>

namespace ND {
  class TTPCPattern;
}

/// The constituents of the pattern are dynamic, they change at each step.
/// For example after the matching across the vertical MM gap of two patterns,
/// we remove some TReconPIDs, and add the junctions of the second pattern
/// to the first.
class ND::TTPCPattern : public TReconVertex {
public:
  TTPCPattern();
  TTPCPattern(ND::THandle<TTPCPath> path);
  virtual ~TTPCPattern();

  void SetId(unsigned int theId);
  unsigned int GetId();

  /// The pattern is constructed by adding junctions.
  /// Junctions are the base object, not paths.
  /// The paths connected to this junction are added
  /// to the list of constituents of the pattern automatically.
  void AddJunction(ND::THandle<TTPCJunction> junction);

  /// Setup things like the detector bit for the path.
  /// This should be called whenever the pattern has all its Junctions and Paths.
  void InitialSetup();

  /// Count the number of TTPCPaths in the list of constituents
  unsigned int GetNPaths();
  /// Count the number of TTPCJunctions in the list of constituents
  unsigned int GetNJunctions();

  void SetUsable(bool usable);
  bool IsUsable();

  /// Return the TTPCPaths in the list of constituents
  std::vector< ND::THandle<ND::TTPCPath> > GetPaths();
  /// Return the TTPCJunctions in the list of constituents
  std::vector< ND::THandle<ND::TTPCJunction> > GetJunctions();

  /// Convenient method to get all the hits of this pattern.
  ND::THandle<ND::THitSelection> GetConstituentHits();

  void SetT0(TTPCT0 &T0);
  bool HasT0();
  TTRExT0Source GetT0Source();
  double GetT0();

  ND::THandle<ND::TReconBase> ConvertToOAEvent();

private:
  bool fUsable;

  unsigned int fId;
  TTPCT0 fT0;

  ClassDef(TTPCPattern,1);
};


#endif
