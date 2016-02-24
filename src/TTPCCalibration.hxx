#ifndef TTPCCalibration_hxx_seen
#define TTPCCalibration_hxx_seen

#define TREXSTDOUTIDOFFSET 1000

#include <TND280Event.hxx>

namespace ND {
  class TTPCCalibration;
}

class ND::TTPCCalibration {
public:

  virtual ~TTPCCalibration() {};

  /// Get a pointer to the singleton instance of the calibration information.
  static TTPCCalibration& Get(void);

  void ReadCalibration(const ND::TND280Event& Event);

  bool IsMC();

  double GetDriftVelocity();

  double GetTimeOffset();

  /// Get the next Path Id available and automatically increment.
  unsigned int GetPathId();
  /// Get the next Junction Id available and automatically increment.
  unsigned int GetJunctionId();
  /// Get the next Pattern Id available and automatically increment.
  unsigned int GetPatternId();

  /// Reset all the Ids to a large offset to mark the difference
  /// between gas output patterns/paths/junctions and those modified
  /// for the standard output
  void OffsetIds();

  /// Set the default T0.
  void SetDefaultT0(double DftT0);
  /// Get the default T0.
  double GetDefaultT0();
  

private:
  TTPCCalibration();

  static TTPCCalibration* _tpcCalibration;

  // Save last event calibrated to avoid calling the calibration table
  // multiple times for the same event.
  int fLastEvent;
  int fLastTime;

  /// Store if the event is MC or not
  bool fEvtIsMC;

  /// Time offset applied to the T0 of each track.
  double fTimeOffset; 
  double fTimeOffset_Data; 
  double fTimeOffset_MC; 
  /// Value of the DriftVelocity from the data base.
  double fDriftVelocity; 
  // Length of the TPC drift volume
  double fMaxDrift;
  /// Indicates if the measured value of the drift velocity is to be used, or the default tpcRecon.TPCgas.DriftSpeed
  bool fUseMeasDriftSpeed;

  /// Default T0 for spill events.
  double fDefaultSpillT0;
  /// Default T0 for cosmic events.
  double fDefaultCosmT0;
  /// Default T0
  double fDefaultT0;

  /// Use this variable to define a unique ID for the path
  /// with a global index.
  unsigned int fPathIdOffset;
  /// Use this variable to define a unique ID for the junction
  /// with a global index.
  unsigned int fJunctionIdOffset;
  /// Use this variable to define a unique ID for the pattern
  /// with a global index.
  unsigned int fPatternIdOffset;
};


namespace ND{
  TTPCCalibration& tpcCalibration();
};

#endif
