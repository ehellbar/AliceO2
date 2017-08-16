/// \file Digitizer.h
/// \brief Definition of the ALICE TPC digitizer task
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef ALICEO2_TPC_DigitizerTask_H_
#define ALICEO2_TPC_DigitizerTask_H_

#include <cstdio>
#include <string>
#include "FairTask.h"
#include "FairLogger.h"
#include "TPCSimulation/Digitizer.h"
#include "TPCSimulation/SCContainer.h"
#include <TClonesArray.h>

class TH3;

namespace o2 {
namespace TPC { 

class AliTPCSpaceCharge3DDriftLine;
class Digitizer;

/// \class DigitizerTask
/// This task steers the digitization process and takes care of the input and output
/// Furthermore, it allows for switching debug output on/off
    
class DigitizerTask : public FairTask{
public:
      
  /// Default constructor
  DigitizerTask();
      
  /// Destructor
  ~DigitizerTask() override;
      
  /// Inititializes the digitizer and connects input and output container
  InitStatus Init() override;

  void setHitFileName(std::string name) { mHitFileName = name; }

  /// Sets the debug flags for the sub-tasks
  /// \param debugsString String containing the debug flags
  ///        o PRFdebug - Debug output after application of the PRF
  void setDebugOutput(TString debugString);

  /// Switch for triggered / continuous readout
  /// \param isContinuous - false for triggered readout, true for continuous readout
  void setContinuousReadout(bool isContinuous);

  /// Set the maximal number of written out time bins
  /// \param nTimeBinsMax Maximal number of time bins to be written out
  void setMaximalTimeBinWriteOut(int i) { mTimeBinMax = i; }

  /// Set space-charge model used in simulation
  /// \param scModel type o2::TPC::SCContainer::SCDistModel (SCDistOff,ConstSCDist,RealSCDist)
  void setSCDistortionsModel(o2::TPC::SCContainer::SCDistModel scModel);

  // Setters for initial space charge distribution
  void setInitialSpaceCharge(TH3 *scDensity);
  void setInitialSpaceCharge(AliTPCSpaceCharge3DDriftLine *spaceCharge3D);
      
  /// Digitization
  /// \param option Option
  void Exec(Option_t *option) override;
      
  void FinishTask() override;

private:
  void fillHitArrayFromFile();

  Digitizer           *mDigitizer;    ///< Digitization process
  DigitContainer      *mDigitContainer;
      
  TClonesArray        *mPointsArray;  ///< Array of detector hits, passed to the digitization
  TClonesArray        *mDigitsArray;  ///< Array of the Digits, passed from the digitization
    
  std::string         mHitFileName;  ///< External hit file exported from AliRoot

  int                 mTimeBinMax;   ///< Maximum time bin to be written out
  bool                mIsContinuousReadout; ///< Switch for continuous readout

  ClassDefOverride(DigitizerTask, 1);
};

inline
void DigitizerTask::setDebugOutput(TString debugString)
{
  LOG(INFO) << "TPC - Debug output enabled for: ";
  if (debugString.Contains("PRFdebug")) {
    LOG(INFO) << "Pad response function, ";
    o2::TPC::Digitizer::setPRFDebug();
  }
  LOG(INFO) << "\n";
}
  
inline
void DigitizerTask::setContinuousReadout(bool isContinuous)
{
  mIsContinuousReadout = isContinuous;
  o2::TPC::Digitizer::setContinuousReadout(isContinuous);
}

inline
void DigitizerTask::setSCDistortionsModel(o2::TPC::SCContainer::SCDistModel scModel)
{
  o2::TPC::SCContainer::setSCDistortionsModel(scModel);
}

inline
void DigitizerTask::setInitialSpaceCharge(TH3 *scDensity)
{
  /// Set initial space-charge density
  /// \param scDensity TH3, format (phi,r,z)
  mDigitizer->setInitialSpaceCharge(scDensity);
}

inline
void DigitizerTask::setInitialSpaceCharge(AliTPCSpaceCharge3DDriftLine *spaceCharge3D)
{
  /// Set precalculated lookup tables as initial space charge
  /// \param spaceCharge3D AliTPCSpaceCharge3DDriftLine object with precalculated lookup tables
  mDigitizer->setInitialSpaceCharge(spaceCharge3D);
}
  
}
}

#endif // ALICEO2_TPC_DigitizerTask_H_
