/// \file SCContainer.h
/// \brief Definition of the ALICE TPC space-charge container
/// \author Ernst Hellbaer, Goethe-Universitaet Frankfurt, hellbaer@ikf.uni-frankfurt.de

#ifndef O2_TPC_SCContainer_H_
#define O2_TPC_SCContainer_H_

#include "TPCSimulation/AliTPCPoissonSolver.h"
#include "TPCSimulation/AliTPCSpaceCharge3DDriftLine.h"

class TH3;

namespace o2 {
namespace TPC {

class Defs;
  
class SCContainer
{
public:

  /// Enumerator for different types of space-charge distortion models
  enum SCDistModel {
    SCDistOff=0,		// no space-charge distortions
    ConstSCDist=1,	// space-charge distortions constant over time
    RealSCDist=2	// realistic space-charge distortions
  };

  // Constructors
  SCContainer();
  SCContainer(Int_t nRRows, Int_t nZColumns, Int_t nPhiSlices, Int_t interpolationorder);
  
  // Destructor 
  ~SCContainer() = default;

  //
  void setSpaceCharge3D(AliTPCSpaceCharge3DDriftLine &spaceCharge);
  const AliTPCSpaceCharge3DDriftLine &getSpaceCharge3D() const {return mSpaceCharge;}

  // Calculation of distortion/correction lookup tables
  void calculateLookupTables();
  void recalculateLookupTables();

  // Distortion/correction operations
  void distortPoint(GlobalPosition3D &point);
  void correctPoint(GlobalPosition3D &point);

  /// Define the space-charge distotion model
  /// \param scDistType space-charge distortion model applied in simulation 
  static void setSCDistortionsModel(o2::TPC::SCContainer::SCDistModel scModel) {sSCDistModel = scModel;}
  /// Get the space-charge distotion model
  /// \return space-charge distortion model applied in simulation 
  static SCDistModel getSCDistortionsModel() {return sSCDistModel;}

  //
  void setInitialSCDensity(TH3 *scDensity);
 
private:
  
  static SCDistModel sSCDistModel; ///< variable to define type of space-charge distortions applied
  Bool_t mLookupTables;		   ///< flag to keep track if distortion/correction lookup tables have been calculated

  TH3 *mInitialSCDensityHistogram;

  AliTPCPoissonSolver mPoissonSolver;
  AliTPCSpaceCharge3DDriftLine mSpaceCharge;

};

}
}

#endif	// O2_TPC_SCContainer_H_

