// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SpaceChargeContainer.h
/// \brief Definition of the space-charge density container for the ALICE TPC
/// \author Ernst Hellbär, Goethe-Universität Frankfurt, ernst.hellbar@cern.ch

/*
 * TODO:
 */

#ifndef ALICEO2_TPC_SPACECHARGEINTERFACE_H
#define ALICEO2_TPC_SPACECHARGEINTERFACE_H

#include "AliTPCSpaceCharge3DCalc.h"
#include "DataFormatsTPC/Defs.h"
#include "TPCSimulation/SpaceChargeContainer.h"

class TH3;

namespace o2 {
namespace TPC {

class SpaceChargeInterface
{
  public:

    /// Enumerator for setting the space-charge distortion mode
    enum SCDistortionType {
      SCDistortionsConstant=0,  // space-charge distortions constant over time
      SCDistortionsRealistic=1  // realistic evolution of space-charge distortions over time
    };

    // Constructors
    /// Default constructor using a grid size of (129 z bins, 180 phi bins, 129 r bins)
    SpaceChargeInterface();
    /// Constructor with grid size specified by user
    /// \param nZSlices number of grid points in z, must be (2**N)+1
    /// \param nPhiBins number of grid points in phi
    /// \param nRBins number of grid points in r, must be (2**N)+1
    SpaceChargeInterface(int nZSlices, int nPhiBins, int nRBins);

    // Destructor
    ~SpaceChargeInterface() = default;

    /// Calculate lookup tables if initial space-charge density is provided
    void init();

    /// Set omega*tau and T1, T2 tensor terms in Langevin-equation solution
    /// \param omegaTau omega*tau
    /// \param t1 T1 tensor term
    /// \param t2 T2 tensor term
    void SetOmegaTauT1T2(float omegaTau, float t1, float t2);

    /// Calculate distortion and correction lookup tables using AliTPCSpaceChargeCalc class
    void calculateLookupTables();
    /// Update distortion and correction lookup tables by current space-charge density stored in mSpaceChargeContainer
    void updateLookupTables();

    /// Correct point using correction lookup tables
    /// \param point 3D coordinates to be corrected
    void correctPoint(GlobalPosition3D &point);
    /// Distort point using distortion lookup tables
    /// \param point 3D coordinates to be distorted
    void distortPoint(GlobalPosition3D &point);

    /// Set the space-charge distortions model
    /// \param distortionType distortion type (constant or realistic)
    void setSCDistortionType(SCDistortionType distortionType) {mSCDistortionType = distortionType;}

    /// Set an initial space-charge density
    /// \param hisSCDensity 3D space-charge density histogram, expected format (phi,r,z)
    void setInitialSpaceChargeDensity(TH3 *hisSCDensity);

  private:
    const int mNZSlices;
    const int mNPhiBins;
    const int mNRBins;

    bool mUseInitialSCDensity;  ///< Flag for the use of an initial space-charge density at the beginning of the simulation
    bool mInitLookUpTables; ///< Flag to indicate if lookup tables have been calculated
    SCDistortionType mSCDistortionType;  ///< Type of space-charge distortions

    AliTPCSpaceCharge3DCalc mLookUpTableCalculator;  ///< object to calculate and store correction and distortion lookup tables
    SpaceChargeContainer mSpaceChargeContainer;  ///<  object to handle space-charge density and ion transport
};

} // namespace TPC
} // namespace o2

#endif // ALICEO2_TPC_SPACECHARGEINTERFACE_H
