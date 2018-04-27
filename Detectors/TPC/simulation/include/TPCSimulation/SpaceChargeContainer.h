// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SpaceChargeInterface.h
/// \brief Definition of the interface for the ALICE TPC space-charge distortions calculations
/// \author Ernst Hellbär, Goethe-Universität Frankfurt, ernst.hellbar@cern.ch

/*
 * TODO:
 *   - fix constants (more precise values, export into TPCBase/Constants)
 *   - pad granularity in r, rphi?
 *   - accumulate and add next slice
 *     - event based: propagate charge(ievent-1), add charge(ievent)
 *     - time based: mTime0, mEffectiveTime
 *                   addIon(eventtime, drifttime, r, phi)
 *                   time in us, 50 kHz = <one event / 20 us>
 *                   if (ev.time+dr.time-mTime0 < mLengthTimebin) => add2NextSlice
 *                   if  (mLengthTimebin < ev.time+dr.time-mTime0 < mLengthTimebin+100us) add2NextToNextSlice
 *                   apply updated distortions to ions in NextToNextSlice when space charge is propagated; need to store exact positions (e.g. std::vector<std::vector<float>>)!
 *   - include next slice and transport ions along straight lines
 *   - ion transport along the E field -> Jacobi matrices?
 *   - irregular bin sizes in r and rphi
 *   - timebins or z bins?
 */

#ifndef ALICEO2_TPC_SPACECHARGECONTAINER_H
#define ALICEO2_TPC_SPACECHARGECONTAINER_H

namespace o2 {
namespace TPC {

class SpaceChargeContainer
{
  public:
    /// Default constructor with max. container size
    SpaceChargeContainer();

    /// Constructor to specify container size manually
    /// \param nZBins number of container z bins
    /// \param nR number of container r bins
    /// \param nPhi number of container phi bins
    SpaceChargeContainer(int nZSlices, int nPhiBins, int nRBins);

    /// Store space-charge density in a TMatrixD for A and C side, respectively
    /// \param spaceChargeA output space-charge density A side
    /// \param spaceChargeC output space-charge density C side
    /// \param nZSlices number of z slices
    /// \param nPhiBins number of phi bins
    /// \param nRBins number of r bins
    void getSpaceChargeDensity(TMatrixD **spaceChargeA, TMatrixD **spaceChargeC, int nZSlices, int nPhiBins, int nRBins);

    /// Set an initial space-charge density
    /// \param matricesChargeA space-charge density (C/m^3) on the A side as TMatrixD**[phi](r,z)
    /// \param matricesChargeC space-charge density (C/m^3) on the C side as TMatrixD**[phi](r,z)
    /// \param nZSlices number of z slices
    /// \param nPhiBins number of phi bins
    /// \param nRBins number of r bins
    void setInitialSpaceChargeDensity(TMatrixD **matricesChargeA, TMatrixD **matricesChargeC, int nZSlices, int nPhiBins, int nRBins);

  private:
    /// Convert amount of ions into charge density C/m^3
    /// \param nIons number of ions
    /// \return space-charge density (C/m^3)
    float ions2Charge(int nIons);

    /// # of TPC timebins, 1 TPC timebin = 0.2 us, default number of z bins
    static constexpr int mTPCTIMEBINS = 500;  //!
    /// default number of phi bins
    static constexpr int mMAXPHIBINS = 360;  //!
    /// ion dirft time (TPC timebins)
    static constexpr int mDRIFTTIMEIONS = 8e5;  //!
    ///  drift length of the TPC in (cm)
    static constexpr float mDRIFTLENGTH = 250.;  //!
    /// inner radius of the TPC active area
    static constexpr float mRADIUSINNER = 85.;  //!
    /// outer radius of the TPC active area
    static constexpr float mRADIUSOUTER = 245.;  //!

    /// length of one z bin (cm)
    const float mLengthZSlice; //!
    /// width of one phi bin (radians)
    const float mWidthPhiBin;  //!
    /// length of one r bin (cm)
    const float mLengthRBin;  //!

    std::vector <std::vector<float>> mSpaceChargeDensityA;  ///< current space-charge density on the A side, stored in C/m^3 (z)(phi*r), ordering: z=[0,250], ir+iphi*nRBins
    std::vector <std::vector<float>> mSpaceChargeDensityC;  ///< current space-charge density on the C side, stored in C/m^3 (z)(phi*r), ordering: z=[0,-250], ir+iphi*nRBins
    std::vector<int> mNextSpaceChargeSliceA;  ///< next space-charge density slice on the A side, stored in # ions and accumulated for duration of one z / time bin
    std::vector<int> mNextSpaceChargeSliceC;  ///< next space-charge density slice on the C side, stored in # ions and accumulated for duration of one z / time bin
    std::vector<std::vector<float>> mNextToNextIon;  ///< next to next space-charge density slice on the A side, stored in # ions, containing ions between the end of mNextSpaceChargeSliceA and mNextSpaceChargeSliceA+100us
};

} // namespace TPC
} // namespace o2
#endif // ALICEO2_TPC_SPACECHARGECONTAINER_H
