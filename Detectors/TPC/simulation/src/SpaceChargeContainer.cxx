// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SpaceChargeContainer.cxx
/// \brief Implementation of the space-charge density container for the ALICE TPC
/// \author Ernst Hellbär, Goethe-Universität Frankfurt, ernst.hellbar@cern.ch

#include "TH3.h"
#include "TMath.h"
#include "TMatrixD.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsTPC/Constants.h"
#include "DataFormatsTPC/Defs.h"

#include "TPCSimulation/SpaceChargeContainer.h"

using namespace o2::TPC;

SpaceChargeContainer::SpaceChargeContainer()
  : mLengthZSlice(mDRIFTLENGTH / mTPCTIMEBINS),
    mWidthPhiBin(o2::constants::math::TwoPI / mMAXPHIBINS),
    mLengthRBin((mRADIUSOUTER - mRADIUSINNER) / Constants::MAXGLOBALPADROW),
    mSpaceChargeDensityA(mTPCTIMEBINS),
    mSpaceChargeDensityC(mTPCTIMEBINS),
    mNextSpaceChargeSliceA(mMAXPHIBINS * Constants::MAXGLOBALPADROW),
    mNextSpaceChargeSliceC(mMAXPHIBINS * Constants::MAXGLOBALPADROW)
{
  for (auto i = 0; i < SpaceChargeContainer::mTPCTIMEBINS; ++i) {
    mSpaceChargeDensityA[i].resize(mMAXPHIBINS * Constants::MAXGLOBALPADROW);
    mSpaceChargeDensityC[i].resize(mMAXPHIBINS * Constants::MAXGLOBALPADROW);
  }
}

SpaceChargeContainer::SpaceChargeContainer(int nZSlices, int nPhiBins, int nRBins)
  : mLengthZSlice(mDRIFTLENGTH / nZSlices),
    mWidthPhiBin(TWOPI / nPhiBins),
    mLengthRBin((mRADIUSOUTER - mRADIUSINNER) / nRBins),
    mSpaceChargeDensityA(nZSlices),
    mSpaceChargeDensityC(nZSlices),
    mNextSpaceChargeSliceA(nPhiBins * nRBins),
    mNextSpaceChargeSliceC(nPhiBins * nRBins)
{
  for (auto i = 0; i < nZSlices; ++i) {
    mSpaceChargeDensityA[i].resize(nPhiBins * nRBins);
    mSpaceChargeDensityC[i].resize(nPhiBins * nRBins);
  }
}

void SpaceChargeContainer::getSpaceChargeDensity(TMatrixD **spaceChargeA, TMatrixD **spaceChargeC, int nZSlices, int nPhiBins, int nRBins)
{
  for (int iside=0;iside<2;++iside){
    for (int iphi=0;iphi<nPhiBins;++iphi){
      TMatrixD &chargeDensity = iside==0 ? *spaceChargeA[iphi] : *spaceChargeC[iphi];
      for (int ir=0;ir<nRBins;++ir){
        for (int iz=0;iz<nZSlices;++iz){
          if (iside==0){
            chargeDensity(ir,iz) = mSpaceChargeDensityA[iz][ir+iphi*nRBins];
          } else {
            chargeDensity(ir,iz) = mSpaceChargeDensityC[iz][ir+iphi*nRBins];
          }
        }
      }
    }
  }
}

void SpaceChargeContainer::setInitialSpaceChargeDensity(TMatrixD **matricesChargeA, TMatrixD **matricesChargeC, int nZSlices, int nPhiBins, int nRBins)
{
  for (int iside=0;iside<2;++iside){
    for (int iphi=0;iphi<nPhiBins;++iphi){
      TMatrixD &chargeDensity = iside==0 ? *matricesChargeA[iphi] : *matricesChargeC[iphi];
      for (int ir=0;ir<nRBins;++ir){
        for (int iz=0;iz<nZSlices;++iz){
          if (iside==0){
            mSpaceChargeDensityA[iz][ir+iphi*nRBins] = chargeDensity(ir,iz);
          } else {
            mSpaceChargeDensityC[iz][ir+iphi*nRBins] = chargeDensity(ir,iz);
          }
        }
      }
    }
  }
}

float SpaceChargeContainer::ions2Charge(int nIons)
{
  return 1.e6f * nIons * TMath::Qe() / (mLengthZSlice * (mWidthPhiBin * mLengthRBin) * mLengthRBin);
}
