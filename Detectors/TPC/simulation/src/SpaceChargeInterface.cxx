// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SpaceChargeInterface.cxx
/// \brief Implementation of the interface for the ALICE TPC space-charge distortions calculations
/// \author Ernst Hellbär, Goethe-Universität Frankfurt, ernst.hellbar@cern.ch

#include "TH3.h"
#include "TMath.h"
#include "TMatrixD.h"

#include "FairRunAna.h"

#include "Field/MagneticField.h"
#include "TPCBase/ParameterGas.h"
#include "TPCSimulation/SpaceChargeInterface.h"

using namespace o2::TPC;

SpaceChargeInterface::SpaceChargeInterface()
  : mNZSlices(129),
    mNPhiBins(180),
    mNRBins(129),
    mUseInitialSCDensity(false),
    mInitLookUpTables(false),
    mSCDistortionType(SpaceChargeInterface::SCDistortionType::SCDistortionsRealistic),
    mLookUpTableCalculator(129, 129, 180, 5, 3, 0),
    mSpaceChargeContainer(129, 180, 129)
{
}

SpaceChargeInterface::SpaceChargeInterface(int nZSlices, int nPhiBins, int nRBins)
  : mNZSlices(nZSlices),
    mNPhiBins(nPhiBins),
    mNRBins(nRBins),
    mUseInitialSCDensity(false),
    mInitLookUpTables(false),
    mSCDistortionType(SpaceChargeInterface::SCDistortionType::SCDistortionsRealistic),
    mLookUpTableCalculator(nRBins, nZSlices, nPhiBins, 5, 3, 0),
    mSpaceChargeContainer(nZSlices, nPhiBins, nRBins)
{
}

void SpaceChargeInterface::init()
{
  /// TODO how to get the magnetic field when FairRunAna is suspended in O2 simulation?
  o2::field::MagneticField *magField = (o2::field::MagneticField*)FairRunAna::Instance()->GetField();
  float bzField = magField->solenoidField();  // magnetic field in kGauss
  /// TODO is there a faster way to get the drift velocity
  const static ParameterGas &gasParam = ParameterGas::defaultInstance();
  float vDrift = gasParam.getVdrift();  // drift velocity in cm/us
  /// TODO fix hard coded values (ezField, t1, t2): export to Constants.h or get from somewhere?
  float ezField = 400.;  // electric field in z in V/cm
  float t1 = 1.;
  float t2 = 1.;
  /// TODO use this parameterization or fixed value(s) from Magboltz calculations?
  float omegaTau = -10. * bzField * vDrift / ezField;
  SetOmegaTauT1T2(omegaTau, t1, t2);
  if (mUseInitialSCDensity) calculateLookupTables();
}

void SpaceChargeInterface::SetOmegaTauT1T2(float omegaTau, float t1, float t2)
{
  mLookUpTableCalculator.SetOmegaTauT1T2(omegaTau, t1, t2);
}

void SpaceChargeInterface::calculateLookupTables()
{
  mLookUpTableCalculator.ForceInitSpaceCharge3DPoissonIntegralDz(mNRBins, mNZSlices, mNPhiBins, 300, 1e-8);
  mInitLookUpTables = true;
}

void SpaceChargeInterface::updateLookupTables()
{
  TMatrixD *spaceChargeA[mNPhiBins];
  TMatrixD *spaceChargeC[mNPhiBins];
  for (int iphi = 0; iphi < mNPhiBins; ++iphi) {
    spaceChargeA[iphi] = new TMatrixD(mNRBins, mNZSlices);
    spaceChargeC[iphi] = new TMatrixD(mNRBins, mNZSlices);
  }
  mSpaceChargeContainer.getSpaceChargeDensity(spaceChargeA, spaceChargeC, mNZSlices, mNPhiBins, mNRBins);
  mLookUpTableCalculator.SetInputSpaceChargeA(spaceChargeA);
  mLookUpTableCalculator.SetInputSpaceChargeC(spaceChargeC);
  calculateLookupTables();
}

void SpaceChargeInterface::correctPoint(GlobalPosition3D &point)
{
  if (!mInitLookUpTables) return;
  const float x[3] = {point.X(), point.Y(), point.Z()};
  float dx[3] = {0.f, 0.f, 0.f};
  float phi = TMath::ATan2(x[1], x[0]);
  if (phi < 0)
    phi += TMath::TwoPi();
  int roc = phi / TMath::Pi() * 9;
  if (x[2] < 0)
    roc += 18;
  mLookUpTableCalculator.GetCorrection(x, roc, dx);
  point = GlobalPosition3D(x[0] + dx[0], x[1] + dx[1], x[2] + dx[2]);
}

void SpaceChargeInterface::distortPoint(GlobalPosition3D &point)
{
  if (!mInitLookUpTables) return;
  const float x[3] = {point.X(), point.Y(), point.Z()};
  float dx[3] = {0.f, 0.f, 0.f};
  float phi = TMath::ATan2(x[1], x[0]);
  if (phi < 0)
    phi += TMath::TwoPi();
  int roc = phi / TMath::Pi() * 9;
  if (x[2] < 0)
    roc += 18;
  mLookUpTableCalculator.GetDistortion(x, roc, dx);
  point = GlobalPosition3D(x[0] + dx[0], x[1] + dx[1], x[2] + dx[2]);
}

void SpaceChargeInterface::setInitialSpaceChargeDensity(TH3 *hisSCDensity)
{
  TMatrixD *spaceChargeA[mNPhiBins];
  TMatrixD *spaceChargeC[mNPhiBins];
  for (int iphi = 0; iphi < mNPhiBins; ++iphi) {
    spaceChargeA[iphi] = new TMatrixD(mNRBins, mNZSlices);
    spaceChargeC[iphi] = new TMatrixD(mNRBins, mNZSlices);
  }
  mLookUpTableCalculator.GetChargeDensity(spaceChargeA, spaceChargeC, hisSCDensity, mNRBins, mNZSlices, mNPhiBins);
  mLookUpTableCalculator.SetInputSpaceChargeA(spaceChargeA);
  mLookUpTableCalculator.SetInputSpaceChargeC(spaceChargeC);
  mSpaceChargeContainer.setInitialSpaceChargeDensity(spaceChargeA, spaceChargeC, mNZSlices, mNPhiBins, mNRBins);
  mUseInitialSCDensity = true;
}
