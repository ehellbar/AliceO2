/// \file SCContainer.cxx
/// \brief Definition of the ALICE TPC space-charge container
/// \author Ernst Hellbaer, Goethe-Universitaet Frankfurt, hellbaer@ikf.uni-frankfurt.de

#include "TMath.h"

#include "FairLogger.h"

#include "TPCBase/Defs.h"
#include "TPCSimulation/AliTPCPoissonSolver.h"
#include "TPCSimulation/AliTPCSpaceCharge3DDriftLine.h"
#include "TPCSimulation/Constants.h"
#include "TPCSimulation/SCContainer.h"

using namespace o2::TPC;

o2::TPC::SCContainer::SCDistModel o2::TPC::SCContainer::sSCDistModel = o2::TPC::SCContainer::SCDistModel::SCDistOff;

SCContainer::SCContainer() :
  mInitialSCDensityHistogram(nullptr),
  mLookupTables(kFALSE)
{
/// Default constructer
  mPoissonSolver = new AliTPCPoissonSolver();
  mSpaceCharge = new AliTPCSpaceCharge3DDriftLine(65,65,72,5,3,1,3); // TODO: defaults (so far,desired) nRRows=nZColumns=(17,129), nPhiSlices=(18,144), interpolationorder=(1,5), irregulargridsize=(3,?), strategy=1(kUseInterpolator), rbfKernelType(no default)
  mSpaceCharge->SetOmegaTauT1T2(OMEGATAU,1.,1.);
  mSpaceCharge->SetPoissonSolver(mPoissonSolver);
}

SCContainer::SCContainer(Int_t nRRows, Int_t nZColumns, Int_t nPhiSlices, Int_t interpolationorder) :
  mInitialSCDensityHistogram(nullptr),
  mLookupTables(kFALSE)
{
/// Constructer with input options
/// \param nRRows number of rows (radial direction) of the distortion/correction lookup tables  
/// \param nZColumns number of columns (z direction) of the distortion/correction lookup tables
/// \param nPhiSlices number of slices (phi direction) of the distortion/correction lookup tables
/// \param interpolationorder order of cubic spline interpolation functions which are used to interpolate the lookup tables
  mPoissonSolver = new AliTPCPoissonSolver();
  mSpaceCharge = new AliTPCSpaceCharge3DDriftLine(nRRows,nZColumns,nPhiSlices,interpolationorder,3,1,3); // TODO: default values for irregulargridsize, strategy (should be 1=kUseInterpolator), rbfKernelType
  mSpaceCharge->SetOmegaTauT1T2(OMEGATAU,1.,1.);
  mSpaceCharge->SetPoissonSolver(mPoissonSolver);
}

SCContainer::~SCContainer()
{
  delete mInitialSCDensityHistogram;
  delete mPoissonSolver;
  delete mSpaceCharge;
}

void SCContainer::calculateLookupTables()
{
/// Calculate the distortion/correction lookup tables
  if (sSCDistModel==o2::TPC::SCContainer::SCDistModel::ConstSCDist && !mInitialSCDensityHistogram && !mLookupTables) {
    LOG(INFO) << "ATTENTION: Neither an initial space-charge density distribution nor precalculated lookup tables are provided, but required for constant space-charge distortions! Simulation of space-charge distortions is skipped..." << FairLogger::endl;
    return;
  }
  if (!mLookupTables){
    Int_t nRRows = mSpaceCharge->GetNRRows();
    Int_t nZColumns = mSpaceCharge->GetNZColumns();
    Int_t nPhiSlices = mSpaceCharge->GetNPhiSlices();
    mSpaceCharge->InitSpaceCharge3DPoissonIntegralDz(nRRows,nZColumns,nPhiSlices,300,1e-8);
    mLookupTables = kTRUE;
  }
}

void SCContainer::recalculateLookupTables()
{
/// Force a recalculation of the distortion/correction lookup tables
  mLookupTables = kFALSE;
  calculateLookupTables();
}

void SCContainer::distortPoint(GlobalPosition3D &point)
{
/// Distort space-point using precalculated distortion lookup tables
/// \param point space-point coordinates (x,y,z) to be distorted
  const Float_t x[3] = {point.getX(),point.getY(),point.getZ()};
  Float_t dx[3] = {0.,0.,0.};
  if (mLookupTables){
    Float_t phi = TMath::ATan2(x[1],x[0]);
    if (phi<0) phi += TMath::TwoPi();
    Int_t roc = phi/TMath::Pi()*9;
    if (x[2]<0) roc += 18;
    mSpaceCharge->GetDistortion(x,roc,dx);
  }
  point = GlobalPosition3D(x[0]+dx[0],x[1]+dx[1],x[2]+dx[2]);
}

void SCContainer::correctPoint(GlobalPosition3D &point)
{
/// Correct space-point using precalculated correction lookup tables
/// \param point space-point coordinates (x,y,z) to be corrected
  const Float_t x[3] = {point.getX(),point.getY(),point.getZ()};
  Float_t dx[3] = {0.,0.,0.};
  if (mLookupTables){
    Float_t phi = TMath::ATan2(x[1],x[0]);
    if (phi<0) phi += TMath::TwoPi();
    Int_t roc = phi/TMath::Pi()*9;
    if (x[2]<0) roc += 18;
    mSpaceCharge->GetCorrection(x,roc,dx);
  }
  point = GlobalPosition3D(x[0]+dx[0],x[1]+dx[1],x[2]+dx[2]);
}

void SCContainer::setInitialSCDensity(TH3 *scDensity)
{
  /// Set a static space-charge density distribution used for the calculation of the distortion/correction lookup tables
  /// \param scDensity space-charge density distribution in TH3 format (phi,r,z)
  mInitialSCDensityHistogram = scDensity;
  mSpaceCharge->SetInputSpaceCharge(mInitialSCDensityHistogram,1);
  mLookupTables = kFALSE;
}

void SCContainer::setSpaceCharge3D(AliTPCSpaceCharge3DDriftLine *spaceCharge)
{
  mSpaceCharge = spaceCharge;
  mInitialSCDensityHistogram = mSpaceCharge->GetInputSpaceChargeHistogram();
  mLookupTables = kTRUE;
}
