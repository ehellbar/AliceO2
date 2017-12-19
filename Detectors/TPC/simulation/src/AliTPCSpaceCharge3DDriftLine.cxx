/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TGeoGlobalMagField.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector.h"
#include "TVector3.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TPCSimulation/AliTPCPoissonSolver.h"
#include "TPCSimulation/AliTPCSpaceCharge3DDriftLine.h"
// added
#include "FairLogger.h"
#include "TPCBase/Mapper.h"

// removed
/*
  #include "AliMagF.h"
  #include "AliTPCcalibDB.h" 
  #include "AliTPCParam.h" 
  #include "AliLog.h"
  #include "AliTPCParam.h"
  #include "AliTPCParamSR.h"
  #include "AliSysInfo.h"
  #include "AliTPCROC.h"
  // #include "AliTPCSpaceCharge3D.h"
  */

/// \cond CLASSIMP
ClassImp(o2::TPC::AliTPCSpaceCharge3DDriftLine)
/// \endcond

using namespace o2::TPC;


/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Default values
/// ~~~
/// fInterpolationOrder = 1; //triliniear 3d //2 there is an error
/// fStrategy = kUseInterpolator;	
/// fNRRows = 72;
/// fNPhiSlices = 180; // the maximum of phi-slices so far = (8 per sector)
/// fNZColumns = 166; // the maximum on column-slices so  ~ 2cm slicing
/// ~~~
///  
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine()
  // : AliTPCCorrection(),fC0(0.),fC1(0.),fCorrectionFactor(1.),  fInitLookUp(kFALSE) // removed
  : TNamed("SpaceCharge3D","SpaceCharge3D"),fT1(1.),fT2(1.),fC0(0.),fC1(0.),fCorrectionFactor(1.),  fInitLookUp(kFALSE) // added
{
	
  fInterpolationOrder = 1; //triliniear 3d //2 there is an error
  fIrregularGridSize = 3;
  fStrategy = kUseInterpolator;	
  fNRRows = 17;
  fNPhiSlices = 18; // the maximum of phi-slices so far = (8 per sector)
  fNZColumns = 17; // the maximum on column-slices so  ~ 2cm slicing
  
  fLookupRList = new Double_t[fNRRows];
  fLookupPhiList = new Double_t[fNPhiSlices];
  fLookupZList = new Double_t[fNZColumns];
  
  fLookupZListA = new Double_t[fNZColumns/2];
  fLookupZListC = new Double_t[fNZColumns/2];
  
  
  
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (fNRRows - 1) ;
  const Float_t gridSizeZ   =  (o2::TPC::AliTPCPoissonSolver::fgkTPCZ0) / (fNZColumns/2  - 1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / fNPhiSlices;
	
  
	
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) fLookupPhiList[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < fNRRows ; i++ )  fLookupRList[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < fNZColumns; j++ ) fLookupZList[j]  = (j * gridSizeZ)  - o2::TPC::AliTPCPoissonSolver::fgkTPCZ0;
	
  for ( Int_t j = 0 ; j < fNZColumns/2; j++ ) 
    {
      fLookupZListA[j]  = (j * gridSizeZ);
		
    }
	
  // negative 
  for ( Int_t j = 0 ; j < fNZColumns/2; j++ ) fLookupZListC[j]  = -(j * gridSizeZ);
	
	
	
	
  // Array which will contain the solution according to the setted charge density distribution
  // see InitSpaceCharge3DDistortion() function
  printf("allocating look up %d slices, (%d,%d), interp order : %d\n",fNPhiSlices,fNRRows,fNZColumns,fInterpolationOrder);
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    fLookUpIntDistDrEz[k]   =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntDistDphiREz[k] =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntDistDz[k]    =  new TMatrixD(fNRRows,fNZColumns);
    
    fLookUpIntDistDrEzA[k]   =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDphiREzA[k] =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDzA[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
    
    fLookUpIntDistDrEzC[k]   =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDphiREzC[k] =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDzC[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
    
    
    fLookUpIntCorrDrEz[k]   =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntCorrDphiREz[k] =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntCorrDz[k]    =  new TMatrixD(fNRRows,fNZColumns);
    
    
    /// no drift
		
    fLookUpErOverEz[k]    =  new TMatrixD(fNRRows,fNZColumns);;
    fLookUpEphiOverEz[k]  =  new TMatrixD(fNRRows,fNZColumns);;
    fLookUpDeltaEz[k]     =  new TMatrixD(fNRRows,fNZColumns);;

    fLookUpErOverEzA[k]    =  new TMatrixD(fNRRows,fNZColumns/2);;
    fLookUpEphiOverEzA[k]  =  new TMatrixD(fNRRows,fNZColumns/2);;
    fLookUpDeltaEzA[k]     =  new TMatrixD(fNRRows,fNZColumns/2);;


    /// no drift
    fLookUpErOverEzC[k]    =  new TMatrixD(fNRRows,fNZColumns/2);;
    fLookUpEphiOverEzC[k]  =  new TMatrixD(fNRRows,fNZColumns/2);;
    fLookUpDeltaEzC[k]     =  new TMatrixD(fNRRows,fNZColumns/2);;

    
    
    fSCdensityDistribution[k] = new TMatrixF(fNRRows,fNZColumns);

    // irregular grid for correction 
	
    // Correction with irregular interpolation
    fLookUpIntCorrDrEzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);   		//[kNPhi]  
    fLookUpIntCorrDphiREzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);   //[kNPhi]  
    fLookUpIntCorrDzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fRListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);;
    fPhiListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);; 
    fZListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);;   

    fLookUpIntCorrDrEzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);   		//[kNPhi]  
    fLookUpIntCorrDphiREzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);   //[kNPhi]  
    fLookUpIntCorrDzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fRListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);;
    fPhiListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);; 
    fZListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);;   

    // lookup for charge
    fLookUpChargeA[k]  = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeC[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeInverseA[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeInverseC[k] = new TMatrixD(fNRRows,fNZColumns/2);


  }

  fLookupIntDist = 
    new AliTPCLookUpTable3DInterpolatorD(	
					 fNRRows,fLookUpIntDistDrEz,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREz,fLookupPhiList,fNZColumns,fLookUpIntDistDz,fLookupZList,fInterpolationOrder);
			
  fLookupIntDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEzA,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREzA,fLookupPhiList,fNZColumns/2,fLookUpIntDistDzA,fLookupZListA,fInterpolationOrder);		
	
  fLookupIntDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEzC,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREzC,fLookupPhiList,fNZColumns/2,fLookUpIntDistDzC,fLookupZListC,fInterpolationOrder);		
			
  fLookupIntCorr = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntCorrDrEz,fLookupRList,fNPhiSlices,fLookUpIntCorrDphiREz,fLookupPhiList,fNZColumns,fLookUpIntCorrDz,fLookupZList,fInterpolationOrder);


  fLookupIntENoDrift =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEz,fLookupRList,fNPhiSlices,fLookUpEphiOverEz,fLookupPhiList,fNZColumns,fLookUpDeltaEz,fLookupZList,fInterpolationOrder);		
	
  fLookupIntENoDriftA =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEzA,fLookupRList,fNPhiSlices,fLookUpEphiOverEzA,fLookupPhiList,fNZColumns/2,fLookUpDeltaEzA,fLookupZListA,fInterpolationOrder);		
	
	
  fLookupIntENoDriftC =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEzC,fLookupRList,fNPhiSlices,fLookUpEphiOverEzC,fLookupPhiList,fNZColumns/2,fLookUpDeltaEzC,fLookupZListC,fInterpolationOrder);		
	
	
  // look up table for correction error based on irregular grid
  fLookupIntCorrIrregularA =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     fNRRows,
     fLookUpIntCorrDrEzIrregularA,
     fRListIrregularA,
     fLookupRList,
     fNPhiSlices,
     fLookUpIntCorrDphiREzIrregularA,
     fPhiListIrregularA,
     fLookupPhiList,
     fNZColumns/2,
     fLookUpIntCorrDzIrregularA,
     fZListIrregularA,
     fLookupZListA,
     2,
     3,
     3,
     3,
     1
     );
		
  fLookupIntCorrIrregularC =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     fNRRows,
     fLookUpIntCorrDrEzIrregularC,
     fRListIrregularC,
     fLookupRList,
     fNPhiSlices,
     fLookUpIntCorrDphiREzIrregularC,
     fPhiListIrregularC,
     fLookupPhiList,
     fNZColumns/2,
     fLookUpIntCorrDzIrregularC,
     fZListIrregularC,
     fLookupZListC,
     2,
     3,
     3,
     3,
     1
     );



  // bind lookup matrices to lookup table

  fInterpolatorChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorChargeC = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeC = new AliTPC3DCylindricalInterpolator();



  // should be in contructor
  fInterpolatorChargeA->SetNR(fNRRows);
  fInterpolatorChargeA->SetNZ(fNZColumns/2);
  fInterpolatorChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorChargeA->SetRList(fLookupRList);
  fInterpolatorChargeA->SetZList(fLookupZListA);
  fInterpolatorChargeA->SetPhiList(fLookupPhiList);
  fInterpolatorChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorChargeC->SetNR(fNRRows);
  fInterpolatorChargeC->SetNZ(fNZColumns/2);
  fInterpolatorChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorChargeC->SetRList(fLookupRList);
  fInterpolatorChargeC->SetZList(fLookupZListC);
  fInterpolatorChargeC->SetPhiList(fLookupPhiList);
  fInterpolatorChargeC->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeA->SetNR(fNRRows);
  fInterpolatorInverseChargeA->SetNZ(fNZColumns/2);
  fInterpolatorInverseChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeA->SetRList(fLookupRList);
  fInterpolatorInverseChargeA->SetZList(fLookupZListA);
  fInterpolatorInverseChargeA->SetPhiList(fLookupPhiList);
  fInterpolatorInverseChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeC->SetNR(fNRRows);
  fInterpolatorInverseChargeC->SetNZ(fNZColumns/2);
  fInterpolatorInverseChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeC->SetRList(fLookupRList);
  fInterpolatorInverseChargeC->SetZList(fLookupZListC);
  fInterpolatorInverseChargeC->SetPhiList(fLookupPhiList);
  fInterpolatorInverseChargeC->SetOrder(fInterpolationOrder);	
  /// local distortion

  fLookupDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListA,fInterpolationOrder);		
	
  fLookupDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListC,fInterpolationOrder);		

  fLookupInverseDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListA,fInterpolationOrder);		
	
  fLookupInverseDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListC,fInterpolationOrder);		

			
}



/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Member values from params
///
/// \param nrrows Int_t number of grid in r direction
/// \param nzcolumns Int_t number of grid in z direction
/// \param nphislices Int_t number of grid in \f$ \phi \f$ direction
/// \param interpolationorder Int_t order of interpolation
/// \param strategy Int_t strategy for global distortion
///
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine
(	
 Int_t nrrows, 
 Int_t nzcolumns, 
 Int_t nphislices,
 Int_t interpolationorder,
 Int_t irregulargridsize,
 Int_t strategy,
 Int_t rbfKernelType
	)
  // : AliTPCCorrection(),fC0(0.),fC1(0.),fCorrectionFactor(1.), fInitLookUp(kFALSE) // removed
  : TNamed("SpaceCharge3D","SpaceCharge3D"),fT1(1.),fT2(1.),fC0(0.),fC1(0.),fCorrectionFactor(1.),  fInitLookUp(kFALSE) // added
{
	
  fInterpolationOrder = interpolationorder; //triliniear 3d //2 there is an error		

  fIrregularGridSize = irregulargridsize; // default 3
  fStrategy = strategy;
  
  fNRRows = nrrows;
  fNPhiSlices = nphislices; // the maximum of phi-slices so far = (8 per sector)
  fNZColumns = nzcolumns; // the maximum on column-slices so  ~ 2cm slicing
  
  fLookupRList = new Double_t[fNRRows];
  fLookupPhiList = new Double_t[fNPhiSlices];
  fLookupZList = new Double_t[fNZColumns];
  
  
  fLookupZListA = new Double_t[fNZColumns/2];
  fLookupZListC = new Double_t[fNZColumns/2];
  
  
  
  Int_t phiSlicesPerSector = fNPhiSlices / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (fNRRows - 1) ;
  const Float_t gridSizeZ   = o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (fNZColumns/2 -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / fNPhiSlices;
	
	
	
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) fLookupPhiList[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < fNRRows ; i++ )  fLookupRList[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < fNZColumns; j++ ) fLookupZList[j]  = (j * gridSizeZ)  - o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; 
	
  for ( Int_t j = 0 ; j < fNZColumns/2; j++ ) {
    fLookupZListA[j]  = (j * gridSizeZ);
    //printf("%d,%f\n",j,fLookupZListA[j]);
  }
	
  // negative 
  for ( Int_t j = 0 ; j < fNZColumns/2; j++ ) fLookupZListC[j]  = -(j * gridSizeZ);
	
	
  
	
  printf("allocating look up %d slices, (%d,%d), interp order : %d\n",fNPhiSlices,fNRRows,fNZColumns,fInterpolationOrder);
  
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    fLookUpIntDistDrEz[k]   =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntDistDphiREz[k] =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntDistDz[k]    =  new TMatrixD(fNRRows,fNZColumns);
    
    
    fLookUpIntDistDrEzA[k]   =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDphiREzA[k] =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDzA[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
  
  
    fLookUpIntDistDrEzC[k]   =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDphiREzC[k] =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDzC[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
  
    fLookUpIntCorrDrEz[k]   =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntCorrDphiREz[k] =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntCorrDz[k]    =  new TMatrixD(fNRRows,fNZColumns);    
    
    // split the side
    
    /// no drift
		
    fLookUpErOverEz[k]    =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpEphiOverEz[k]  =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpDeltaEz[k]     =  new TMatrixD(fNRRows,fNZColumns);

    fLookUpErOverEzA[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpEphiOverEzA[k]  =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpDeltaEzA[k]     =  new TMatrixD(fNRRows,fNZColumns/2);


    /// no drift
    fLookUpErOverEzC[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpEphiOverEzC[k]  =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpDeltaEzC[k]     =  new TMatrixD(fNRRows,fNZColumns/2);
    
    fSCdensityDistribution[k] = new TMatrixF(fNRRows,fNZColumns);


    // Correction with irregular interpolation
    fLookUpIntCorrDrEzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);   		//[kNPhi]  
    fLookUpIntCorrDphiREzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);   //[kNPhi]  
    fLookUpIntCorrDzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fRListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);;
    fPhiListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);; 
    fZListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);;   

    fLookUpIntCorrDrEzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);   		//[kNPhi]  
    fLookUpIntCorrDphiREzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);   //[kNPhi]  
    fLookUpIntCorrDzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fRListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);;
    fPhiListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);; 
    fZListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);;   

    // lookup for charge
    fLookUpChargeA[k]  = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeC[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeInverseA[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeInverseC[k] = new TMatrixD(fNRRows,fNZColumns/2);

  }
  
  fLookupIntDist = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEz,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREz,fLookupPhiList,fNZColumns,fLookUpIntDistDz,fLookupZList,fInterpolationOrder);

  fLookupIntCorr = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntCorrDrEz,fLookupRList,fNPhiSlices,fLookUpIntCorrDphiREz,fLookupPhiList,fNZColumns,fLookUpIntCorrDz,fLookupZList,fInterpolationOrder);
			
  fLookupIntDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEzA,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREzA,fLookupPhiList,fNZColumns/2,fLookUpIntDistDzA,fLookupZListA,fInterpolationOrder);		
	
  fLookupIntDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEzC,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREzC,fLookupPhiList,fNZColumns/2,fLookUpIntDistDzC,fLookupZListC,fInterpolationOrder);		


  fLookupIntENoDrift =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEz,fLookupRList,fNPhiSlices,fLookUpEphiOverEz,fLookupPhiList,fNZColumns,fLookUpDeltaEz,fLookupZList,fInterpolationOrder);		
	
  fLookupIntENoDriftA =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEzA,fLookupRList,fNPhiSlices,fLookUpEphiOverEzA,fLookupPhiList,fNZColumns/2,fLookUpDeltaEzA,fLookupZListA,fInterpolationOrder);		
	
	
  fLookupIntENoDriftC =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEzC,fLookupRList,fNPhiSlices,fLookUpEphiOverEzC,fLookupPhiList,fNZColumns/2,fLookUpDeltaEzC,fLookupZListC,fInterpolationOrder);		
			

  // look up table for correction error based on irregular grid
  fLookupIntCorrIrregularA =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     fNRRows,
     fLookUpIntCorrDrEzIrregularA,
     fRListIrregularA,
     fLookupRList,
     fNPhiSlices,
     fLookUpIntCorrDphiREzIrregularA,
     fPhiListIrregularA,
     fLookupPhiList,
     fNZColumns/2,
     fLookUpIntCorrDzIrregularA,
     fZListIrregularA,
     fLookupZListA,
     2,
     GetIrregularGridSize(),
     GetIrregularGridSize(),
     GetIrregularGridSize(),
     1
     );
		
  fLookupIntCorrIrregularC =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     fNRRows,
     fLookUpIntCorrDrEzIrregularC,
     fRListIrregularC,
     fLookupRList,
     fNPhiSlices,
     fLookUpIntCorrDphiREzIrregularC,
     fPhiListIrregularC,
     fLookupPhiList,
     fNZColumns/2,
     fLookUpIntCorrDzIrregularC,
     fZListIrregularC,
     fLookupZListC,
     2,
     GetIrregularGridSize(),
     GetIrregularGridSize(),
     GetIrregularGridSize(),
     1
     );

  fInterpolatorChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorChargeC = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeC = new AliTPC3DCylindricalInterpolator();

  // should be in contructor
  fInterpolatorChargeA->SetNR(fNRRows);
  fInterpolatorChargeA->SetNZ(fNZColumns/2);
  fInterpolatorChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorChargeA->SetRList(fLookupRList);
  fInterpolatorChargeA->SetZList(fLookupZListA);
  fInterpolatorChargeA->SetPhiList(fLookupPhiList);
  fInterpolatorChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorChargeC->SetNR(fNRRows);
  fInterpolatorChargeC->SetNZ(fNZColumns/2);
  fInterpolatorChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorChargeC->SetRList(fLookupRList);
  fInterpolatorChargeC->SetZList(fLookupZListC);
  fInterpolatorChargeC->SetPhiList(fLookupPhiList);
  fInterpolatorChargeC->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeA->SetNR(fNRRows);
  fInterpolatorInverseChargeA->SetNZ(fNZColumns/2);
  fInterpolatorInverseChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeA->SetRList(fLookupRList);
  fInterpolatorInverseChargeA->SetZList(fLookupZListA);
  fInterpolatorInverseChargeA->SetPhiList(fLookupPhiList);
  fInterpolatorInverseChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeC->SetNR(fNRRows);
  fInterpolatorInverseChargeC->SetNZ(fNZColumns/2);
  fInterpolatorInverseChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeC->SetRList(fLookupRList);
  fInterpolatorInverseChargeC->SetZList(fLookupZListC);
  fInterpolatorInverseChargeC->SetPhiList(fLookupPhiList);
  fInterpolatorInverseChargeC->SetOrder(fInterpolationOrder);

  fLookupDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListA,fInterpolationOrder);		
	
  fLookupDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListA,fInterpolationOrder);		

  fLookupInverseDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListA,fInterpolationOrder);		
	
  fLookupInverseDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,NULL,fLookupRList,fNPhiSlices,NULL,fLookupPhiList,fNZColumns/2,NULL,fLookupZListA,fInterpolationOrder);		

  fRBFKernelType = rbfKernelType;
  fLookupIntCorrIrregularA->SetKernelType(rbfKernelType);
  fLookupIntCorrIrregularC->SetKernelType(rbfKernelType);
}

/// Destruction for AliTPCSpaceCharge3DDriftLine
/// Deallocate memory for lookup table and charge distribution
AliTPCSpaceCharge3DDriftLine::~AliTPCSpaceCharge3DDriftLine() {
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    delete fLookUpIntDistDrEz[k];
    delete fLookUpIntDistDphiREz[k];
    delete fLookUpIntDistDz[k];
    
    delete fLookUpIntDistDrEzA[k];
    delete fLookUpIntDistDphiREzA[k];
    delete fLookUpIntDistDzA[k];
    
    delete fLookUpIntDistDrEzC[k];
    delete fLookUpIntDistDphiREzC[k];
    delete fLookUpIntDistDzC[k];
    
    
    delete fLookUpIntCorrDrEz[k];
    delete fLookUpIntCorrDphiREz[k];
    delete fLookUpIntCorrDz[k];
    
    
    /// no drift
    delete fLookUpErOverEz[k];
    delete fLookUpEphiOverEz[k];
    delete fLookUpDeltaEz[k];

    delete fLookUpErOverEzA[k];
    delete fLookUpEphiOverEzA[k];
    delete fLookUpDeltaEzA[k];


    /// no drift
    delete fLookUpErOverEzC[k];
    delete fLookUpEphiOverEzC[k];
    delete fLookUpDeltaEzC[k];
		
    delete fSCdensityDistribution[k];
    
    delete 		fLookUpIntCorrDrEzIrregularA[k];   		//[kNPhi]  
    delete 		fLookUpIntCorrDphiREzIrregularA[k];   //[kNPhi]  
    delete 		fLookUpIntCorrDzIrregularA[k];
    delete 		fRListIrregularA[k];
    delete 		fPhiListIrregularA[k]; 
    delete 		fZListIrregularA[k];   

    delete 		fLookUpIntCorrDrEzIrregularC[k];   		//[kNPhi]  
    delete 		fLookUpIntCorrDphiREzIrregularC[k];   //[kNPhi]  
    delete 		fLookUpIntCorrDzIrregularC[k];
    delete 		fRListIrregularC[k];
    delete 		fPhiListIrregularC[k] ; 
    delete 		fZListIrregularC[k];   

    // lookup for charge
    delete fLookUpChargeA[k]; 
    delete fLookUpChargeC[k]; 
    delete fLookUpChargeInverseA[k]; 
    delete fLookUpChargeInverseC[k]; 


    
  }
  delete fLookupRList;
  delete fLookupPhiList;
  delete fLookupZList; 
  delete fLookupZListA; 
  delete fLookupZListC;
   
  delete fLookupIntDist;
  delete fLookupIntDistA;
  delete fLookupIntDistC;
  delete fLookupIntENoDriftA;
	
	
  delete fLookupIntENoDrift;
	
  delete fLookupIntENoDriftC;
	
  
  delete fLookupIntCorr;

  delete fLookupIntCorrIrregularA;
		
  delete fLookupIntCorrIrregularC;

  delete fInterpolatorChargeA;
  delete 	fInterpolatorChargeC;
  delete 	fInterpolatorInverseChargeA;
  delete fInterpolatorInverseChargeC;

  delete fLookupDistA;
  delete fLookupDistC;
  delete fLookupInverseDistA;
  delete fLookupInverseDistC;



}


// should be called after construction
void AliTPCSpaceCharge3DDriftLine::Recreate(
					    Int_t nrrows, 
					    Int_t nzcolumns, 
					    Int_t nphislices,
					    Int_t interpolationorder,
					    Int_t strategy
					    ) 
{
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    delete fLookUpIntDistDrEz[k];
    delete fLookUpIntDistDphiREz[k];
    delete fLookUpIntDistDz[k];
    
    delete fLookUpIntDistDrEzA[k];
    delete fLookUpIntDistDphiREzA[k];
    delete fLookUpIntDistDzA[k];
    
    delete fLookUpIntDistDrEzC[k];
    delete fLookUpIntDistDphiREzC[k];
    delete fLookUpIntDistDzC[k];
    
    
    delete fLookUpIntCorrDrEz[k];
    delete fLookUpIntCorrDphiREz[k];
    delete fLookUpIntCorrDz[k];
    
    
    /// no drift
    delete fLookUpErOverEz[k];
    delete fLookUpEphiOverEz[k];
    delete fLookUpDeltaEz[k];

    delete fLookUpErOverEzA[k];
    delete fLookUpEphiOverEzA[k];
    delete fLookUpDeltaEzA[k];


    /// no drift
    delete fLookUpErOverEzC[k];
    delete fLookUpEphiOverEzC[k];
    delete fLookUpDeltaEzC[k];
		
    delete fSCdensityDistribution[k];
    
    delete 		fLookUpIntCorrDrEzIrregularA[k];   		//[kNPhi]  
    delete 		fLookUpIntCorrDphiREzIrregularA[k];   //[kNPhi]  
    delete 		fLookUpIntCorrDzIrregularA[k];
    delete 		fRListIrregularA[k];
    delete 		fPhiListIrregularA[k]; 
    delete 		fZListIrregularA[k];   

    delete 		fLookUpIntCorrDrEzIrregularC[k];   		//[kNPhi]  
    delete 		fLookUpIntCorrDphiREzIrregularC[k];   //[kNPhi]  
    delete 		fLookUpIntCorrDzIrregularC[k];
    delete 		fRListIrregularC[k];
    delete 		fPhiListIrregularC[k] ; 
    delete 		fZListIrregularC[k];   

    // lookup for charge
    delete fLookUpChargeA[k]; 
    delete fLookUpChargeC[k]; 
    delete fLookUpChargeInverseA[k]; 
    delete fLookUpChargeInverseC[k]; 


    
  }
  delete fLookupRList;
  delete fLookupPhiList;
  delete fLookupZList; 
  delete fLookupZListA; 
  delete fLookupZListC;
   
  delete fLookupIntDist;
  delete fLookupIntDistA;
  delete fLookupIntDistC;
  delete fLookupIntENoDriftA;
	
	
  delete fLookupIntENoDrift;
	
  delete fLookupIntENoDriftC;
	
  
  delete fLookupIntCorr;

  delete fLookupIntCorrIrregularA;
		
  delete fLookupIntCorrIrregularC;

  delete fInterpolatorChargeA;
  delete 	fInterpolatorChargeC;
  delete 	fInterpolatorInverseChargeA;
  delete fInterpolatorInverseChargeC;



  // allocate again
  fInterpolationOrder = interpolationorder; //triliniear 3d //2 there is an error		
  fStrategy = strategy;

  fNRRows = nrrows;
  fNPhiSlices = nphislices; // the maximum of phi-slices so far = (8 per sector)
  fNZColumns = nzcolumns; // the maximum on column-slices so  ~ 2cm slicing
  
  fLookupRList = new Double_t[fNRRows];
  fLookupPhiList = new Double_t[fNPhiSlices];
  fLookupZList = new Double_t[fNZColumns];
  
  
  fLookupZListA = new Double_t[fNZColumns/2];
  fLookupZListC = new Double_t[fNZColumns/2];
  
  
  
  Int_t phiSlicesPerSector = fNPhiSlices / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (fNRRows - 1) ;
  const Float_t gridSizeZ   = o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (fNZColumns/2 -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / fNPhiSlices;
	
	
	
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) fLookupPhiList[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < fNRRows ; i++ )  fLookupRList[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < fNZColumns; j++ ) fLookupZList[j]  = (j * gridSizeZ)  - o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; 
	
  for ( Int_t j = 0 ; j < fNZColumns/2; j++ ) {
    fLookupZListA[j]  = (j * gridSizeZ);
    //printf("%d,%f\n",j,fLookupZListA[j]);
  }
	
  // negative 
  for ( Int_t j = 0 ; j < fNZColumns/2; j++ ) fLookupZListC[j]  = -(j * gridSizeZ);
	
	
  
	
  printf("allocating look up %d slices, (%d,%d), interp order : %d\n",fNPhiSlices,fNRRows,fNZColumns,fInterpolationOrder);
  
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    fLookUpIntDistDrEz[k]   =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntDistDphiREz[k] =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntDistDz[k]    =  new TMatrixD(fNRRows,fNZColumns);
    
    
    fLookUpIntDistDrEzA[k]   =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDphiREzA[k] =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDzA[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
  
  
    fLookUpIntDistDrEzC[k]   =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDphiREzC[k] =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpIntDistDzC[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
  
    fLookUpIntCorrDrEz[k]   =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntCorrDphiREz[k] =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpIntCorrDz[k]    =  new TMatrixD(fNRRows,fNZColumns);    
    
    // split the side
    
    /// no drift
		
    fLookUpErOverEz[k]    =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpEphiOverEz[k]  =  new TMatrixD(fNRRows,fNZColumns);
    fLookUpDeltaEz[k]     =  new TMatrixD(fNRRows,fNZColumns);

    fLookUpErOverEzA[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpEphiOverEzA[k]  =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpDeltaEzA[k]     =  new TMatrixD(fNRRows,fNZColumns/2);


    /// no drift
    fLookUpErOverEzC[k]    =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpEphiOverEzC[k]  =  new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpDeltaEzC[k]     =  new TMatrixD(fNRRows,fNZColumns/2);
    
    fSCdensityDistribution[k] = new TMatrixF(fNRRows,fNZColumns);


    // Correction with irregular interpolation
    fLookUpIntCorrDrEzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);   		//[kNPhi]  
    fLookUpIntCorrDphiREzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);   //[kNPhi]  
    fLookUpIntCorrDzIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fRListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);;
    fPhiListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);; 
    fZListIrregularA[k] = new TMatrixD(fNRRows,fNZColumns/2);;   

    fLookUpIntCorrDrEzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);   		//[kNPhi]  
    fLookUpIntCorrDphiREzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);   //[kNPhi]  
    fLookUpIntCorrDzIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fRListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);;
    fPhiListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);; 
    fZListIrregularC[k] = new TMatrixD(fNRRows,fNZColumns/2);;   

    // lookup for charge
    fLookUpChargeA[k]  = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeC[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeInverseA[k] = new TMatrixD(fNRRows,fNZColumns/2);
    fLookUpChargeInverseC[k] = new TMatrixD(fNRRows,fNZColumns/2);

  }
  
  fLookupIntDist = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEz,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREz,fLookupPhiList,fNZColumns,fLookUpIntDistDz,fLookupZList,fInterpolationOrder);

  fLookupIntCorr = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntCorrDrEz,fLookupRList,fNPhiSlices,fLookUpIntCorrDphiREz,fLookupPhiList,fNZColumns,fLookUpIntCorrDz,fLookupZList,fInterpolationOrder);
			
  fLookupIntDistA = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEzA,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREzA,fLookupPhiList,fNZColumns/2,fLookUpIntDistDzA,fLookupZListA,fInterpolationOrder);		
	
  fLookupIntDistC = 
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpIntDistDrEzC,fLookupRList,fNPhiSlices,fLookUpIntDistDphiREzC,fLookupPhiList,fNZColumns/2,fLookUpIntDistDzC,fLookupZListC,fInterpolationOrder);		


  fLookupIntENoDrift =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEz,fLookupRList,fNPhiSlices,fLookUpEphiOverEz,fLookupPhiList,fNZColumns,fLookUpDeltaEz,fLookupZList,fInterpolationOrder);		
	
  fLookupIntENoDriftA =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEzA,fLookupRList,fNPhiSlices,fLookUpEphiOverEzA,fLookupPhiList,fNZColumns/2,fLookUpDeltaEzA,fLookupZListA,fInterpolationOrder);		
	
	
  fLookupIntENoDriftC =
    new AliTPCLookUpTable3DInterpolatorD(
					 fNRRows,fLookUpErOverEzC,fLookupRList,fNPhiSlices,fLookUpEphiOverEzC,fLookupPhiList,fNZColumns/2,fLookUpDeltaEzC,fLookupZListC,fInterpolationOrder);		
			

  // look up table for correction error based on irregular grid
  fLookupIntCorrIrregularA =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     fNRRows,
     fLookUpIntCorrDrEzIrregularA,
     fRListIrregularA,
     fLookupRList,
     fNPhiSlices,
     fLookUpIntCorrDphiREzIrregularA,
     fPhiListIrregularA,
     fLookupPhiList,
     fNZColumns/2,
     fLookUpIntCorrDzIrregularA,
     fZListIrregularA,
     fLookupZListA,
     2,
     3,
     3,
     3,
     1
     );
		
  fLookupIntCorrIrregularC =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     fNRRows,
     fLookUpIntCorrDrEzIrregularC,
     fRListIrregularC,
     fLookupRList,
     fNPhiSlices,
     fLookUpIntCorrDphiREzIrregularC,
     fPhiListIrregularC,
     fLookupPhiList,
     fNZColumns/2,
     fLookUpIntCorrDzIrregularC,
     fZListIrregularC,
     fLookupZListC,
     2,
     3,
     3,
     3,
     1
     );

  fInterpolatorChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorChargeC = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeC = new AliTPC3DCylindricalInterpolator();

  // should be in contructor
  fInterpolatorChargeA->SetNR(fNRRows);
  fInterpolatorChargeA->SetNZ(fNZColumns/2);
  fInterpolatorChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorChargeA->SetRList(fLookupRList);
  fInterpolatorChargeA->SetZList(fLookupZListA);
  fInterpolatorChargeA->SetPhiList(fLookupPhiList);
  fInterpolatorChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorChargeC->SetNR(fNRRows);
  fInterpolatorChargeC->SetNZ(fNZColumns/2);
  fInterpolatorChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorChargeC->SetRList(fLookupRList);
  fInterpolatorChargeC->SetZList(fLookupZListC);
  fInterpolatorChargeC->SetPhiList(fLookupPhiList);
  fInterpolatorChargeC->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeA->SetNR(fNRRows);
  fInterpolatorInverseChargeA->SetNZ(fNZColumns/2);
  fInterpolatorInverseChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeA->SetRList(fLookupRList);
  fInterpolatorInverseChargeA->SetZList(fLookupZListA);
  fInterpolatorInverseChargeA->SetPhiList(fLookupPhiList);
  fInterpolatorInverseChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeC->SetNR(fNRRows);
  fInterpolatorInverseChargeC->SetNZ(fNZColumns/2);
  fInterpolatorInverseChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeC->SetRList(fLookupRList);
  fInterpolatorInverseChargeC->SetZList(fLookupZListC);
  fInterpolatorInverseChargeC->SetPhiList(fLookupPhiList);
  fInterpolatorInverseChargeC->SetOrder(fInterpolationOrder);

}



/// Creating look-up tables of Correction/Distortion by integration following 
/// drift line, input from space charge 3d histogram (fSpaceCharge3D) and boundary values are filled with zeroes
///
/// TODO: provide an interface for setting boundary values
/// 
/// The algorithm and implementations of this function is the following:
///
/// Do for each side A,C
///
/// 1) Solving \f$ \nabla^2 \Phi(r,\phi,z) = -  \rho(r,\phi,z)\f$
/// ~~~ Calling poisson solver
/// fPoissonSolver->PoissonSolver3D( matricesV, matricesCharge, rrow, zcolumn, phiSlice, maxIteration, symmetry ) ;
/// ~~~
///
/// 2) Get the electric field \f$ \vec{E} = - \nabla \Phi(r,\phi,z) \f$		
/// ~~~ 
/// ElectricField( matricesV, matricesEr,  matricesEphi, matricesEz, rrow, zcolumn, phiSlice, 
/// gridSizeR, gridSizePhi ,gridSizeZ,symmetry, o2::TPC::AliTPCPoissonSolver::fgkIFCRadius); 
/// ~~~
///
/// 3) Calculate local distortion and correction, useing Langevin formula
/// ~~~ cxx
/// LocalDistCorrDz (matricesEr, matricesEphi, 	matricesEz,  
///	matricesDistDrDz,  matricesDistDphiRDz, matricesDistDz,
///	matricesCorrDrDz,  matricesCorrDphiRDz, matricesCorrDz,
///	rrow,  zcolumn, phiSlice, gridSizeZ, ezField);	
/// ~~~
///
/// 4) Integrate distortion by following the drift line
///
/// 5) Fill look up table for Integral distortion 	
///		
/// 6) Fill look up table for Integral correction 	
///
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
/// \param maxIteration Int_t Maximum iteration for poisson solver
/// \param stoppingConv Convergence error stopping conditioin for poisson solver
///
/// \post Lookup tables for distortion: 
/// ~~~ 
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz  
/// ~~~ 
/// and correction: 
/// ~~~ 
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz  
/// ~~~ 
/// are initialized 
/// 
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoissonIntegralDz
(
 Int_t rrow, 
 Int_t zcolumn, 
 Int_t phiSlice,
 Int_t maxIteration,
 Double_t stoppingConv
 )
{
  Int_t phiSlicesPerSector = phiSlice / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	

	
  // local variables
  Float_t radius0, phi0, z0;
		
  // memory allocation for temporary matrices:
  // potential (boundary values), charge distribution
  TMatrixD *matricesV[phiSlice], *matricesCharge[phiSlice];
  TMatrixD *matricesEr[phiSlice], *matricesEphi[phiSlice],  *matricesEz[phiSlice];
  TMatrixD *matricesDistDrDz[phiSlice], *matricesDistDphiRDz[phiSlice],  *matricesDistDz[phiSlice];
  TMatrixD *matricesCorrDrDz[phiSlice], *matricesCorrDphiRDz[phiSlice],  *matricesCorrDz[phiSlice];
  TMatrixD *matricesGDistDrDz[phiSlice], *matricesGDistDphiRDz[phiSlice],  *matricesGDistDz[phiSlice];
  TMatrixD *matricesGCorrDrDz[phiSlice], *matricesGCorrDphiRDz[phiSlice],  *matricesGCorrDz[phiSlice];
	
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    matricesV[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesCharge[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesEr[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesEphi[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesEz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDz[k]     	=   new TMatrixD(rrow,zcolumn) ;		
    matricesCorrDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesCorrDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesCorrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesGDistDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesGDistDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesGDistDz[k]     	=   new TMatrixD(rrow,zcolumn) ;		
    matricesGCorrDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesGCorrDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesGCorrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
  }
	
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t  rlist[rrow], zedlist[zcolumn] , philist[phiSlice];
	
	
		
  TStopwatch w;

		
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) philist[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow ; i++ )  rlist[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn ; j++ ) zedlist[j]  = j * gridSizeZ ;
	
	

  AliTPCLookUpTable3DInterpolatorD *lookupLocalDist = 
    new AliTPCLookUpTable3DInterpolatorD(
					 rrow,matricesDistDrDz,rlist,phiSlice,matricesDistDphiRDz,philist,zcolumn,matricesDistDz,zedlist,fInterpolationOrder);
			
  AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr = 
    new AliTPCLookUpTable3DInterpolatorD(
					 rrow,matricesCorrDrDz,rlist,phiSlice,matricesCorrDphiRDz,philist,zcolumn,matricesCorrDz,zedlist,fInterpolationOrder);
			
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalDist = 
    new AliTPCLookUpTable3DInterpolatorD(
					 rrow,matricesGDistDrDz,rlist,phiSlice,matricesGDistDphiRDz,philist,zcolumn,matricesGDistDz,zedlist,fInterpolationOrder);
			
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalCorr = 
    new AliTPCLookUpTable3DInterpolatorD(
					 rrow,matricesGCorrDrDz,rlist,phiSlice,matricesGCorrDphiRDz,philist,zcolumn,matricesGCorrDz,zedlist,fInterpolationOrder);
			
	
	
  // should be set, in another place
  const Int_t   symmetry = 0;

  // for irregular
  TMatrixD** matricesIrregularDrDz = NULL;
  TMatrixD** matricesIrregularDphiRDz = NULL;
  TMatrixD** matricesIrregularDz = NULL;
  TMatrixD** matricesPhiIrregular = NULL;
  TMatrixD** matricesRIrregular = NULL;
  TMatrixD** matricesZIrregular = NULL;

  // for charge
  TMatrixD** matricesLookUpCharge = NULL;
  AliTPC3DCylindricalInterpolator *chargeInterpolator;
	
	
  // do if look up table haven't be initialized
  if ( !fInitLookUp ) {
    //* do for 2 sides (z direction)
    // solve Poisson3D twice; once for +Z and once for -Z	
		
    for ( Int_t side = 0 ; side < 2 ; side++ ) {  					
      // zeroing global distortion/correction
      for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
	matricesGDistDrDz[k]->Zero();
	matricesGDistDphiRDz[k]->Zero();
	matricesGDistDz[k]->Zero();		
	matricesGCorrDrDz[k]->Zero();
	matricesGCorrDphiRDz[k]->Zero();
	matricesGCorrDz[k]->Zero();						

	if (side == 0)
	  {

								

	    matricesIrregularDrDz = fLookUpIntCorrDrEzIrregularA;
	    matricesIrregularDphiRDz = fLookUpIntCorrDphiREzIrregularA;
	    matricesIrregularDz = fLookUpIntCorrDzIrregularA;
	    matricesPhiIrregular = fPhiListIrregularA;
	    matricesRIrregular = fRListIrregularA;
	    matricesZIrregular = fZListIrregularA;

	    matricesLookUpCharge = fLookUpChargeA;
	    chargeInterpolator = fInterpolatorChargeA;


	    fLookupDistA->SetLookUpR(matricesDistDrDz);
	    fLookupDistA->SetLookUpPhi(matricesDistDphiRDz);
	    fLookupDistA->SetLookUpZ(matricesDistDz);

					
	  }	else
	  {
	    matricesIrregularDrDz = fLookUpIntCorrDrEzIrregularC;
	    matricesIrregularDphiRDz = fLookUpIntCorrDphiREzIrregularC;
	    matricesIrregularDz = fLookUpIntCorrDzIrregularC;

	    matricesPhiIrregular = fPhiListIrregularC;
	    matricesRIrregular = fRListIrregularC;
	    matricesZIrregular = fZListIrregularC;

	    matricesLookUpCharge = fLookUpChargeC;
	    chargeInterpolator = fInterpolatorChargeC;

	    fLookupDistC->SetLookUpR(matricesDistDrDz);
	    fLookupDistC->SetLookUpPhi(matricesDistDphiRDz);
	    fLookupDistC->SetLookUpZ(matricesDistDz);

	  }

      }
			
      for ( Int_t k = 0 ; k < phiSlice ; k++ )  {				
	TMatrixD *mV    =  matricesV[k] ;
	TMatrixD *mCharge    =  matricesCharge[k] ;				
	TMatrixD *mLookUpCharge    =  matricesLookUpCharge[k] ;				
	phi0    = philist[k];	
	for ( Int_t i = 0 ; i < rrow ; i++ ) {
	  radius0 = rlist[i];					
	  for ( Int_t j = 0 ; j < zcolumn ; j++ ) {						
	    z0 = zedlist[j];						
	    if (side==1) z0= -TMath::Abs(zedlist[j]);							
	    if (fSpaceChargeHistogram3D != NULL) {
	      // * Boundary values and charge distribution setup
	      (*mV)(i,j) = 0.0 ;
	      (*mCharge)(i,j) = -1 * InterpolatePhi(fSpaceChargeHistogram3D,phi0,radius0,z0);

	      // Ernst: modified to remove index outside of range error
	      // (*mLookUpCharge)(i,j) = (*mCharge)(i,j);
	      if (j<(zcolumn/2)){
		z0 = zedlist[j*2];						
		if (side==1) z0= -TMath::Abs(zedlist[j*2]);
		(*mLookUpCharge)(i,j) = -1 * InterpolatePhi(fSpaceChargeHistogram3D,phi0,radius0,z0);
	      }
	    }
	  }
	}				
      }	
			
			
      // setvals for lookup charge
      w.Start();

      chargeInterpolator->SetVals(matricesLookUpCharge);
      chargeInterpolator->InitCubicSpline();
      w.Stop();

      // AliInfo(Form("Step 0: Preparing Charge interpolator: %f\n",w.CpuTime())); // removed
      LOG(INFO) << Form("Step 0: Preparing Charge interpolator: %f\n",w.CpuTime()) << FairLogger::endl;
      //
			
      o2::TPC::AliTPCPoissonSolver::fgConvErr = stoppingConv;
      fPoissonSolver->SetStrategy(o2::TPC::AliTPCPoissonSolver::kMultigrid);
      (fPoissonSolver->fMgParameters).cycleType = o2::TPC::AliTPCPoissonSolver::kFCycle;	
      (fPoissonSolver->fMgParameters).isFull3D  = kFALSE;
      (fPoissonSolver->fMgParameters).NMGCYCLE  = maxIteration;
      (fPoissonSolver->fMgParameters).MAXLOOP   = 6;
      w.Start();
      fPoissonSolver->PoissonSolver3D( matricesV, matricesCharge, rrow, zcolumn, phiSlice, maxIteration, symmetry ) ;
      w.Stop();
      printf("%f\t",w.CpuTime());	
      
			
      // AliInfo(Form("Step 1: Poisson solver: %f\n",w.CpuTime())); // removed
      LOG(INFO) << Form("Step 1: Poisson solver: %f\n",w.CpuTime()) << FairLogger::endl;
      w.Start();
      ElectricField( matricesV, 
		     matricesEr,  matricesEphi, matricesEz, rrow, zcolumn, phiSlice, 
		     gridSizeR, gridSizePhi ,gridSizeZ,symmetry, o2::TPC::AliTPCPoissonSolver::fgkIFCRadius); 
      w.Stop();
      // AliInfo(Form("Step 2: Electric Field Calculation: %f\n",w.CpuTime())); // removed
      LOG(INFO) << Form("Step 2: Electric Field Calculation: %f\n",w.CpuTime()) << FairLogger::endl;

			
      //		AliInfo("Step 3: Calculate local distortion");
			
      w.Start();
			
      LocalDistCorrDz (matricesEr, matricesEphi, 	matricesEz,  
		       matricesDistDrDz,  matricesDistDphiRDz, matricesDistDz,
		       matricesCorrDrDz,  matricesCorrDphiRDz, matricesCorrDz,
		       rrow,  zcolumn, phiSlice, gridSizeZ, ezField);	
			
      w.Stop();
			
			
      //// copy to 1D interpolator /////
      lookupLocalDist->CopyVals();
      lookupLocalCorr->CopyVals();

      if (side == 0) 
	fLookupDistA->CopyVals();
      else
	fLookupDistC->CopyVals();

      //// 
			
      // AliInfo(Form("Step 3: Local distortion and correction: %f\n",w.CpuTime())); //  removed
      LOG(INFO) << Form("Step 3: Local distortion and correction: %f\n",w.CpuTime()) << FairLogger::endl;

			
      w.Start();
		
		
      //		AliInfo("Step 4a: Integrate distortion and correction the drift line");
			
      if (fStrategy == kNaive) 
	IntegrateDistCorrDriftLineDz(
				     lookupLocalDist, 
				     matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
				     lookupLocalCorr, 
				     matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  
				     matricesIrregularDrDz, matricesIrregularDphiRDz,  matricesIrregularDz,
				     matricesRIrregular, matricesPhiIrregular,matricesZIrregular,	
				     rrow,  zcolumn, phiSlice, rlist, philist, zedlist
				     );  
      else
	IntegrateDistCorrDriftLineDzOpt2(
					 lookupLocalDist, 
					 matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
					 lookupLocalCorr, 
					 matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  
					 rrow,  zcolumn, phiSlice, rlist, philist, zedlist
					 );  
			
			
			
      w.Stop();
      // AliInfo(Form("Step 4: Global correction/distortion: %f\n",w.CpuTime())); // removed
      LOG(INFO) << Form("Step 4: Global correction/distortion: %f\n",w.CpuTime()) << FairLogger::endl;
      w.Start();
			
			
      //// copy to 1D interpolator /////
      lookupGlobalDist->CopyVals();
      lookupGlobalCorr->CopyVals();
      //// 
			
      FillLookUpTable(lookupGlobalDist, 
      	fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz,
      	rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  
		
			
      FillLookUpTable(lookupGlobalCorr, 
		      fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz,
		      rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  
			
      w.Stop();
      // AliInfo(Form("Step 5: Filling up the look up: %f\n",w.CpuTime())); // removed
      LOG(INFO) << Form("Step 5: Filling up the look up: %f\n",w.CpuTime()) << FairLogger::endl;

		
		
		
      if ( side == 0 ) {
	FillLookUpTableA(lookupGlobalDist,
			 fLookUpIntDistDrEzA,fLookUpIntDistDphiREzA,fLookUpIntDistDzA,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntDistA->CopyVals();		
				
	fLookupIntCorrIrregularA->CopyVals(); 


	// AliInfo(" A side done"); // removed
	LOG(INFO) << " A side done" << FairLogger::endl;
      }
      if ( side == 1 ) {
	FillLookUpTableC(lookupGlobalDist,
			 fLookUpIntDistDrEzC,fLookUpIntDistDphiREzC,fLookUpIntDistDzC,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntDistC->CopyVals();				
				

	fLookupIntCorrIrregularC->CopyVals(); 
	// AliInfo(" C side done"); // removed
	LOG(INFO) << " C side done" << FairLogger::endl;
      }
			
      //if ( side == 0 ) AliInfo(" A side done");
      //if ( side == 1 ) AliInfo(" C side done");		
    }
	
    //// copy to 1D interpolator /////
    fLookupIntDist->CopyVals();
    fLookupIntCorr->CopyVals();
    //// 
		
    fInitLookUp = kTRUE;
  }

 
 
  // memory deallocation for temporary matrices
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    delete matricesV[k];
    delete matricesCharge[k];		
    delete matricesEr[k];		
    delete matricesEphi[k];		
    delete matricesEz[k];
    delete matricesDistDrDz[k] ;
    delete matricesDistDphiRDz[k]	;
    delete matricesDistDz[k] ;
		
    delete matricesCorrDrDz[k] ;
    delete matricesCorrDphiRDz[k]	;
    delete matricesCorrDz[k] ;
    delete matricesGDistDrDz[k] ;
    delete matricesGDistDphiRDz[k]	;
    delete matricesGDistDz[k] ;
		
    delete matricesGCorrDrDz[k] ;
    delete matricesGCorrDphiRDz[k]	;
    delete matricesGCorrDz[k] ;
		
  }
  delete lookupLocalDist;
  delete lookupLocalCorr;
  delete lookupGlobalDist;
  delete lookupGlobalCorr;
}



/// Creating look-up tables of Correction/Distortion by integration following 
/// drift line with known distributions for potential and spacecharge.
/// 
/// \param  matricesVA TMatrixD** potential distribution in side A (output)
/// \param  matricesChargeA TMatrixD** charge distribution in side A (input)
/// \param  matricesVC  TMatrixD** potential distribution in side C (output)
/// \param  matricesChargeC TMatrixD** charge distribution in side C (input)
/// \param rrow	Int_t  number of grid in row direction
///	\param zcolumn Int_t number of grid in z direction
/// \param phiSlice 	Int_t number of slicees in phi direction
/// \param maxIteration Int_t max iteration for convergence
/// \param stoppingConv Double_t stopping criteria for convergence
/// \param matricesErA TMatrixD** electric field r direction distribution side A
/// \param matricesEphiA TMatrixD** electric field phi direction distribution side A
/// \param matricesEzA TMatrixD** electric field z direction distribution side A
/// \param matricesErC TMatrixD** electric field r distribution side C
/// \param matricesEphiC TMatrixD** electric field phi distribution side C
/// \param matricesEzC TMatrixD** electric field z distribution side C
/// \param matricesDistDrDzA TMatrixD**  local r distortion (output) A side
/// \param matricesDistDphiRDzA TMatrixD** local r phi distortion (output) A side
/// \param matricesDistDzA TMatrixD**  local z distortion (output) A side
/// \param matricesCorrDrDzA TMatrixD** local r correction (output) A side
/// \param matricesCorrDphiRDzA TMatrixD** local r phi correction (output) A side
/// \param matricesCorrDzA 	TMatrixD** local z correction (output) A side
/// \param matricesDistDrDzC 	TMatrixD**   local r distortion (output) C side
/// \param matricesDistDphiRDzC 	TMatrixD**  local r phi distortion (output) C side
/// \param matricesDistDzC TMatrixD** local z distortion (output) C side
/// \param matricesCorrDrDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDphiRDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDzC	TMatrixD** local z correction (output) C side
///
/// \post Lookup tables for distortion: 
/// ~~~ 
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz  
/// ~~~ fo
/// and correction: 
/// ~~~ 
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz  
/// ~~~ 
/// are initialized 
/// 
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoisson
(
 Int_t rrow, 
 Int_t zcolumn, 
 Int_t phiSlice,
 Int_t maxIteration,
 Double_t stoppingConv
 )
{
  // Compute grid size for all direction 
  Int_t phiSlicesPerSector = phiSlice / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	
  printf("gridSizeZ in init: %f\n",gridSizeZ);
	
  // local variables
  Float_t radius0, phi0, z0;
	
	
  // memory allocation for temporary matrices:
  // potential (boundary values), charge distribution
	
  TMatrixD *matricesV[phiSlice], *matricesCharge[phiSlice];
  TMatrixD *matricesEr[phiSlice], *matricesEphi[phiSlice],  *matricesEz[phiSlice];
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    matricesEr[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesEphi[k]			=   new TMatrixD(rrow,zcolumn) ;
    matricesEz[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    //matricesErC[k]     	=   new TMatrixD(rrow,zcolumn) ;
    //matricesEphiC[k]	=   new TMatrixD(rrow,zcolumn) ;
    //matricesEzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesV[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesCharge[k]	=   new TMatrixD(rrow,zcolumn) ;
		
		
  }
	
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t  rlist[rrow], zedlist[zcolumn] , philist[phiSlice];
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) philist[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow ; i++ )  rlist[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn ; j++ ) zedlist[j]  = j * gridSizeZ ;
	
	
	
	
  // should be set, in another place
  const Int_t   symmetry = 0;
	
  // do if look up table haven't be initialized
  if ( !fInitLookUp ) {
    //* do for 2 sides (z direction)
			
		
    for ( Int_t side = 0 ; side < 2 ; side++ ) {  				
      /**
	 if (side == 0) {
	 matricesV = matricesVA;
	 matricesCharge = matricesChargeA;
	 matricesEr = matricesErA;
	 matricesEphi = matricesEphiA;
	 matricesEz = matricesEzA;
		
				
	 } else {
	 matricesV = matricesVC;
	 matricesCharge = matricesChargeC;
	 matricesEr = matricesErC;
	 matricesEphi = matricesEphiC;
	 matricesEz = matricesEzC;
		
	 } 
      **/
			
      for ( Int_t k = 0 ; k < phiSlice ; k++ )  {				
	TMatrixD *mV    =  matricesV[k] ;
	TMatrixD *mCharge    =  matricesCharge[k] ;				
	phi0    = philist[k];	
	for ( Int_t i = 0 ; i < rrow ; i++ ) {
	  radius0 = rlist[i];					
	  for ( Int_t j = 0 ; j < zcolumn ; j++ ) {						
	    z0 = zedlist[j];						
	    if (side==1) z0= -TMath::Abs(zedlist[j]);							
	    if (fSpaceChargeHistogram3D != NULL) {
	      // * Boundary values and charge distribution setup
	      (*mV)(i,j) = 0.0 ;
	      (*mCharge)(i,j) = -1 * InterpolatePhi(fSpaceChargeHistogram3D,phi0,radius0,z0);
	    }
	  }
	}				
      }
			
				
      AliTPCLookUpTable3DInterpolatorD *lookupEField = 
	new AliTPCLookUpTable3DInterpolatorD(
					     rrow,
					     matricesEr,
					     rlist,phiSlice,
					     matricesEphi,
					     philist,zcolumn,
					     matricesEz,
					     zedlist,
					     fInterpolationOrder
					     );
			
			
	
      // AliInfo("Step 1: Solving poisson solver");
      LOG(INFO) << "Step 1: Solving poisson solver" << FairLogger::endl;

			
			
			
	
      fPoissonSolver->PoissonSolver3D( matricesV, matricesCharge, rrow, zcolumn, phiSlice, maxIteration, symmetry ) ;
	
	
		
      // AliInfo("Step 2: Calculate electric field");
      LOG(INFO) << "Step 2: Calculate electric field" << FairLogger::endl;
      CalculateEField(
		      matricesV, 
		      matricesEr,  
		      matricesEphi, 
		      matricesEz, 
		      rrow, 
		      zcolumn, 
		      phiSlice, 
		      maxIteration,   
		      symmetry 
		      );
		
      lookupEField->CopyVals();		
      // AliInfo("Step 3: Fill the ");
      LOG(INFO) << "Step 3: Fill the " << FairLogger::endl;

      FillLookUpTable(lookupEField, 
		      fLookUpErOverEz,fLookUpEphiOverEz,fLookUpDeltaEz,
		      rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  
			
			
      if ( side == 0 ) {
	FillLookUpTableA(lookupEField,
			 fLookUpErOverEzA,fLookUpEphiOverEzA,fLookUpDeltaEzA,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
				
	fLookupIntENoDriftA->CopyVals();		
	// AliInfo(" A side done");
	LOG(INFO) << " A side done" << FairLogger::endl;
      }
      if ( side == 1 ) {
	FillLookUpTableC(lookupEField,
			 fLookUpErOverEzC,fLookUpEphiOverEzC,fLookUpDeltaEzC,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntENoDriftC->CopyVals();		
				
	// AliInfo(" C side done");		
	LOG(INFO) << " C side done" << FairLogger::endl;
      }
			
			
      delete lookupEField;

			
    }
		
    fInitLookUp = kTRUE;
		
    fLookupIntENoDrift->CopyVals();		
		
  }
 
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    delete matricesV[k];
    delete matricesCharge[k];		
    delete matricesEr[k];		
    delete matricesEphi[k];		
    delete matricesEz[k];		
  }
	
	
}
/// Creating look-up tables of Correction/Distortion by integration following 
/// drift line with known distributions for potential and spacecharge.
/// 
///
/// \param rrow	Int_t  number of grid in row direction
///	\param zcolumn Int_t number of grid in z direction
/// \param phiSlice 	Int_t number of slicees in phi direction
/// \param maxIteration Int_t max iteration for convergence
/// \param stoppingConv Double_t stopping criteria for convergence
/// \param matricesDistDrDzA TMatrixD**  local r distortion (output) A side
/// \param matricesDistDphiRDzA TMatrixD** local r phi distortion (output) A side
/// \param matricesDistDzA TMatrixD**  local z distortion (output) A side
/// \param matricesCorrDrDzA TMatrixD** local r correction (output) A side
/// \param matricesCorrDphiRDzA TMatrixD** local r phi correction (output) A side
/// \param matricesCorrDzA 	TMatrixD** local z correction (output) A side
/// \param matricesDistDrDzC 	TMatrixD**   local r distortion (output) C side
/// \param matricesDistDphiRDzC 	TMatrixD**  local r phi distortion (output) C side
/// \param matricesDistDzC TMatrixD** local z distortion (output) C side
/// \param matricesCorrDrDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDphiRDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDzC	TMatrixD** local z correction (output) C side
///
/// \post Lookup tables for distortion: 
/// ~~~ 
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz  
/// ~~~ fo
/// and correction: 
/// ~~~ 
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz  
/// ~~~ 
/// are initialized 
/// 
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoissonIntegralDz
(
 Int_t rrow, 
 Int_t zcolumn, 
 Int_t phiSlice,
 Int_t maxIteration,
 Double_t stoppingConv,
 TMatrixD** matricesDistDrDzA,
 TMatrixD** matricesDistDphiRDzA,
 TMatrixD** matricesDistDzA,
 TMatrixD** matricesCorrDrDzA,
 TMatrixD** matricesCorrDphiRDzA,
 TMatrixD** matricesCorrDzA,
 TMatrixD** matricesDistDrDzC,
 TMatrixD** matricesDistDphiRDzC,
 TMatrixD** matricesDistDzC,
 TMatrixD** matricesCorrDrDzC,
 TMatrixD** matricesCorrDphiRDzC,
 TMatrixD** matricesCorrDzC,
 TFormula * intErDzTestFunction,
 TFormula * intEphiRDzTestFunction,
 TFormula * intDzTestFunction
 )
{
  // Compute grid size for all direction 
  Int_t phiSlicesPerSector = phiSlice / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	
	
  // local variables
  Float_t radius0, phi0, z0;	
	
  TMatrixD *matricesGDistDrDz[phiSlice], *matricesGDistDphiRDz[phiSlice],  *matricesGDistDz[phiSlice];
  TMatrixD *matricesGCorrDrDz[phiSlice], *matricesGCorrDphiRDz[phiSlice],  *matricesGCorrDz[phiSlice];
	
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
		
    matricesGDistDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesGDistDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesGDistDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesGCorrDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesGCorrDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesGCorrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
  }
	
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t  rlist[rrow], zedlist[zcolumn] , philist[phiSlice];
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) philist[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow ; i++ )  rlist[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn ; j++ ) zedlist[j]  = j * gridSizeZ ;
  TMatrixD** matricesDistDrDz = NULL;
  TMatrixD** matricesDistDphiRDz = NULL;
  TMatrixD** matricesDistDz = NULL;
  TMatrixD** matricesCorrDrDz = NULL;
  TMatrixD** matricesCorrDphiRDz = NULL;
  TMatrixD** matricesCorrDz = NULL;
	
	
  AliTPCLookUpTable3DInterpolatorD *lookupLocalDist;
			
  AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr; 
					
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalDist = 
    new AliTPCLookUpTable3DInterpolatorD(
					 rrow,matricesGDistDrDz,rlist,phiSlice,matricesGDistDphiRDz,philist,zcolumn,matricesGDistDz,zedlist,fInterpolationOrder);
			
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalCorr = 
    new AliTPCLookUpTable3DInterpolatorD(
					 rrow,matricesGCorrDrDz,rlist,phiSlice,matricesGCorrDphiRDz,philist,zcolumn,matricesGCorrDz,zedlist,fInterpolationOrder);
	
	
  // for irregular
  TMatrixD** matricesIrregularDrDz = NULL;
  TMatrixD** matricesIrregularDphiRDz = NULL;
  TMatrixD** matricesIrregularDz = NULL;
  TMatrixD** matricesPhiIrregular = NULL;
  TMatrixD** matricesRIrregular = NULL;
  TMatrixD** matricesZIrregular = NULL;
	
  // should be set, in another place
  const Int_t   symmetry = 0;
	
  // do if look up table haven't be initialized
  if ( !fInitLookUp ) {
    for ( Int_t side = 0 ; side < 2 ; side++ ) {  				
      if (side == 0) {
	matricesDistDrDz = matricesDistDrDzA;
	matricesDistDphiRDz = matricesDistDphiRDzA;
	matricesDistDz =matricesDistDzA;
	matricesCorrDrDz = matricesCorrDrDzA;
	matricesCorrDphiRDz =matricesCorrDphiRDzA;
	matricesCorrDz =matricesCorrDzA;			



	matricesIrregularDrDz = fLookUpIntCorrDrEzIrregularA;
	matricesIrregularDphiRDz = fLookUpIntCorrDphiREzIrregularA;
	matricesIrregularDz = fLookUpIntCorrDzIrregularA;

	matricesPhiIrregular = fPhiListIrregularA;
	matricesRIrregular = fRListIrregularA;
	matricesZIrregular = fZListIrregularA;


	fLookupDistA->SetLookUpR(matricesDistDrDz);
	fLookupDistA->SetLookUpPhi(matricesDistDphiRDz);
	fLookupDistA->SetLookUpZ(matricesDistDz);


				
      } else {
	matricesDistDrDz = matricesDistDrDzC;
	matricesDistDphiRDz = matricesDistDphiRDzC;
	matricesDistDz =matricesDistDzC;
	matricesCorrDrDz = matricesCorrDrDzC;
	matricesCorrDphiRDz =matricesCorrDphiRDzC;
	matricesCorrDz =matricesCorrDzC;			



	matricesIrregularDrDz = fLookUpIntCorrDrEzIrregularC;
	matricesIrregularDphiRDz = fLookUpIntCorrDphiREzIrregularC;
	matricesIrregularDz = fLookUpIntCorrDzIrregularC;

	matricesPhiIrregular = fPhiListIrregularC;
	matricesRIrregular = fRListIrregularC;
	matricesZIrregular = fZListIrregularC;




	fLookupDistC->SetLookUpR(matricesDistDrDz);
	fLookupDistC->SetLookUpPhi(matricesDistDphiRDz);
	fLookupDistC->SetLookUpZ(matricesDistDz);
      } 
			
      lookupLocalDist = new AliTPCLookUpTable3DInterpolatorD(
							     rrow,matricesDistDrDz,rlist,phiSlice,matricesDistDphiRDz,philist,zcolumn,matricesDistDz,zedlist,fInterpolationOrder);
      lookupLocalCorr = new AliTPCLookUpTable3DInterpolatorD(
							     rrow,matricesCorrDrDz,rlist,phiSlice,matricesCorrDphiRDz,philist,zcolumn,matricesCorrDz,zedlist,fInterpolationOrder);
			
      //// copy to 1D interpolator /////
      lookupLocalDist->CopyVals();
      lookupLocalCorr->CopyVals();


      if (side == 0) 
	fLookupDistA->CopyVals();
      else
	fLookupDistC->CopyVals();

      //// 
			
			
      // AliInfo("Step 4a: Integrate distortion and correction the drift line -- testCorrectness");
      LOG(INFO) <<"Step 4a: Integrate distortion and correction the drift line -- testCorrectness"  << FairLogger::endl;

      if (fStrategy == kNaive) 
	//IntegrateDistCorrDriftLineDz(
	//	lookupLocalDist, 
	//	matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
	//	lookupLocalCorr, 
	//	matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  
	//	rrow,  zcolumn, phiSlice, rlist, philist, zedlist
	//);  
	IntegrateDistCorrDriftLineDz(
				     //lookupLocalDist, 
				     intErDzTestFunction,
				     intEphiRDzTestFunction,
				     intDzTestFunction,
				     ezField,
				     matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
				     //lookupLocalCorr, 
				     matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  


				     matricesIrregularDrDz, matricesIrregularDphiRDz,  matricesIrregularDz,
				     matricesRIrregular, matricesPhiIrregular,matricesZIrregular,	


				     rrow,  zcolumn, phiSlice, rlist, philist, zedlist
				     );  
				
      else
	IntegrateDistCorrDriftLineDzOpt2(matricesDistDrDz,  matricesDistDphiRDz, matricesDistDz, 
					 matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
					 matricesCorrDrDz,  matricesCorrDphiRDz, matricesCorrDz, 
					 matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  
					 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  
			
	
      //// copy to 1D interpolator /////
      lookupGlobalDist->CopyVals();
      lookupGlobalCorr->CopyVals();
      //// 
			
			
			
      // AliInfo("Step 5: Fill look up table for distortion");	
      LOG(INFO) << "Step 5: Fill look up table for distortion" << FairLogger::endl;

      //FillLookUpTable(lookupGlobalDist, 
      //	fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz,
      //	rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  


      if ( side == 0 ) {
	FillLookUpTableA(lookupGlobalDist,
			 fLookUpIntDistDrEzA,fLookUpIntDistDphiREzA,fLookUpIntDistDzA,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntDistA->CopyVals();		
	fLookupIntCorrIrregularA->CopyVals();
	// AliInfo(" A side done");
	LOG(INFO) << " A side done" << FairLogger::endl;
      }
      if ( side == 1 ) {
	FillLookUpTableC(lookupGlobalDist,
			 fLookUpIntDistDrEzC,fLookUpIntDistDphiREzC,fLookUpIntDistDzC,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntDistC->CopyVals();				

	fLookupIntCorrIrregularC->CopyVals();
				
	// AliInfo(" C side done");		
	LOG(INFO) << " C side done" << FairLogger::endl;
      }
			
		
			
      // AliInfo("Step 6: Fill look up table for correction");	
      LOG(INFO) << "Step 6: Fill look up table for correction" << FairLogger::endl;

      FillLookUpTable(lookupGlobalCorr, 
		      fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz,
		      rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  
			
		
      // if ( side == 0 ) AliInfo(" A side done");
      // if ( side == 1 ) AliInfo(" C side done");
      if ( side == 0 ) LOG(INFO) << " A side done" << FairLogger::endl;
      if ( side == 1 ) LOG(INFO) << " C side done" << FairLogger::endl;

      delete lookupLocalDist;
      delete lookupLocalCorr;
    }
	
    //// copy to 1D interpolator /////
    fLookupIntDist->CopyVals();
    fLookupIntCorr->CopyVals();
		
    fInitLookUp = kTRUE;
  }

  //matricesErRet = matricesEr;
 
  // memory deallocation for temporary matrices
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    delete matricesGDistDrDz[k] ;
    delete matricesGDistDphiRDz[k]	;
    delete matricesGDistDz[k] ;
		
    delete matricesGCorrDrDz[k] ;
    delete matricesGCorrDphiRDz[k]	;
    delete matricesGCorrDz[k] ;
		
  }
  delete lookupGlobalDist;
  delete lookupGlobalCorr;
	
}


/// Creating look-up tables of Correction/Distortion by integration following 
/// drift line with known distributions for potential and spacecharge.
/// 
/// \param  matricesVA TMatrixD** potential distribution in side A (output)
/// \param  matricesChargeA TMatrixD** charge distribution in side A (input)
/// \param  matricesVC  TMatrixD** potential distribution in side C (output)
/// \param  matricesChargeC TMatrixD** charge distribution in side C (input)
/// \param rrow	Int_t  number of grid in row direction
///	\param zcolumn Int_t number of grid in z direction
/// \param phiSlice 	Int_t number of slicees in phi direction
/// \param maxIteration Int_t max iteration for convergence
/// \param stoppingConv Double_t stopping criteria for convergence
/// \param matricesErA TMatrixD** electric field r direction distribution side A
/// \param matricesEphiA TMatrixD** electric field phi direction distribution side A
/// \param matricesEzA TMatrixD** electric field z direction distribution side A
/// \param matricesErC TMatrixD** electric field r distribution side C
/// \param matricesEphiC TMatrixD** electric field phi distribution side C
/// \param matricesEzC TMatrixD** electric field z distribution side C
/// \param matricesDistDrDzA TMatrixD**  local r distortion (output) A side
/// \param matricesDistDphiRDzA TMatrixD** local r phi distortion (output) A side
/// \param matricesDistDzA TMatrixD**  local z distortion (output) A side
/// \param matricesCorrDrDzA TMatrixD** local r correction (output) A side
/// \param matricesCorrDphiRDzA TMatrixD** local r phi correction (output) A side
/// \param matricesCorrDzA 	TMatrixD** local z correction (output) A side
/// \param matricesDistDrDzC 	TMatrixD**   local r distortion (output) C side
/// \param matricesDistDphiRDzC 	TMatrixD**  local r phi distortion (output) C side
/// \param matricesDistDzC TMatrixD** local z distortion (output) C side
/// \param matricesCorrDrDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDphiRDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDzC	TMatrixD** local z correction (output) C side
///
/// \post Lookup tables for distortion: 
/// ~~~ 
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz  
/// ~~~ fo
/// and correction: 
/// ~~~ 
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz  
/// ~~~ 
/// are initialized 
/// 
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoissonIntegralDz
(
 TMatrixD** matricesVA,  
 TMatrixD** matricesChargeA,  
 TMatrixD** matricesVC,  
 TMatrixD** matricesChargeC,  				
 Int_t rrow, 
 Int_t zcolumn, 
 Int_t phiSlice,
 Int_t maxIteration,
 Double_t stoppingConv,
 TMatrixD** matricesErA,
 TMatrixD** matricesEphiA,
 TMatrixD** matricesEzA,
 TMatrixD** matricesErC,
 TMatrixD** matricesEphiC,
 TMatrixD** matricesEzC,
 TMatrixD** matricesDistDrDzA,
 TMatrixD** matricesDistDphiRDzA,
 TMatrixD** matricesDistDzA,
 TMatrixD** matricesCorrDrDzA,
 TMatrixD** matricesCorrDphiRDzA,
 TMatrixD** matricesCorrDzA,
 TMatrixD** matricesDistDrDzC,
 TMatrixD** matricesDistDphiRDzC,
 TMatrixD** matricesDistDzC,
 TMatrixD** matricesCorrDrDzC,
 TMatrixD** matricesCorrDphiRDzC,
 TMatrixD** matricesCorrDzC,
 TMatrixD** matricesLocalIntErDzA,
 TMatrixD** matricesLocalIntEphiDzA,
 TMatrixD** matricesLocalIntEzA,
 TMatrixD** matricesLocalIntErDzC,
 TMatrixD** matricesLocalIntEphiDzC,
 TMatrixD** matricesLocalIntEzC,
 TMatrixD** matricesGDistDrDzA,
 TMatrixD** matricesGDistDphiRDzA,
 TMatrixD** matricesGDistDzA,
 TMatrixD** matricesGDistDrDzC,
 TMatrixD** matricesGDistDphiRDzC,
 TMatrixD** matricesGDistDzC
 )
{
  // Compute grid size for all direction 
  Int_t phiSlicesPerSector = phiSlice / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	
  printf("gridSizeZ in init: %f\n",gridSizeZ);
	
  // local variables
  Float_t radius0, phi0, z0;
	
	
  // memory allocation for temporary matrices:
  // potential (boundary values), charge distribution
	
  TMatrixD** matricesV, **matricesCharge;
  TMatrixD** matricesEr, **matricesEphi, **matricesEz;
  TMatrixD** matricesDistDrDz;
  TMatrixD** matricesDistDphiRDz;
  TMatrixD** matricesDistDz;
  TMatrixD** matricesCorrDrDz;
  TMatrixD** matricesCorrDphiRDz;
  TMatrixD** matricesCorrDz;
  TMatrixD** matricesLocalIntErDz;
  TMatrixD** matricesLocalIntEphiDz;
  TMatrixD** matricesLocalIntEz;
  TMatrixD **matricesGDistDrDz;
  TMatrixD **matricesGDistDphiRDz;  
  TMatrixD **matricesGDistDz;
	
  //TMatrixD *matricesDistDrDz[phiSlice], *matricesDistDphiRDz[phiSlice],  *matricesDistDz[phiSlice];
  //TMatrixD *matricesCorrDrDz[phiSlice], *matricesCorrDphiRDz[phiSlice],  *matricesCorrDz[phiSlice];
  //TMatrixD *matricesGDistDrDz[phiSlice], *matricesGDistDphiRDz[phiSlice],  *matricesGDistDz[phiSlice];
  TMatrixD *matricesGCorrDrDz[phiSlice], *matricesGCorrDphiRDz[phiSlice],  *matricesGCorrDz[phiSlice];
	
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    matricesErA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesEphiA[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesEzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesErC[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesEphiC[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesEzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesDistDrDzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDphiRDzA[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesCorrDrDzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesCorrDphiRDzA[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesCorrDzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesDistDrDzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDphiRDzC[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesDistDzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesCorrDrDzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesCorrDphiRDzC[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesCorrDzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
		
    //matricesGDistDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    //matricesGDistDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    //matricesGDistDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesGCorrDrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesGCorrDphiRDz[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesGCorrDz[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
		
    //matricesLocalIntErDzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    //matricesLocalIntEphiDzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    //matricesLocalIntEzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    //matricesLocalIntErDzC[k]    	=   new TMatrixD(rrow,zcolumn) ;
    //matricesLocalIntEphiDzC[k]    	=   new TMatrixD(rrow,zcolumn) ;
    //dmatricesLocalIntEzC[k]    	=   new TMatrixD(rrow,zcolumn) ;
		
		
  }
	
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t  rlist[rrow], zedlist[zcolumn] , philist[phiSlice];
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) philist[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow ; i++ )  rlist[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn ; j++ ) zedlist[j]  = j * gridSizeZ ;

  // for irregular
  TMatrixD** matricesIrregularDrDz = NULL;
  TMatrixD** matricesIrregularDphiRDz = NULL;
  TMatrixD** matricesIrregularDz = NULL;
  TMatrixD** matricesPhiIrregular = NULL;
  TMatrixD** matricesRIrregular = NULL;
  TMatrixD** matricesZIrregular = NULL;
	
	
	
	
  // should be set, in another place
  const Int_t   symmetry = 0;

  // for charge
  TMatrixD** matricesLookUpCharge = NULL;
  AliTPC3DCylindricalInterpolator *chargeInterpolator;

	
  // do if look up table haven't be initialized
  if ( !fInitLookUp ) {
    //* do for 2 sides (z direction)
    // solve Poisson3D twice; once for +Z and once for -Z
		
		
    for ( Int_t side = 0 ; side < 2 ; side++ ) {  				
      if (side == 0) {
	matricesV = matricesVA;
	matricesCharge = matricesChargeA;
	matricesEr = matricesErA;
	matricesEphi = matricesEphiA;
	matricesEz = matricesEzA;
	matricesDistDrDz = matricesDistDrDzA;
	matricesDistDphiRDz = matricesDistDphiRDzA;
	matricesDistDz =matricesDistDzA;
	matricesCorrDrDz = matricesCorrDrDzA;
	matricesCorrDphiRDz =matricesCorrDphiRDzA;
	matricesCorrDz =matricesCorrDzA;			
	matricesLocalIntErDz     	=   matricesLocalIntErDzA;
	matricesLocalIntEphiDz     	=   matricesLocalIntEphiDzA ;
	matricesLocalIntEz    	=   matricesLocalIntEzA;
	matricesGDistDrDz = matricesGDistDrDzA;
	matricesGDistDphiRDz = matricesGDistDphiRDzA;
	matricesGDistDz = matricesGDistDzA;
	//lookupGlobalDist = fLookupIntDistA;
	//lookupGlobalCorr = fLookupIntCorrA;;


	matricesIrregularDrDz = fLookUpIntCorrDrEzIrregularA;
	matricesIrregularDphiRDz = fLookUpIntCorrDphiREzIrregularA;
	matricesIrregularDz = fLookUpIntCorrDzIrregularA;

	matricesPhiIrregular = fPhiListIrregularA;
	matricesRIrregular = fRListIrregularA;
	matricesZIrregular = fZListIrregularA;

	matricesLookUpCharge = fLookUpChargeA;
	chargeInterpolator = fInterpolatorChargeA;
			



	fLookupDistA->SetLookUpR(matricesDistDrDz);
	fLookupDistA->SetLookUpPhi(matricesDistDphiRDz);
	fLookupDistA->SetLookUpZ(matricesDistDz);

		
				
      } else {
	matricesV = matricesVC;
	matricesCharge = matricesChargeC;
	matricesEr = matricesErC;
	matricesEphi = matricesEphiC;
	matricesEz = matricesEzC;
	matricesDistDrDz = matricesDistDrDzC;
	matricesDistDphiRDz = matricesDistDphiRDzC;
	matricesDistDz =matricesDistDzC;
	matricesCorrDrDz = matricesCorrDrDzC;
	matricesCorrDphiRDz =matricesCorrDphiRDzC;
	matricesCorrDz =matricesCorrDzC;			
				
	matricesLocalIntErDz     	=   matricesLocalIntErDzC;
	matricesLocalIntEphiDz     	=   matricesLocalIntEphiDzC;
	matricesLocalIntEz    	=   matricesLocalIntEzC;
	matricesGDistDrDz = matricesGDistDrDzC;
	matricesGDistDphiRDz = matricesGDistDphiRDzC;
	matricesGDistDz = matricesGDistDzC;
		
	//lookupGlobalDist = fLookupIntDistC;
	//lookupGlobalCorr = fLookupIntCorrC;
	matricesIrregularDrDz = fLookUpIntCorrDrEzIrregularC;
	matricesIrregularDphiRDz = fLookUpIntCorrDphiREzIrregularC;
	matricesIrregularDz = fLookUpIntCorrDzIrregularC;

	matricesPhiIrregular = fPhiListIrregularC;
	matricesRIrregular = fRListIrregularC;
	matricesZIrregular = fZListIrregularC;

	matricesLookUpCharge = fLookUpChargeC;
	chargeInterpolator = fInterpolatorChargeC;


	fLookupDistC->SetLookUpR(matricesDistDrDz);
	fLookupDistC->SetLookUpPhi(matricesDistDphiRDz);
	fLookupDistC->SetLookUpZ(matricesDistDz);
			
		
      } 
			
			

      // setvals for lookup charge

      for ( Int_t k = 0 ; k < phiSlice ; k++ ) 	*(matricesLookUpCharge[k]) = *(matricesCharge[k]);

      chargeInterpolator->SetVals(matricesLookUpCharge);
      chargeInterpolator->InitCubicSpline();

      // AliInfo(Form("Step 0: Preparing Charge interpolator: \n"));
      LOG(INFO) << Form("Step 0: Preparing Charge interpolator: \n") << FairLogger::endl;

			
      AliTPCLookUpTable3DInterpolatorD *lookupLocalDist = 
	new AliTPCLookUpTable3DInterpolatorD(
					     rrow,matricesDistDrDz,rlist,phiSlice,matricesDistDphiRDz,philist,zcolumn,matricesDistDz,zedlist,fInterpolationOrder);
			
      AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr = 
	new AliTPCLookUpTable3DInterpolatorD(
					     rrow,matricesCorrDrDz,rlist,phiSlice,matricesCorrDphiRDz,philist,zcolumn,matricesCorrDz,zedlist,fInterpolationOrder);
	
	
      AliTPCLookUpTable3DInterpolatorD *lookupGlobalDist = 
	new AliTPCLookUpTable3DInterpolatorD(
					     rrow,matricesGDistDrDz,rlist,phiSlice,matricesGDistDphiRDz,philist,zcolumn,matricesGDistDz,zedlist,fInterpolationOrder);
			
      AliTPCLookUpTable3DInterpolatorD *lookupGlobalCorr = 
	new AliTPCLookUpTable3DInterpolatorD(
					     rrow,matricesGCorrDrDz,rlist,phiSlice,matricesGCorrDphiRDz,philist,zcolumn,matricesGCorrDz,zedlist,fInterpolationOrder);
		
      for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
	matricesGDistDrDz[k]->Zero();
	matricesGDistDphiRDz[k]->Zero();
	matricesGDistDz[k]->Zero();
		
	matricesGCorrDrDz[k]->Zero();
	matricesGCorrDphiRDz[k]->Zero();
	matricesGCorrDz[k]->Zero();
		
      }
	
      // AliInfo("Step 1: Solving poisson solver");
      LOG(INFO) << "Step 1: Solving poisson solver" << FairLogger::endl;

			
			
			
	
      fPoissonSolver->PoissonSolver3D( matricesV, matricesCharge, rrow, zcolumn, phiSlice, maxIteration, symmetry ) ;
	
		
      // AliInfo("Step 2: Calculate electric field");
      LOG(INFO) << "Step 2: Calculate electric field" << FairLogger::endl;
      ElectricField( matricesV, 
		     matricesEr,  matricesEphi, matricesEz, rrow, zcolumn, phiSlice, 
		     gridSizeR, gridSizePhi ,gridSizeZ,symmetry, o2::TPC::AliTPCPoissonSolver::fgkIFCRadius); 
			
			
      // AliInfo("Step 3: Calculate local distortion");
      LOG(INFO) << "Step 3: Calculate local distortion" << FairLogger::endl;

			
      LocalDistCorrDz (matricesEr, matricesEphi, 	matricesEz,  
		       matricesDistDrDz,  matricesDistDphiRDz, matricesDistDz,
		       matricesCorrDrDz,  matricesCorrDphiRDz, matricesCorrDz,
		       matricesLocalIntErDz, matricesLocalIntEphiDz, matricesLocalIntEz,		
		       rrow,  zcolumn, phiSlice, gridSizeZ, ezField);	
			
			
      //// copy to 1D interpolator /////
      lookupLocalDist->CopyVals();
      lookupLocalCorr->CopyVals();


      if (side == 0) 
	fLookupDistA->CopyVals();
      else
	fLookupDistC->CopyVals();

      //// 
		
      // AliInfo("Step 4a: Integrate distortion and correction the drift line --- test Correctness == ");
      LOG(INFO) << "Step 4a: Integrate distortion and correction the drift line --- test Correctness == " << FairLogger::endl;

      if (fStrategy == kNaive) 
	IntegrateDistCorrDriftLineDz(
				     lookupLocalDist, 
				     matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
				     lookupLocalCorr, 
				     matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  

				     matricesIrregularDrDz, matricesIrregularDphiRDz,  matricesIrregularDz,
				     matricesRIrregular, matricesPhiIrregular,matricesZIrregular,	


				     rrow,  zcolumn, phiSlice, rlist, philist, zedlist
				     );  
      else
	IntegrateDistCorrDriftLineDzOpt2(matricesDistDrDz,  matricesDistDphiRDz, matricesDistDz, 
					 matricesGDistDrDz,  matricesGDistDphiRDz, matricesGDistDz,  
					 matricesCorrDrDz,  matricesCorrDphiRDz, matricesCorrDz, 
					 matricesGCorrDrDz,  matricesGCorrDphiRDz, matricesGCorrDz,  
					 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  
			
			
			
      // AliInfo("Step 5: Fill look up table for distortion");	
      LOG(INFO) << "Step 5: Fill look up table for distortion" << FairLogger::endl;

      //// copy to 1D interpolator /////
      lookupGlobalDist->CopyVals();
      lookupGlobalCorr->CopyVals();
      //// 
		
		
		
		
      //FillLookUpTable(lookupGlobalDist,
      //	fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz,
      //	rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  
		
			
      // AliInfo("Step 6: Fill look up table for correction");	
      LOG(INFO) << "Step 6: Fill look up table for correction" << FairLogger::endl;

      //FillLookUpTable(lookupGlobalCorr,
      //	fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz,
      //	rrow,  zcolumn, phiSlice, rlist, philist, zedlist, side);  
			

		
      if ( side == 0 ) {
	FillLookUpTableA(lookupGlobalDist,
			 fLookUpIntDistDrEzA,fLookUpIntDistDphiREzA,fLookUpIntDistDzA,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntDistA->CopyVals();		

	fLookupIntCorrIrregularA->CopyVals();
	// AliInfo(" A side done");
	LOG(INFO) << " A side done" << FairLogger::endl;
      }
      if ( side == 1 ) {
	FillLookUpTableC(lookupGlobalDist,
			 fLookUpIntDistDrEzC,fLookUpIntDistDphiREzC,fLookUpIntDistDzC,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntDistC->CopyVals();				

	fLookupIntCorrIrregularC->CopyVals();
				
	// AliInfo(" C side done");		
	LOG(INFO) << " C side done" << FairLogger::endl;
      }
			
			
      delete lookupLocalDist;
      delete lookupLocalCorr;
      delete lookupGlobalDist;
      delete lookupGlobalCorr;

			
    }
		
    fInitLookUp = kTRUE;
    //// copy to 1D interpolator /////
    //fLookupIntDist->CopyVals();		
		
		
    fLookupIntCorr->CopyVals();
		
		
  }
 
  // memory deallocation for temporary matrices
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    //delete matricesGDistDrDz[k] ;
    //delete matricesGDistDphiRDz[k]	;
    //delete matricesGDistDz[k] ;
		
    delete matricesGCorrDrDz[k] ;
    delete matricesGCorrDphiRDz[k]	;
    delete matricesGCorrDz[k] ;
		
  }
	
}


/// Creating look-up tables of Correction/Distortion by integration following 
/// drift line with known distributions for potential and spacecharge.
/// 
/// \param  matricesVA TMatrixD** potential distribution in side A (output)
/// \param  matricesChargeA TMatrixD** charge distribution in side A (input)
/// \param  matricesVC  TMatrixD** potential distribution in side C (output)
/// \param  matricesChargeC TMatrixD** charge distribution in side C (input)
/// \param rrow	Int_t  number of grid in row direction
///	\param zcolumn Int_t number of grid in z direction
/// \param phiSlice 	Int_t number of slicees in phi direction
/// \param maxIteration Int_t max iteration for convergence
/// \param stoppingConv Double_t stopping criteria for convergence
/// \param matricesErA TMatrixD** electric field r direction distribution side A
/// \param matricesEphiA TMatrixD** electric field phi direction distribution side A
/// \param matricesEzA TMatrixD** electric field z direction distribution side A
/// \param matricesErC TMatrixD** electric field r distribution side C
/// \param matricesEphiC TMatrixD** electric field phi distribution side C
/// \param matricesEzC TMatrixD** electric field z distribution side C
/// \param matricesDistDrDzA TMatrixD**  local r distortion (output) A side
/// \param matricesDistDphiRDzA TMatrixD** local r phi distortion (output) A side
/// \param matricesDistDzA TMatrixD**  local z distortion (output) A side
/// \param matricesCorrDrDzA TMatrixD** local r correction (output) A side
/// \param matricesCorrDphiRDzA TMatrixD** local r phi correction (output) A side
/// \param matricesCorrDzA 	TMatrixD** local z correction (output) A side
/// \param matricesDistDrDzC 	TMatrixD**   local r distortion (output) C side
/// \param matricesDistDphiRDzC 	TMatrixD**  local r phi distortion (output) C side
/// \param matricesDistDzC TMatrixD** local z distortion (output) C side
/// \param matricesCorrDrDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDphiRDzC TMatrixD** local r phi correction (output) C side
/// \param matricesCorrDzC	TMatrixD** local z correction (output) C side
///
/// \post Lookup tables for distortion: 
/// ~~~ 
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz  
/// ~~~ fo
/// and correction: 
/// ~~~ 
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz  
/// ~~~ 
/// are initialized 
/// 
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoisson
(
 TMatrixD** matricesVA,  
 TMatrixD** matricesChargeA,  
 TMatrixD** matricesVC,  
 TMatrixD** matricesChargeC,  				
 Int_t rrow, 
 Int_t zcolumn, 
 Int_t phiSlice,
 Int_t maxIteration,
 Double_t stoppingConv,
 TMatrixD** matricesErA,
 TMatrixD** matricesEphiA,
 TMatrixD** matricesEzA,
 TMatrixD** matricesErC,
 TMatrixD** matricesEphiC,
 TMatrixD** matricesEzC
 )
{
  // Compute grid size for all direction 
  Int_t phiSlicesPerSector = phiSlice / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	
  printf("gridSizeZ in init: %f\n",gridSizeZ);
	
  // local variables
  Float_t radius0, phi0, z0;
	
	
  // memory allocation for temporary matrices:
  // potential (boundary values), charge distribution
	
  TMatrixD** matricesV, **matricesCharge;
  TMatrixD** matricesEr, **matricesEphi, **matricesEz;
	
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    matricesErA[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesEphiA[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesEzA[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
    matricesErC[k]     	=   new TMatrixD(rrow,zcolumn) ;
    matricesEphiC[k]	=   new TMatrixD(rrow,zcolumn) ;
    matricesEzC[k]     	=   new TMatrixD(rrow,zcolumn) ;
		
		
  }
	
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t  rlist[rrow], zedlist[zcolumn] , philist[phiSlice];
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) philist[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow ; i++ )  rlist[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn ; j++ ) zedlist[j]  = j * gridSizeZ ;
	
	
	
	
  // should be set, in another place
  const Int_t   symmetry = 0;
	
  // do if look up table haven't be initialized
  if ( !fInitLookUp ) {
    //* do for 2 sides (z direction)
			
		
    for ( Int_t side = 0 ; side < 2 ; side++ ) {  				
      if (side == 0) {
	matricesV = matricesVA;
	matricesCharge = matricesChargeA;
	matricesEr = matricesErA;
	matricesEphi = matricesEphiA;
	matricesEz = matricesEzA;
		
				
      } else {
	matricesV = matricesVC;
	matricesCharge = matricesChargeC;
	matricesEr = matricesErC;
	matricesEphi = matricesEphiC;
	matricesEz = matricesEzC;
		
      } 
			
			
      AliTPCLookUpTable3DInterpolatorD *lookupEField = 
	new AliTPCLookUpTable3DInterpolatorD(
					     rrow,
					     matricesEr,
					     rlist,phiSlice,
					     matricesEphi,
					     philist,zcolumn,
					     matricesEz,
					     zedlist,
					     fInterpolationOrder
					     );
			
			
	
      // AliInfo("Step 1: Solving poisson solver");
      LOG(INFO) << "Step 1: Solving poisson solver" << FairLogger::endl;

			
			
			
	
      fPoissonSolver->PoissonSolver3D( matricesV, matricesCharge, rrow, zcolumn, phiSlice, maxIteration, symmetry ) ;
	
	
		
      // AliInfo("Step 2: Calculate electric field");
      LOG(INFO) << "Step 2: Calculate electric field" << FairLogger::endl;
      CalculateEField(
		      matricesV, 
		      matricesEr,  
		      matricesEphi, 
		      matricesEz, 
		      rrow, 
		      zcolumn, 
		      phiSlice, 
		      maxIteration,   
		      symmetry 
		      );
		
      lookupEField->CopyVals();		
      // AliInfo("Step 3: Fill the ");
      LOG(INFO) << "Step 3: Fill the " << FairLogger::endl;

			
			
      if ( side == 0 ) {
	FillLookUpTableA(lookupEField,
			 fLookUpErOverEzA,fLookUpEphiOverEzA,fLookUpDeltaEzA,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
				
	fLookupIntENoDriftA->CopyVals();		
	// AliInfo(" A side done");
	LOG(INFO) << " A side done" << FairLogger::endl;
      }
      if ( side == 1 ) {
	FillLookUpTableC(lookupEField,
			 fLookUpErOverEzC,fLookUpEphiOverEzC,fLookUpDeltaEzC,
			 rrow,  zcolumn, phiSlice, rlist, philist, zedlist);  			
	fLookupIntENoDriftC->CopyVals();		
				
	// AliInfo(" C side done");		
	LOG(INFO) << " C side done" << FairLogger::endl;
      }
			
			
      delete lookupEField;

			
    }
		
    fInitLookUp = kTRUE;
		
		
		
  }
 
	
	
}
/// Force creating look-up table of Correction/Distortion by integration following 
/// drift line.
///
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
/// \param maxIteration Int_t Maximum iteration for poisson solver
/// \param stoppingConv Convergence error stopping conditioin for poisson solver
/// 
void AliTPCSpaceCharge3DDriftLine::ForceInitSpaceCharge3DPoissonIntegralDz
(
 Int_t rrow, 
 Int_t zcolumn, 
 Int_t phiSlice,
 Int_t maxIteration,
 Double_t stoppingConv
 )
{
  fInitLookUp = kFALSE;
  InitSpaceCharge3DPoissonIntegralDz(
				     rrow, 
				     zcolumn, 
				     phiSlice,
				     maxIteration,
				     stoppingConv
				     );
	
}

///
/// Electricfield Calculation: 
///
/// \param arrayofArrayV TMatrixD** 3D matrix representing calculated potential 
/// \param arrayofArrayEr TMatrix** 3D matrix representing e-field at Er
/// \param arrayofArrayEz TMatrix** 3D matrix representing e-field at Ez
/// \param arrayofArrayEphi TMatrix** 3D matrix representing e-field at Ephi
/// \param rows Int_t number of rows of discritization (in R direction)
/// \param columns Int_t number of columns  of discritization (in Z direction)
/// \param phislices Int_t number of (phislices in phi direction) 
/// \param symmetry Int_t symmetry?
/// 
/// \pre   Matrix arrayofArrayV is assumed had been calculated  by Poisson solver
/// \post  Results of  E-fields are calculated by measuring gradient at potential distribution
///
///
///	* Differentiate potential on all direction (r,z and phi)
/// * Non-boundary -> Central difference (3 stencil) TODO: 5 Stencil
///
///   \f$  \nabla_{r} V(r_{i},\phi_{j},z_{k}) \approx -( V_{i+1,j,k} - V_{i-1,j,k}) / (2* h_{r}) \f$
///
///   \f$ -\nabla_{\phi} V(r_{i},\phi_{j},z_{k}) \approx -( V_{i,j-1,k} - V_{i,j+1,k}) / (2* r_{j} * h_{\phi}) \f$
///
///   \f$ -\nabla_{z} V(r_{i},\phi_{j},z_{k}) \approx -( V_{i,j,k+1} - V_{i,j,k-1}) / (2* h_{z}) \f$
///
///   ~~~ cxx
///   arrayEr(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR); // r direction				
///		arrayEz(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ; // z direction				
///		arrayEphi(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ; // phi didrection			      		
///   ~~~
///
/// * Boundary -> Forward/Backward difference (3 stencil) TODO: 5 Stencil
///
///   \f$ -\nabla_{r} V(r_{0},\phi_{j},z_{k}) \approx -( -0.5 V_{2,j,k} + 2 V_{1,j,k} - 1.5 * V_{0,j,k}) /  h_{r} \f$
///
///   \f$ -\nabla_{r} V(r_{rrow - 1},\phi_{j},z_{k}) \approx -( 1.5 V_{rrow-1,j,k} - 2.0 V_{rrow-2,j,k} + 0.5 V_{rrow -3,j,k}) / h_{\phi} \f$
///
void AliTPCSpaceCharge3DDriftLine::ElectricField
(
 TMatrixD** matricesV,  
 TMatrixD** matricesEr,  
 TMatrixD** matricesEphi, 
 TMatrixD** matricesEz,
 const Int_t rrow, 
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Float_t gridSizeR, 
 const Float_t  gridSizePhi ,
 const Float_t  gridSizeZ,
 const Int_t symmetry, 
 const Float_t  innerRadius 
 )
{
					
	
  Float_t radius;
  Int_t mplus, mminus, signplus, signminus  ;
  // iterate over phislices
	
	
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {		
    mplus  = m + 1;   signplus  = 1 ;
    mminus = m - 1 ;  signminus = 1 ;
    if (symmetry==1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if ( mplus  > phiSlice-1 ) mplus  = phiSlice - 2 ;
      if ( mminus < 0 )           mminus = 1 ;
    }
    else if (symmetry==-1) {       // Anti-symmetry in phi
      if ( mplus  > phiSlice-1 ) { mplus  = phiSlice - 2 ;  signplus  = -1 ; }
      if ( mminus < 0 )           { mminus = 1 ;	            signminus = -1 ; }
    }
    else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if ( mplus  > phiSlice-1 ) mplus  = m + 1 - phiSlice ;
      if ( mminus < 0 )           mminus = m - 1 + phiSlice ;
    }
		
    TMatrixD &arrayVP     =  *matricesV[mplus] ;
    TMatrixD &arrayVM     =  *matricesV[mminus] ;    
    TMatrixD& arrayV    =  *matricesV[m] ;
    TMatrixD& arrayEr	  =  *matricesEr[m] ;
    TMatrixD& arrayEz	  =  *matricesEz[m] ;
    TMatrixD& arrayEphi	 =  *matricesEphi[m] ;
	
    // for non-boundary V
    for ( Int_t i = 1 ; i < rrow-1; i++) {
      radius = innerRadius  + i*gridSizeR ;
      for ( Int_t j = 1; j < zcolumn-1 ; j++ ) {				
	arrayEr(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR); // r direction				
	arrayEz(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ; // z direction				
	arrayEphi(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ; // phi didrection			      		
      }
    }
		
    // for boundary-r
    for ( Int_t j = 0; j < zcolumn ; j++ )  {
      arrayEr(0,j)      =  -1 * ( -0.5*arrayV(2,j) + 2.0*arrayV(1,j) - 1.5*arrayV(0,j) ) / gridSizeR ; // forward difference
      arrayEr(rrow-1,j) =  -1 * ( 1.5*arrayV(rrow-1,j) - 2.0*arrayV(rrow-2,j) + 0.5*arrayV(rrow-3,j) ) / gridSizeR ; // backward difference
    }

    for ( Int_t i = 0 ; i < rrow; i += rrow-1) {
      radius = innerRadius + i*gridSizeR ;
      for ( Int_t j = 1; j < zcolumn-1 ; j++ ) {
	arrayEz(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ; // z direction				
	arrayEphi(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ; // phi didrection			      		
      }
    }
		
    // for boundary-z
    for ( Int_t i = 0; i < rrow ; i++ )  {
      arrayEz(i,0)         =  -1 * ( -0.5*arrayV(i,2) + 2.0*arrayV(i,1) - 1.5*arrayV(i,0) ) / gridSizeZ ;
      arrayEz(i,zcolumn-1) =  -1 * ( 1.5*arrayV(i,zcolumn-1) - 2.0*arrayV(i,zcolumn-2) + 0.5*arrayV(i,zcolumn-3) ) / gridSizeZ ;			
    }

    for ( Int_t i = 1 ; i < rrow-1; i ++) {
      radius = innerRadius + i*gridSizeR ;
      for ( Int_t j = 0; j < zcolumn ; j += zcolumn-1 ) {
	arrayEr(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR); // r direction								
	arrayEphi(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ; // phi didrection			      		
      }
    }
		
    // corner points for Ephi
    for ( Int_t i = 0 ; i < rrow; i += rrow - 1) {
      radius = innerRadius  + i*gridSizeR ;
      for ( Int_t j = 0; j < zcolumn ; j += zcolumn-1 ) {
	arrayEphi(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ; // phi didrection			      		
      }
    }		
  }		
}




///
/// Local distortion and correction, calculate local distortion/correction 
/// based on simplified langevin equation, see internal note ALICE-INT-2010-016.
///
/// <b> Local Distortion </b> 
///
/// Local distortion is calculated based on formulation in ALICE-INT-2010-016, this function assume that
/// electric field \f$\vec{E}(r_{i},z_{j},\phi_{m})\f$ is provided.
///
/// First, we calculate integration of the Electric field in z-direction for all direction. 
/// Assumption: \f$ z_{0} \f$ is location of CE (Central Electrode) and \f$ z_{zcolumn - 1} \f$ is location of End Plate.
///
/// This integration is in \f$z\f$ direction. If no change in granurality, then we can only use trapezoidal rule.
///
/// Let suppose we want to calculate local distortion at \f$(r_{i},z_{j},\phi_{m})\f$. 
/// Assume \f$\vec{E}(r_{i},z_{j+1},\phi_{m}) \f$ and \f$\vec{E}(r_{i},z_{j},\phi_{m}) \f$ are known, see Figure \ref fig1 (a),
///
/// \anchor fig1
/// ![Local Distortion](localdist.png)
///
/// Than we can calculate definite integrations for each directions in respect of $z$  from \f$ z_{j} \f$  to \f$ z_{j + 1} \f$   as follows:
/// 
/// \f$  \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{r}(r_{i},z_{j},\phi_{m}) + E_{r}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx  \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{\phi}(r_{i},z_{j},\phi_{m}) + E_{\phi}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j+1}}_{z_{j}} E_{z}(r_{i},z_{j},\phi_{m}) dz  \approx   \frac{h_{z}}{2.0} \left( E_{z}(r_{i},z_{j},\phi_{m}) + E_{z}(r_{i},z_{j+1},\phi_{m}) \right) \f$
///
/// The code sniplet at \ref impllocaldist is an implementation of the local integration of electric field.
///
/// \anchor impllocaldist
/// ~~~
/// Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;
///
/// localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1))/(-1*ezField) ;
/// localIntEphiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1))/(-1*ezField) ;
/// localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;				
/// ~~~
///
///
/// After we have local integrations for ellectric fields in each direction, 
/// local distortion \f$\hat{\delta}(r_{i},z_{j},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref1 (b) for illustration):
/// 
/// \f$ \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) = c_{0} \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
/// ~~~ 
///
/// \f$ r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m})  = - c_{1} \int^{z_{j+1}}_{z_{j}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j+1}}_{j_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
/// ~~~ 
///
/// \f$ \hat{\delta}_{z}(r_{i},z_{j},\phi_{m})  = \int_{z_{j}}^{z_{j+1}} \frac{v^{\prime}(E)}{v_{0}} (E - E_{0}) dz\f$ 
///
/// ~~~
/// (*distDz)(i,j) = localIntDeltaEz*o2::TPC::AliTPCPoissonSolver::fgkdvdE; 
/// ~~~
///
/// Where \f$c_{0}\f$ and \f$c_{1}\f$ are constants (see the ALICE-INT-2010-016 for further details).
///
/// <b> Local correction </b>
///
/// Local correction is computed as local distortion where the electric fields are in opposite direction (see Figure \ref fig2 (a)).
/// 
/// \anchor fig2
/// ![Local Correction](localcorr.png)
///
/// Let suppose we want to calculate local correction at \f$(r_{i},\mathbf{z_{j+1}},\phi_{m})\f$. 
/// Assume \f$\vec{E}(r_{i},z_{j+1},\phi_{m}) \f$ and \f$\vec{E}(r_{i},z_{j},\phi_{m}) \f$ are known.
///
/// Than we can calculate definite integrations for each directions in respect of \f$z\f$  from \f$ z_{j+1} \f$  to \f$ z_{j} \f$   as follows:
/// 
/// \f$  \int^{z_{j}}_{z_{j+1}} \frac{E_{r}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx -1 * \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{r}(r_{i},z_{j},\phi_{m}) + E_{r}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j}}_{z_{j+1}} \frac{E_{\phi}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx  -1 *  \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{\phi}(r_{i},z_{j},\phi_{m}) + E_{\phi}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j}}_{z_{j+1}} E_{z}(r_{i},z_{j},\phi_{m}) dz  \approx  -1 *   \frac{h_{z}}{2.0} \left( E_{z}(r_{i},z_{j},\phi_{m}) + E_{z}(r_{i},z_{j+1},\phi_{m}) \right) \f$
///
/// Local correction at \f$\hat{\delta'}(r_{i},\mathbf{z_{j+1}},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref fig2 (b) for illustration):
///
/// \f$ \hat{\delta'}_{rE}(r_{i},z_{j+1},\phi_{m}) = c_{0} \int^{z_{j}}_{z_{j+1}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j-1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// \f$ r\hat{\delta'}_{\phi E}(r_{i},z_{j+1},\phi_{m})  = - c_{1} \int^{z_{j}}_{z_{j+1}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j-1}}_{j_{k}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// \f$ \hat{\delta'}_{z}(r_{i},z_{j+1},\phi_{m})  = \int_{z_{j}}^{z_{j+1}} \frac{v^{\prime}(E)}{v_{0}} (E - E_{0}) dz\f$ 
///
/// For implementation, we use the fact that  
///
/// \f$ \hat{\delta'}_{rE}(r_{i},z_{j+1},\phi_{m}) = -1 * \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) \f$
///
/// \f$ r\hat{\delta'}_{\phi E}(r_{i},z_{j+1},\phi_{m}) =  -1 *  r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m}) \f$
///
/// \f$ \hat{\delta'}_{z}(r_{i},z_{j+1},\phi_{m}) =  -1 *  \hat{\delta}_{z}(r_{i},z_{j},\phi_{m}) \f$
///
/// ~~~
///	(*corrDrDz)(i,j+1) 		= -1* (*distDrDz)(i,j) ;				
/// (*corrDphiRDz)(i,j+1) = -1* (*distDphiRDz)(i,j);
/// (*corrDz)(i,j+1)      = -1* (*distDz)(i,j);				
/// ~~~
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesDistDrDz TMatrixD**  local distortion \f$\hat{\delta}_{r}\f$
///	\param matricesDistDphiRDz TMatrixD** local distortion \f$r \hat{\delta}_{\phi}\f$
///	\param matricesDistDz TMatrixD**   local distortion \f$ \hat{\delta}_{z}\f$
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez assummed already been calculated
/// \post Local distortion and correction are computed according simplified Langevin equation
/// ~~~ 
/// matricesDistDrDz,matricesDistDphiRDz,matricesDistDz
/// ~~~ 
/// and correction: 
/// ~~~ 
/// matricesCorrDrDz,matricesCorrDphiRDz,matricesCorrDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::LocalDistCorrDz
( 
 TMatrixD** matricesEr, 	
 TMatrixD** matricesEphi, 	
 TMatrixD** matricesEz, 	
 TMatrixD** matricesDistDrDz, 
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 		
 TMatrixD** matricesCorrDrDz, 
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 const Int_t rrow, 
 const Int_t zcolumn, 
 const Int_t phiSlice,  
 const Float_t gridSizeZ,
 const Double_t ezField
  )
{
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEphiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;
	
  // pointer declaration
  TMatrixD * eR;
  TMatrixD * ePhi;
  TMatrixD * eZ;
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
	
	
  // Initialization for j == colomn-1 integration is 0.0
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
    for (Int_t i=0; i< rrow;i++) {
      (*distDrDz)(i,zcolumn-1) =  0.0 ;			
      (*distDphiRDz)(i,zcolumn-1) =  0.0 ;			
      (*distDz)(i,zcolumn-1) =  0.0 ;			
			
      (*corrDrDz)(i,0) =  0.0 ;			
      (*corrDphiRDz)(i,0) =  0.0 ;			
      (*corrDz)(i,0) =  0.0 ;			
    }			
  }
	
  // for this case
  // use trapezoidal rule assume no ROC displacement
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    eR  =  matricesEr[m] ;
    ePhi  =  matricesEphi[m] ;
    eZ  =  matricesEz[m] ;
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
    for (Int_t j=0; j< zcolumn-1;j++) {
      for (Int_t i=0; i< rrow;i++)  {
	localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1))/(-1*ezField) ;
	localIntEphiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1))/(-1*ezField) ;				
	localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;				
				
				
	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
	(*distDz)(i,j)      = localIntDeltaEz*o2::TPC::AliTPCPoissonSolver::fgkdvdE*o2::TPC::AliTPCPoissonSolver::fgkdvdE;// two times?				
				
				
	(*corrDrDz)(i,j+1) 		= -1* (*distDrDz)(i,j) ;				
	(*corrDphiRDz)(i,j+1) = -1* (*distDphiRDz)(i,j);
	(*corrDz)(i,j+1)      = -1* (*distDz)(i,j);							
				
      }
    }				
  }		
}


///
/// Local distortion and correction, calculate local distortion/correction 
/// based on simplified langevin equation, see internal note ALICE-INT-2010-016.
///
/// <b> Local Distortion </b> 
///
/// Local distortion is calculated based on formulation in ALICE-INT-2010-016, this function assume that
/// electric field \f$\vec{E}(r_{i},z_{j},\phi_{m})\f$ is provided.
///
/// First, we calculate integration of the Electric field in z-direction for all direction. 
/// Assumption: \f$ z_{0} \f$ is location of CE (Central Electrode) and \f$ z_{zcolumn - 1} \f$ is location of End Plate.
///
/// This integration is in \f$z\f$ direction. If no change in granurality, then we can only use trapezoidal rule.
///
/// Let suppose we want to calculate local distortion at \f$(r_{i},z_{j},\phi_{m})\f$. 
/// Assume \f$\vec{E}(r_{i},z_{j+1},\phi_{m}) \f$ and \f$\vec{E}(r_{i},z_{j},\phi_{m}) \f$ are known, see Figure \ref fig1 (a),
///
/// \anchor fig1
/// ![Local Distortion](localdist.png)
///
/// Than we can calculate definite integrations for each directions in respect of $z$  from \f$ z_{j} \f$  to \f$ z_{j + 1} \f$   as follows:
/// 
/// \f$  \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{r}(r_{i},z_{j},\phi_{m}) + E_{r}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx  \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{\phi}(r_{i},z_{j},\phi_{m}) + E_{\phi}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j+1}}_{z_{j}} E_{z}(r_{i},z_{j},\phi_{m}) dz  \approx   \frac{h_{z}}{2.0} \left( E_{z}(r_{i},z_{j},\phi_{m}) + E_{z}(r_{i},z_{j+1},\phi_{m}) \right) \f$
///
/// The code sniplet at \ref impllocaldist is an impplementation of the local integration of electric field.
///
/// \anchor impllocaldist
/// ~~~
/// Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;
///
/// localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1))/(-1*ezField) ;
/// localIntEphiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1))/(-1*ezField) ;
/// localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;				
/// ~~~
///
///
/// After we have local integrations for ellectric fields in each direction, 
/// local distortion \f$\hat{\delta}(r_{i},z_{j},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref1 (b) for illustration):
/// 
/// \f$ \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) = c_{0} \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
/// ~~~ 
///
/// \f$ r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m})  = - c_{1} \int^{z_{j+1}}_{z_{j}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j+1}}_{j_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
/// ~~~ 
///
/// \f$ \hat{\delta}_{z}(r_{i},z_{j},\phi_{m})  = \int_{z_{j}}^{z_{j+1}} \frac{v^{\prime}(E)}{v_{0}} (E - E_{0}) dz\f$ 
///
/// ~~~
/// (*distDz)(i,j) = localIntDeltaEz*o2::TPC::AliTPCPoissonSolver::fgkdvdE; 
/// ~~~
///
/// Where \f$c_{0}\f$ and \f$c_{1}\f$ are constants (see the ALICE-INT-2010-016 for further details).
///
/// <b> Local correction </b>
///
/// Local correction is computed as local distortion where the electric fields are in opposite direction (see Figure \ref fig2 (a)).
/// 
/// \anchor fig2
/// ![Local Correction](localcorr.png)
///
/// Let suppose we want to calculate local correction at \f$(r_{i},\mathbf{z_{j+1}},\phi_{m})\f$. 
/// Assume \f$\vec{E}(r_{i},z_{j+1},\phi_{m}) \f$ and \f$\vec{E}(r_{i},z_{j},\phi_{m}) \f$ are known.
///
/// Than we can calculate definite integrations for each directions in respect of \f$z\f$  from \f$ z_{j+1} \f$  to \f$ z_{j} \f$   as follows:
/// 
/// \f$  \int^{z_{j}}_{z_{j+1}} \frac{E_{r}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx -1 * \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{r}(r_{i},z_{j},\phi_{m}) + E_{r}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j}}_{z_{j+1}} \frac{E_{\phi}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx  -1 *  \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{\phi}(r_{i},z_{j},\phi_{m}) + E_{\phi}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j}}_{z_{j+1}} E_{z}(r_{i},z_{j},\phi_{m}) dz  \approx  -1 *   \frac{h_{z}}{2.0} \left( E_{z}(r_{i},z_{j},\phi_{m}) + E_{z}(r_{i},z_{j+1},\phi_{m}) \right) \f$
///
/// Local correction at \f$\hat{\delta'}(r_{i},\mathbf{z_{j+1}},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref fig2 (b) for illustration):
///
/// \f$ \hat{\delta'}_{rE}(r_{i},z_{j+1},\phi_{m}) = c_{0} \int^{z_{j}}_{z_{j+1}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j-1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// \f$ r\hat{\delta'}_{\phi E}(r_{i},z_{j+1},\phi_{m})  = - c_{1} \int^{z_{j}}_{z_{j+1}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j-1}}_{j_{k}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// \f$ \hat{\delta'}_{z}(r_{i},z_{j+1},\phi_{m})  = \int_{z_{j}}^{z_{j+1}} \frac{v^{\prime}(E)}{v_{0}} (E - E_{0}) dz\f$ 
///
/// For implementation, we use the fact that  
///
/// \f$ \hat{\delta'}_{rE}(r_{i},z_{j+1},\phi_{m}) = -1 * \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) \f$
///
/// \f$ r\hat{\delta'}_{\phi E}(r_{i},z_{j+1},\phi_{m}) =  -1 *  r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m}) \f$
///
/// \f$ \hat{\delta'}_{z}(r_{i},z_{j+1},\phi_{m}) =  -1 *  \hat{\delta}_{z}(r_{i},z_{j},\phi_{m}) \f$
///
/// ~~~
///	(*corrDrDz)(i,j+1) 		= -1* (*distDrDz)(i,j) ;				
/// (*corrDphiRDz)(i,j+1) = -1* (*distDphiRDz)(i,j);
/// (*corrDz)(i,j+1)      = -1* (*distDz)(i,j);				
/// ~~~
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesDistDrDz TMatrixD**  local distortion \f$\hat{\delta}_{r}\f$
///	\param matricesDistDphiRDz TMatrixD** local distortion \f$r \hat{\delta}_{\phi}\f$
///	\param matricesDistDz TMatrixD**   local distortion \f$ \hat{\delta}_{z}\f$
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez assummed already been calculated
/// \post Local distortion and correction are computed according simplified Langevin equation
/// ~~~ 
/// matricesDistDrDz,matricesDistDphiRDz,matricesDistDz
/// ~~~ 
/// and correction: 
/// ~~~ 
/// matricesCorrDrDz,matricesCorrDphiRDz,matricesCorrDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::LocalDistCorrDz
( 
 TMatrixD** matricesEr, 	
 TMatrixD** matricesEphi, 	
 TMatrixD** matricesEz, 	
 TMatrixD** matricesDistDrDz, 
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 		
 TMatrixD** matricesCorrDrDz, 
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 TMatrixD** matricesLocalIntErDz, 
 TMatrixD** matricesLocalIntEphiDz, 
 TMatrixD** matricesLocalIntEz,		
 const Int_t rrow, 
 const Int_t zcolumn, 
 const Int_t phiSlice,  
 const Float_t gridSizeZ,
 const Double_t ezField
  )
{
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEphiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;
	
  // pointer declaration
  TMatrixD * eR;
  TMatrixD * ePhi;
  TMatrixD * eZ;
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
	
  TMatrixD * localIntErDz;
  TMatrixD * localIntEphiDz;
  TMatrixD * localIntEz;
	
	
  // Initialization for j == colomn-1 integration is 0.0
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
		
	
		
    for (Int_t i=0; i< rrow;i++) {
      (*distDrDz)(i,zcolumn-1) =  0.0 ;			
      (*distDphiRDz)(i,zcolumn-1) =  0.0 ;			
      (*distDz)(i,zcolumn-1) =  0.0 ;			
			
      (*corrDrDz)(i,0) =  0.0 ;			
      (*corrDphiRDz)(i,0) =  0.0 ;			
      (*corrDz)(i,0) =  0.0 ;			
    }			
  }
	
  // for this case
  // use trapezoidal rule assume no ROC displacement
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    eR  =  matricesEr[m] ;
    ePhi  =  matricesEphi[m] ;
    eZ  =  matricesEz[m] ;
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
		
    localIntErDz  =  matricesLocalIntErDz[m] ; 
    localIntEphiDz =  matricesLocalIntEphiDz[m] ; 
    localIntEz =  matricesLocalIntEz[m] ; ;
    for (Int_t j=0; j< zcolumn-1;j++) {
      for (Int_t i=0; i< rrow;i++)  {
	localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1)) / (-1*ezField) ;
	localIntEphiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1)) / (-1*ezField) ;				
	localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;				
				
	(*localIntErDz)(i,j) = localIntErOverEz;
	(*localIntEphiDz)(i,j) = localIntEphiOverEz;
	(*localIntEz)(i,j) = localIntDeltaEz;
				
				
				
				
	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz;
	(*distDz)(i,j)      = localIntDeltaEz*o2::TPC::AliTPCPoissonSolver::fgkdvdE*o2::TPC::AliTPCPoissonSolver::fgkdvdE;// two times?				
				
				
	(*corrDrDz)(i,j+1) 		= -1* (*distDrDz)(i,j) ;				
	(*corrDphiRDz)(i,j+1) = -1* (*distDphiRDz)(i,j);
	(*corrDz)(i,j+1)      = -1* (*distDz)(i,j);							
				
      }
    }				
  }		
}


///
/// Local Distortion only
/// See explanation at LocalDistCorrDz
///
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesDistDrDz TMatrixD**  local distortion \f$\hat{\delta}_{r}\f$
///	\param matricesDistDphiRDz TMatrixD** local distortion \f$r \hat{\delta}_{\phi}\f$
///	\param matricesDistDz TMatrixD**   local distortion \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
/// \param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez assummed already been calculated
/// \post Local distortion are computed according simplified Langevin equation
/// ~~~ 
/// matricesDistDrDz,matricesDistDphiRDz,matricesDistDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::LocalDistDz
(
 TMatrixD** matricesEr, 	
 TMatrixD** matricesEphi, 	
 TMatrixD** matricesEz, 	
 TMatrixD** matricesDistDrDz, 
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 
 const Int_t rrow, 
 const Int_t zcolumn, 
 const Int_t phiSlice,  
 const Float_t gridSizeZ,
 const Double_t ezField
 )
{
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEphiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;
	
  // pointer declaration
  TMatrixD * eR;
  TMatrixD * ePhi;
  TMatrixD * eZ;
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;
		
	
	
  // Initialization for j == colomn-1 integration is 0.0
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
    for (Int_t i=0; i< rrow;i++) {
      (*distDrDz)(i,0) =  0.0 ;			
      (*distDphiRDz)(i,0) =  0.0 ;			
      (*distDz)(i,0) =  0.0 ;			
    }
			
  }
	
  // for this case
  // use trapezoidal rule assume no ROC displacement
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    eR  =  matricesEr[m] ;
    ePhi  =  matricesEphi[m] ;
    eZ  =  matricesEz[m] ;
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    for (Int_t j=1; j< zcolumn;j++) {
      for (Int_t i=0; i< rrow;i++)  {
				 
	localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j-1)+(*eR)(i,j))/(-1*ezField) ;
	localIntEphiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j-1)+(*ePhi)(i,j))/(-1*ezField) ;
	localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j-1)+(*eZ)(i,j)) ;				
				
	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
	(*distDz)(i,j)      = localIntDeltaEz*o2::TPC::AliTPCPoissonSolver::fgkdvdE*o2::TPC::AliTPCPoissonSolver::fgkdvdE; // two times?				
      }
    }				
  }		
}


///
/// Local Correction only
/// See explanation at LocalDistCorrDz
///
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez are provided
/// \post Local correction are computed according simplified Langevin equation
/// ~~~ 
/// matricesCorrDz,matricesCorrDphiRDz,matricesDistDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::LocalCorrDz
(
 TMatrixD** matricesEr, 	
 TMatrixD** matricesEphi, 	
 TMatrixD** matricesEz, 	
 TMatrixD** matricesCorrDrDz, 
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 const Int_t rrow, 
 const Int_t zcolumn, 
 const Int_t phiSlice,  
 const Float_t gridSizeZ,
 const Double_t ezField
 )
{
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEphiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;
	
  // pointer declaration
  TMatrixD * eR;
  TMatrixD * ePhi;
  TMatrixD * eZ;
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
		
	
	
  // Initialization for j == colomn-1 integration is 0.0
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
    for (Int_t i=0; i< rrow;i++) {
      (*corrDrDz)(i,zcolumn-1) =  0.0 ;			
      (*corrDphiRDz)(i,zcolumn-1) =  0.0 ;			
      (*corrDz)(i,zcolumn-1) =  0.0 ;			
    }
			
  }
	
  // for this case
  // use trapezoidal rule assume no ROC displacement
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    eR  =  matricesEr[m] ;
    ePhi  =  matricesEphi[m] ;
    eZ  =  matricesEz[m] ;
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
    for (Int_t j=0; j< zcolumn-1;j++) {
      for (Int_t i=0; i< rrow;i++)  {
				
	localIntErOverEz = -1 * (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1))/(-1*ezField) ;
	localIntEphiOverEz = -1 *  (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1))/(-1*ezField) ;
	localIntDeltaEz = -1 * (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;				
	// Scale the Ez distortions with the drift velocity pertubation -> delivers cm
	// from langevin eq (local distortion)
	(*corrDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	(*corrDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
	(*corrDz)(i,j)      = localIntDeltaEz*o2::TPC::AliTPCPoissonSolver::fgkdvdE*o2::TPC::AliTPCPoissonSolver::fgkdvdE; //two times?
				
      }
    }		
    //corrDrDz->Print();
  }		
}

/// Global correction and distortion (Naive algorithm)
/// 
///
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez are provided
/// \post Local correction are computed according simplified Langevin equation
/// ~~~ 
/// matricesCorrDz,matricesCorrDphiRDz,matricesDistDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDz 
(
 TMatrixD** matricesDistDrDz,  
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 
 TMatrixD** matricesGDistDrDz,  
 TMatrixD** matricesGDistDphiRDz, 
 TMatrixD** matricesGDistDz, 
 TMatrixD** matricesCorrDrDz,  
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 TMatrixD** matricesGCorrDrDz,  
 TMatrixD** matricesGCorrDphiRDz, 
 TMatrixD** matricesGCorrDz, 
		
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
	
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0, radius,phi,z,radiusc;
  radiusc = 0.0;
  radius = 0.0;
  TMatrixD * mDistDrDz;
  TMatrixD * mDistDphiRDz;
  TMatrixD * mDistDz;
	
  TMatrixD * mCorrDrDz;
  TMatrixD * mCorrDphiRDz;
  TMatrixD * mCorrDz;
	

  TMatrixD * mCorrIrregularDrDz;
  TMatrixD * mCorrIrregularDphiRDz;
  TMatrixD * mCorrIrregularDz;

  TMatrixD * mRDistortedPoint;
  TMatrixD * mPhiDistortedPoint;
  TMatrixD * mZDistortedPoint;

			
	
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    phi0 = philist[m];
				
    mDistDrDz = matricesGDistDrDz[m];
    mDistDphiRDz = matricesGDistDphiRDz[m];
    mDistDz = matricesGDistDz[m];
    mCorrDrDz = matricesGCorrDrDz[m];
    mCorrDphiRDz = matricesGCorrDphiRDz[m];
    mCorrDz = matricesGCorrDz[m];


    /// irregular grid
    mCorrIrregularDrDz;
    mCorrIrregularDphiRDz;
    mCorrIrregularDz;

    mRDistortedPoint;
    mPhiDistortedPoint;
    mZDistortedPoint;
	



    for (Int_t j=zcolumn-1; j >= 0;j--) {
      z0 = zlist[j] ;		
      printf("global dist j:%d\n",j);
			
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	phi = phi0;
	radius = radius0;
				
	dr = 0.0;
	dphir = 0.0;
	dz = 0.0;
	ddrphi = 0.0;
	for (Int_t jj = j; jj < zcolumn;jj++) {
	  // interpolation the local distortion for current position
	  phi += ddrphi/radius;
	  radius = radius0 + dr;
	  z = zlist[jj] + dz;
					
	  if (phi < 0.0) phi = TMath::TwoPi() + phi;
	  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

	  ddr    = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
					 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDrDz);
	  ddrphi = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
					 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDphiRDz);
	  ddz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDz);
						
	  dr += ddr;
	  dphir += ddrphi;
	  dz += ddz;
				
					
	}
				
	(*mDistDrDz)(i,j) = dr;
	(*mDistDphiRDz)(i,j) = dphir;
	(*mDistDz)(i,j) = dz;

	//////////////// use irregular grid look up table for correction







	///////////////				

	if (j == 1) radiusc = radius0;
				
	dr = (*mCorrDrDz)(i,j-1);
	dphir = (*mCorrDphiRDz)(i,j-1);
	dz = (*mCorrDz)(i,j-1) ;
				
				
				
	radiusc = radius0 + dr;				
	phi = phi0 + dphir/radiusc;
	z = zlist[j - 1] + dz;
				
					
					
	ddr    = Interpolate3DTableCyl(fInterpolationOrder, radiusc, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDrDz);
	ddrphi = Interpolate3DTableCyl(fInterpolationOrder, radiusc, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDphiRDz);
	ddz		 = Interpolate3DTableCyl(fInterpolationOrder, radiusc, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDz);
						
	dr += ddr;
	dz += ddz;
	dphir += ddrphi;
					
				
	(*mCorrDrDz)(i,j) = dr;
	(*mCorrDphiRDz)(i,j) = dphir;
	(*mCorrDz)(i,j) = dz;
				
				
      }
    }
  }
}


/// Local Correction only
/// See explanation at LocalDistCorrDz
///
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez are provided
/// \post Local correction are computed according simplified Langevin equation
/// ~~~ 
/// matricesCorrDz,matricesCorrDphiRDz,matricesDistDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDz 
(
 AliTPCLookUpTable3DInterpolatorD *lookupLocalDist, 
 TMatrixD** matricesGDistDrDz,  
 TMatrixD** matricesGDistDphiRDz, 
 TMatrixD** matricesGDistDz,
 AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr, 
 TMatrixD** matricesGCorrDrDz,  
 TMatrixD** matricesGCorrDphiRDz, 
 TMatrixD** matricesGCorrDz, 


 TMatrixD** matricesGCorrIrregularDrDz,  
 TMatrixD** matricesGCorrIrregularDphiRDz, 
 TMatrixD** matricesGCorrIrregularDz, 

 TMatrixD** matricesRIrregular,  
 TMatrixD** matricesPhiIrregular, 
 TMatrixD** matricesZIrregular, 
		 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
		
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0, radius,phi,z,radiusc;
  radiusc = 0.0;
  radius = 0.0;
  TMatrixD * mDistDrDz;
  TMatrixD * mDistDphiRDz;
  TMatrixD * mDistDz;	
  TMatrixD * mCorrDrDz;
  TMatrixD * mCorrDphiRDz;
  TMatrixD * mCorrDz;

  TMatrixD * mCorrIrregularDrDz;
  TMatrixD * mCorrIrregularDphiRDz;
  TMatrixD * mCorrIrregularDz;

  TMatrixD * mRIrregular;
  TMatrixD * mPhiIrregular;
  TMatrixD * mZIrregular;
	
  // initialized for zcolumn - 1
  Int_t j = zcolumn  - 1;
  z0 = zlist[j] ;		
	
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    phi0 = philist[m];
			
    mDistDrDz = matricesGDistDrDz[m];
    mDistDphiRDz = matricesGDistDphiRDz[m];
    mDistDz = matricesGDistDz[m];
	
    // 
    mCorrDrDz = matricesGCorrDrDz[m];
    mCorrDphiRDz = matricesGCorrDphiRDz[m];
    mCorrDz = matricesGCorrDz[m];


    mCorrIrregularDrDz  = matricesGCorrIrregularDrDz[m];
    mCorrIrregularDphiRDz = matricesGCorrIrregularDphiRDz[m];
    mCorrIrregularDz = matricesGCorrIrregularDz[m];

    mRIrregular = matricesRIrregular[m];
    mPhiIrregular = matricesPhiIrregular[m];
    mZIrregular= matricesZIrregular[m];

		
    for (Int_t i=0; i< rrow;i++)  {
      // do from j to 0
      // follow the drift
      radius0 = rlist[i];
      phi = phi0;
      radius = radius0;
			
      dr = 0.0;
      dphir = 0.0;
      dz = 0.0;
      ddrphi = 0.0;
			
				
				
				
      (*mDistDrDz)(i,j) = dr;
      (*mDistDphiRDz)(i,j) = dphir;
      (*mDistDz)(i,j) = dz;




      //////////////// use irregular grid look up table for correction
      // set 
      (*mCorrIrregularDrDz)(i,j) = -dr;
      (*mCorrIrregularDphiRDz)(i,j) = -dphir;
      (*mCorrIrregularDz)(i,j) = -dz;


      // distorted point
      (*mRIrregular)(i,j) = radius0 + dr;
      (*mPhiIrregular)(i,j) =	phi0 + (dphir/radius0);
      (*mZIrregular)(i,j) = z0 + dz;
				


      ///////////////				

		
				
    }
  }
	

  for (j=zcolumn-2; j>=0;j--) {
    z0 = zlist[j] ;		
    //printf("global dist j:%d\n",j);
			
	
    for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
      phi0 = philist[m];
				
      mDistDrDz = matricesGDistDrDz[m];
      mDistDphiRDz = matricesGDistDphiRDz[m];
      mDistDz = matricesGDistDz[m];
	
      // 
      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDphiRDz = matricesGCorrDphiRDz[m];
      mCorrDz = matricesGCorrDz[m];


      mCorrIrregularDrDz  = matricesGCorrIrregularDrDz[m];
      mCorrIrregularDphiRDz = matricesGCorrIrregularDphiRDz[m];
      mCorrIrregularDz = matricesGCorrIrregularDz[m];

      mRIrregular = matricesRIrregular[m];
      mPhiIrregular = matricesPhiIrregular[m];
      mZIrregular= matricesZIrregular[m];

			
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	phi = phi0;
	radius = radius0;
				
	dr = 0.0;
	dphir = 0.0;
	dz = 0.0;
	ddrphi = 0.0;
				
				
				
	for (Int_t jj = j; jj < zcolumn;jj++) {
	  // interpolation the local distortion for current position
	  phi += ddrphi/radius;					
	  radius = radius0 + dr;
	  z = zlist[jj] + dz;
					
	  if (phi < 0.0) phi = TMath::TwoPi() + phi;
	  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;
					
	  //  printf("lookupdist (%f,%f,%f)\n",radius,phi,z);

	  lookupLocalDist->GetValue(radius,phi,z,ddr,ddrphi,ddz);

	  //printf("lookdist     (%f,%f,%f)\n",ddr,ddrphi,dz);
					
					
	  dr += ddr;
	  dphir += ddrphi;
	  dz += ddz;
					
				
					
	}
				
	(*mDistDrDz)(i,j) = dr;
	(*mDistDphiRDz)(i,j) = dphir;
	(*mDistDz)(i,j) = dz;




	//////////////// use irregular grid look up table for correction
	// set 
	(*mCorrIrregularDrDz)(i,j) = -dr;
	(*mCorrIrregularDphiRDz)(i,j) = -dphir;
	(*mCorrIrregularDz)(i,j) = -dz;


	// distorted point
	(*mRIrregular)(i,j) = radius0 + dr;
	(*mPhiIrregular)(i,j) =	phi0 + (dphir/radius0);
	(*mZIrregular)(i,j) = z0 + dz;
				


	///////////////				

		
				
	if (j == zcolumn - 2) radiusc = radius0;
				
	dr = (*mCorrDrDz)(i,j+1);
	dphir = (*mCorrDphiRDz)(i,j+1);
	dz = (*mCorrDz)(i,j+1) ;
				
				
	phi = phi0 + dphir/radiusc;
	radiusc = radius0 + dr;				
	z = zlist[j + 1] + dz;
				
	//if (phi < 0.0) printf("lookupcorr (%f,%f,%f)\n",radiusc,phi,z);
	while (phi < 0.0) phi = TMath::TwoPi() + phi;
	while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;
				


				
	lookupLocalCorr->GetValue(radiusc,phi,z,ddr,ddrphi,ddz);	
	//if (phi < 0.0) printf("lookup (%f,%f,%f)\n",ddr,ddrphi,ddz);
				
	//if (radiusc > 250.0 //printf("lookcorr     (%f,%f,%f)\n",ddr,ddrphi,dz);
						
	dr += ddr;
	dz += ddz;
	dphir += ddrphi;
					
				
	(*mCorrDrDz)(i,j) = dr;
	(*mCorrDphiRDz)(i,j) = dphir;
	(*mCorrDz)(i,j) = dz;
				
				
      }
    }
  }
}


/// Local Correction only
/// See explanation at LocalDistCorrDz
///
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component 
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component 
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param rrow Int_t Number of rows in r-direction
/// \param zcolumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction 
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez are provided
/// \post Local correction are computed according simplified Langevin equation
/// ~~~ 
/// matricesCorrDz,matricesCorrDphiRDz,matricesDistDz
/// ~~~ 
/// 
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDzOpt2 
(
 TMatrixD** matricesDistDrDz,  
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 
 TMatrixD** matricesGDistDrDz,  
 TMatrixD** matricesGDistDphiRDz, 
 TMatrixD** matricesGDistDz, 
 TMatrixD** matricesCorrDrDz,  
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 TMatrixD** matricesGCorrDrDz,  
 TMatrixD** matricesGCorrDphiRDz, 
 TMatrixD** matricesGCorrDz, 
		
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0, radius,phi,z;
	
  radius   = 0.0;
	
  TMatrixD * mDistDrDz;
  TMatrixD * mDistDphiRDz;
  TMatrixD * mDistDz;
	
  TMatrixD * mCorrDrDz;
  TMatrixD * mCorrDphiRDz;
  TMatrixD * mCorrDz;
	

  dr = 0.0;
  dz = 0.0;
  dphir = 0.0;
		
  for (Int_t j=zcolumn-2; j >= 0;j--) {
    z0 = zlist[j] ;		
			
	
    for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
      phi0 = philist[m];
				
      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDphiRDz = matricesGCorrDphiRDz[m];
      mCorrDz = matricesGCorrDz[m];;
			
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	//if (j == zcolumn-2) radius = radius0;
				
	dr = (*mCorrDrDz)(i,j+1);
	dphir = (*mCorrDphiRDz)(i,j+1);
	dz = (*mCorrDz)(i,j+1) ;
				
				
	radius = radius0 + dr;										
	phi = phi0 + dphir/radius;					
	z = zlist[j + 1] + dz;
				
					
					
	ddr    = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDrDz);
	ddrphi = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDphiRDz);
	ddz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDz);
						
	dr += ddr;
	dz += ddz;
	dphir += ddrphi;
					
				
	(*mCorrDrDz)(i,j) = dr;
	(*mCorrDphiRDz)(i,j) = dphir;
	(*mCorrDz)(i,j) = dz;
				
				
      }
    }
  }
	
	
  dr = 0.0;
  dz = 0.0;
  dphir = 0.0;
	
  for (Int_t j=zcolumn-2; j >= 0;j--) {
    z0 = zlist[j] ;		
			
    for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
      phi0 = philist[m];		
      mDistDrDz = matricesGDistDrDz[m];
      mDistDphiRDz = matricesGDistDphiRDz[m];
      mDistDz = matricesGDistDz[m];
		
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	phi = phi0;
	radius = radius0;
	z = z0;
				
	ddr    = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDrDz);
	ddrphi = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDphiRDz);
	ddz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDz);
				
	phi += ddrphi/radius;
	radius = radius0 + dr;
	z = zlist[j+1] + dz;
				
	dr    =  Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesGDistDrDz,j+1);
	dphir =  Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesGDistDphiRDz,j+1);
	dz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesGDistDz,j+1);
					
				
				
	(*mDistDrDz)(i,j) = dr + ddr;
	(*mDistDphiRDz)(i,j) = dphir + ddrphi;
	(*mDistDz)(i,j) = dz  + ddz;
				
				
      }
    }
  }
	
	
}





void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDzOpt2 
(
 AliTPCLookUpTable3DInterpolatorD *lookupLocalDist, 
 TMatrixD** matricesGDistDrDz,  
 TMatrixD** matricesGDistDphiRDz, 
 TMatrixD** matricesGDistDz, 
 AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr, 
 TMatrixD** matricesGCorrDrDz,  
 TMatrixD** matricesGCorrDphiRDz, 
 TMatrixD** matricesGCorrDz, 
		
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0, radius,phi,z;
	
  radius   = 0.0;
	
  TMatrixD * mDistDrDz;
  TMatrixD * mDistDphiRDz;
  TMatrixD * mDistDz;
	
  TMatrixD * mCorrDrDz;
  TMatrixD * mCorrDphiRDz;
  TMatrixD * mCorrDz;
	

  dr = 0.0;
  dz = 0.0;
  dphir = 0.0;
		
  for (Int_t j=zcolumn-2; j >= 0;j--) {
    z0 = zlist[j] ;		
			
	
    for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
      phi0 = philist[m];
				
      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDphiRDz = matricesGCorrDphiRDz[m];
      mCorrDz = matricesGCorrDz[m];;
			
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	//if (j == zcolumn-2) radius = radius0;
				
	dr = (*mCorrDrDz)(i,j+1);
	dphir = (*mCorrDphiRDz)(i,j+1);
	dz = (*mCorrDz)(i,j+1) ;
				
				
	radius = radius0 + dr;		
	phi = phi0 + dphir/radius;					
	z = zlist[j + 1] + dz;
				
					
					
	lookupLocalCorr->GetValue(radius,phi,z,ddr,ddrphi,ddz);	
						
	dr += ddr;
	dz += ddz;
	dphir += ddrphi;
					
				
	(*mCorrDrDz)(i,j) = dr;
	(*mCorrDphiRDz)(i,j) = dphir;
	(*mCorrDz)(i,j) = dz;
				
				
      }
    }
  }
	
	
  dr = 0.0;
  dz = 0.0;
  dphir = 0.0;
		
	
  for (Int_t j=zcolumn-2; j >= 0;j--) {
    z0 = zlist[j] ;		
			
    for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
      phi0 = philist[m];		
      mDistDrDz = matricesGDistDrDz[m];
      mDistDphiRDz = matricesGDistDphiRDz[m];
      mDistDz = matricesGDistDz[m];
		
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	phi = phi0;
	radius = radius0;
	z = z0;				
				
	lookupLocalDist->GetValue(radius,phi,z,ddr,ddrphi,ddz);	
				
	phi += ddrphi/radius;
	radius = radius0 + ddr;
	z = zlist[j+1] + ddz;
				
	dr    =  Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesGDistDrDz,j+1);
	dphir =  Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
				       rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesGDistDphiRDz,j+1);
	dz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesGDistDz,j+1);
					
				
				
	(*mDistDrDz)(i,j) = dr + ddr;
	(*mDistDphiRDz)(i,j) = dphir + ddrphi;
	(*mDistDz)(i,j) = dz  + ddz;
				
				
      }
    }
  }
	
	
}


void AliTPCSpaceCharge3DDriftLine::IntegrateDistDriftLineDz 
(
 TMatrixD** matricesDistDrDz,  
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 
 TMatrixD** matricesGDistDrDz,  
 TMatrixD** matricesGDistDphiRDz, 
 TMatrixD** matricesGDistDz, 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
	
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0, radius,phi,z;
  TMatrixD * mDistDrDz;
  TMatrixD * mDistDphiRDz;
  TMatrixD * mDistDz;
	
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    mDistDrDz = matricesGDistDrDz[m];
    mDistDphiRDz = matricesGDistDphiRDz[m];
    mDistDz = matricesGDistDz[m];
		
	
    for (Int_t j=1; j< zcolumn;j++) {
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	phi0 = philist[m];
	z0 = zlist[j] ;		
	radius0 = rlist[i];
	phi = phi0;
	radius = radius0;
				
	dr = 0.0;
	dphir = 0.0;
	dz = 0.0;
	ddrphi = 0.0;
	for (Int_t jj = j; jj > 0;jj--) {
	  // interpolation the local distortion for current position
	  phi += ddrphi/radius;
	  radius = radius0 + dr;
	  z = zlist[jj] + dz;
					
	  if (phi < 0.0) phi = TMath::TwoPi() + phi;
	  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;
				
	  ddr    = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
					 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDrDz);
	  ddrphi = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
					 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDphiRDz);
	  ddz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesDistDz);
						
	  dr += ddr;
	  dphir += ddrphi;
	  dz += ddz;
				
					
	}
				
	(*mDistDrDz)(i,j) = dr;
	(*mDistDphiRDz)(i,j) = dphir;
	(*mDistDz)(i,j) = dz;
				
				
      }
    }
  }
}



void AliTPCSpaceCharge3DDriftLine::IntegrateCorrDriftLineDz 
(
 TMatrixD** matricesCorrDrDz,  
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 TMatrixD** matricesGCorrDrDz,  
 TMatrixD** matricesGCorrDphiRDz, 
 TMatrixD** matricesGCorrDz, 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
	
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0,radius,phi,z;
	
  TMatrixD * mCorrDrDz;
  TMatrixD * mCorrDphiRDz;
  TMatrixD * mCorrDz;
	
	
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    mCorrDrDz = matricesGCorrDrDz[m];
    mCorrDphiRDz = matricesGCorrDphiRDz[m];
    mCorrDz = matricesGCorrDz[m];;
	
    for (Int_t j=0; j< zcolumn;j++) {
			
      for (Int_t i=0; i< rrow;i++)  {
	// do from 0 to j
	// follow the drift
	phi0 = philist[m];
	z0 = zlist[j] ;		
	radius0 = rlist[i];
	radius = radius0;
	phi = phi0;
				
	dr = 0.0;
	dphir = 0.0;
	dz = 0.0;
	ddrphi = 0.0;
	for (Int_t jj = 0; jj < j;jj++) {
					
	  phi += ddrphi/radius;
	  radius = radius0 + dr;
	  z = zlist[jj] + dz;
					
	  if (phi < 0.0) phi = TMath::TwoPi() + phi;
	  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

					
	  ddr    = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
					 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDrDz);
	  ddrphi = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
					 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDphiRDz);
	  ddz		 = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi, 
						 rrow, zcolumn, phiSlice,  rlist, zlist, philist, matricesCorrDz);
						
	  dr += ddr;
	  dphir += ddrphi;
	  dz += ddz;					
					
					
	}
	// should be fix
	(*mCorrDrDz)(i,j) = dr;
	(*mCorrDphiRDz)(i,j) = dphir;
	(*mCorrDz)(i,j) = dz;
      }
    }
    //if (m == 0) (*mCorrDrDz).Print();	
  }
	
	
}


void AliTPCSpaceCharge3DDriftLine::FillLookUpTable
(
 TMatrixD** matricesGDrDz,  
 TMatrixD** matricesGDphiRDz, 
 TMatrixD** matricesGDz, 
 TMatrixD** lookupRDz,  
 TMatrixD** lookupPhiRDz, 
 TMatrixD** lookupDz, 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist,
 Int_t side
 )
{
  Double_t  r, phi, z ;
  TMatrixD *mR  ;
  TMatrixD *mPhiR ;
  TMatrixD *mDz ;

  /// * Interpolate basicLookup tables; once for each rod, then sum the results			
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    phi = fLookupPhiList[k] ;
    mR   =    lookupRDz[k]  ;
    mPhiR =   lookupPhiRDz[k];
    mDz    =  lookupDz[k];
    for ( Int_t j = 0 ; j < fNZColumns ; j++ ) {
      z = TMath::Abs(fLookupZList[j]) ;  // Symmetric solution in Z that depends only on ABS(Z)
      if ( side == 0 &&  fLookupZList[j] < 0 ) continue; // Skip rest of this loop if on the wrong side
      if ( side == 1 &&  fLookupZList[j] > 0 ) continue; // Skip rest of this loop if on the wrong side
      for ( Int_t i = 0 ; i < fNRRows ; i++ ) {	
	r = fLookupRList[i] ;
	// Interpolate basicLookup tables; once for each rod, then sum the results
	(*mR)(i,j)   		= Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, rrow, zcolumn, phiSlice,
							rlist, zlist, philist, matricesGDrDz );
	(*mPhiR)(i,j) 	= Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, rrow, zcolumn, phiSlice,
						rlist, zlist, philist, matricesGDphiRDz);
	(*mDz)(i,j)    	= Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, rrow, zcolumn, phiSlice,
						rlist, zlist, philist, matricesGDz);
					
				
	if (side == 1)  (*mDz)(i,j) = -  (*mDz)(i,j); // negative coordinate system on C side
      }
    }
  }
}



void AliTPCSpaceCharge3DDriftLine::FillLookUpTable
(
 AliTPCLookUpTable3DInterpolatorD *lookupGlobal, 
 TMatrixD** lookupRDz,  
 TMatrixD** lookupPhiRDz, 
 TMatrixD** lookupDz, 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist,
 Int_t side
 )
{
  Double_t  r, phi, z ;
  TMatrixD *mR  ;
  TMatrixD *mPhiR ;
  TMatrixD *mDz ;

  /// * Interpolate basicLookup tables; once for each rod, then sum the results			
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    phi = fLookupPhiList[k] ;
    mR   =    lookupRDz[k]  ;
    mPhiR =   lookupPhiRDz[k];
    mDz    =  lookupDz[k];
    for ( Int_t j = 0 ; j < fNZColumns ; j++ ) {
      z = TMath::Abs(fLookupZList[j]) ;  // Symmetric solution in Z that depends only on ABS(Z)
      if ( side == 0 &&  fLookupZList[j] < 0 ) continue; // Skip rest of this loop if on the wrong side
      if ( side == 1 &&  fLookupZList[j] > 0 ) continue; // Skip rest of this loop if on the wrong side
      for ( Int_t i = 0 ; i < fNRRows ; i++ ) {	
	r = fLookupRList[i] ;
					
	lookupGlobal->GetValue(r, phi, z, (*mR)(i,j) ,(*mPhiR)(i,j),(*mDz)(i,j));
					
	if (side == 1)  (*mDz)(i,j) = -  (*mDz)(i,j); // negative coordinate system on C side
      }
    }
  }
}

void AliTPCSpaceCharge3DDriftLine::FillLookUpTableA
(
 AliTPCLookUpTable3DInterpolatorD *lookupGlobal, 
 TMatrixD** lookupRDz,  
 TMatrixD** lookupPhiRDz, 
 TMatrixD** lookupDz, 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
  Double_t  r, phi, z ,rl,phil,zl;
  TMatrixD *mR  ;
  TMatrixD *mPhiR ;
  TMatrixD *mDz ;

  /// * Interpolate basicLookup tables; once for each rod, then sum the results			
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    phi = fLookupPhiList[k] ;
    phil = philist[k];
		
    mR   =    lookupRDz[k]  ;
    mPhiR =   lookupPhiRDz[k];
    mDz    =  lookupDz[k];
    for ( Int_t j = 0 ; j < fNZColumns/2 ; j++ ) {
      z = fLookupZListA[j] ;  // Symmetric solution in Z that depends only on ABS(Z)
			
      zl = zlist[j];
      //			printf("(%f,%f)(%f,%f)\n",z,zl);
			
      for ( Int_t i = 0 ; i < fNRRows ; i++ ) {	
	r = fLookupRList[i] ;
	rl = rlist[i];
				
						
	lookupGlobal->GetValue(r, phi, z, (*mR)(i,j) ,(*mPhiR)(i,j),(*mDz)(i,j));
			
      }
    }
  }
}




void AliTPCSpaceCharge3DDriftLine::FillLookUpTableC
(
 AliTPCLookUpTable3DInterpolatorD *lookupGlobal, 
 TMatrixD** lookupRDz,  
 TMatrixD** lookupPhiRDz, 
 TMatrixD** lookupDz, 
 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{
  Double_t  r, phi, z ;
  TMatrixD *mR  ;
  TMatrixD *mPhiR ;
  TMatrixD *mDz ;

  /// * Interpolate basicLookup tables; once for each rod, then sum the results			
  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    phi = fLookupPhiList[k] ;
    mR   =    lookupRDz[k]  ;
    mPhiR =   lookupPhiRDz[k];
    mDz    =  lookupDz[k];
    for ( Int_t j = 0 ; j < fNZColumns/2 ; j++ ) {
      z = TMath::Abs(fLookupZListC[j]) ;  // Symmetric solution in Z that depends only on ABS(Z)
			
      for ( Int_t i = 0 ; i < fNRRows ; i++ ) {	
	r = fLookupRList[i] ;
					
	lookupGlobal->GetValue(r, phi, z, (*mR)(i,j) ,(*mPhiR)(i,j),(*mDz)(i,j));
	(*mDz)(i,j) = -(*mDz)(i,j);
      }
    }
  }
}

void AliTPCSpaceCharge3DDriftLine::GetDistortionCyl(const Float_t x[], Short_t roc,Float_t dx[]) {
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
		
  GetDistCorrFromLookUpTableCyl(x,roc,dx,fLookupIntDist);
	
}



void AliTPCSpaceCharge3DDriftLine::GetDistortionCylAC(const Float_t x[], Short_t roc,Float_t dx[]) {
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
	
	
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
	
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;
    
  if (z > -1e-6) 
    fLookupIntDistA->GetValue(r,phi,z,dR,dPhiR,dZ);
  else 
    fLookupIntDistC->GetValue(r,phi,-z,dR,dPhiR,dZ);
	
		
	
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)
	
}



// Get Correction from irreggular table
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylAC(const Float_t x[], Short_t roc,Float_t dx[]) {
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
	
	
  Double_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;


  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (fNRRows - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (fNZColumns/2 -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / fNPhiSlices;


  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
	
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;


  // get distorion from irregular table


  Int_t ianchor = TMath::FloorNint((r - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR); 
  Int_t kanchor = TMath::FloorNint(phi/gridSizePhi); 
  Int_t zanchor = TMath::FloorNint(z/gridSizeZ);

  if (z > -1e-6) 
    fLookupIntCorrIrregularA->GetValue(r,phi,z,dR,dPhiR,dZ,ianchor,kanchor,zanchor,fNRRows/8  + 1, fNPhiSlices/8 + 1,fNZColumns/8 + 1,0);
  else 
    fLookupIntCorrIrregularC->GetValue(r,phi,z,dR,dPhiR,dZ,ianchor,kanchor,zanchor,fNRRows/8  + 1, fNPhiSlices/8 + 1,fNZColumns/8 + 1,0);
	
		
	
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)
	
}


void AliTPCSpaceCharge3DDriftLine::GetDistortion(const Float_t x[], Short_t roc,Float_t dx[]) {
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
  //printf("calling get distortion spacecharge3dintegraldz\n");
  //GetDistCorrFromLookUpTable(x,roc,dx,fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz);
  GetDistCorrFromLookUpTable(x,roc,dx,fLookupIntDist);
}


void AliTPCSpaceCharge3DDriftLine::GetCorrectionCyl(const Float_t x[], Short_t roc,Float_t dx[]) 
{
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
  //printf("calling get correction spacecharge3dintegraldz\n");
	
	
  //GetDistCorrFromLookUpTableCyl(x,roc,dx,fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz);
  GetDistCorrFromLookUpTableCyl(x,roc,dx,fLookupIntCorr);
}

void AliTPCSpaceCharge3DDriftLine::GetCorrection(const Float_t x[], Short_t roc,Float_t dx[]) 
{
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
  //printf("calling get correction spacecharge3dintegraldz\n");
	
  //GetDistCorrFromLookUpTable(x,roc,dx,fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz);
  GetDistCorrFromLookUpTable(x,roc,dx,fLookupIntCorr);

}


void AliTPCSpaceCharge3DDriftLine::GetDistCorrFromLookUpTable
(
 const Float_t x[], 
 Short_t roc,
 Float_t dx[],
 TMatrixD**lookUpR,
 TMatrixD**lookUpPhiR,
 TMatrixD**lookUpDz
 ) 
{
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  TMath::Sqrt( x[0]*x[0] + x[1]*x[1] ) ;
  phi    =  TMath::ATan2(x[1],x[0]) ;
  
  while (phi>TMath::Pi()) phi-=TMath::TwoPi();
  while (phi<-TMath::Pi()) phi+=TMath::TwoPi();

  
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // remove
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;


	
  // Get distortion
  dR      = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
				  fLookupRList, fLookupZList, fLookupPhiList, lookUpR );
  dPhiR    = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
				   fLookupRList, fLookupZList, fLookupPhiList, lookUpPhiR);
  dZ = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
			     fLookupRList, fLookupZList, fLookupPhiList, lookUpDz   );

  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + fCorrectionFactor * dPhiR / r;
    r   =  r   + fCorrectionFactor * dR;
  }
  dZ =  fCorrectionFactor * dZ ;

  //printf("correction factor %f\n",fCorrectionFactor);
  // Calculate correction in cartesian coordinates
  dx[0] = (r * TMath::Cos(phi) - x[0]);
  dx[1] = (r * TMath::Sin(phi) - x[1]);
  dx[2] = dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}


void AliTPCSpaceCharge3DDriftLine::GetDistCorrFromLookUpTable
(
 const Float_t x[], 
 Short_t roc,
 Float_t dx[],
 AliTPCLookUpTable3DInterpolatorD * lookUp
 ) 
{
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  TMath::Sqrt( x[0]*x[0] + x[1]*x[1] ) ;
  phi    =  TMath::ATan2(x[1],x[0]) ;
  
  while (phi>TMath::Pi()) phi-=TMath::TwoPi();
  while (phi<-TMath::Pi()) phi+=TMath::TwoPi();

  
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;


	
  // Get distortion
  lookUp->GetValue(r,phi,z,dR,dPhiR,dZ);
  
  //dR      = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
  //			  fLookupRList, fLookupZList, fLookupPhiList, lookUpR );
  //dPhiR    = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
  //			  fLookupRList, fLookupZList, fLookupPhiList, lookUpPhiR);
  //dZ = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
  //			  fLookupRList, fLookupZList, fLookupPhiList, lookUpDz   );

  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + fCorrectionFactor * dPhiR / r;
    r   =  r   + fCorrectionFactor * dR;
  }
  dZ =  fCorrectionFactor * dZ ;

  //printf("correction factor %f\n",fCorrectionFactor);
  // Calculate correction in cartesian coordinates
  dx[0] = (r * TMath::Cos(phi) - x[0]);
  dx[1] = (r * TMath::Sin(phi) - x[1]);
  dx[2] = dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}


Double_t AliTPCSpaceCharge3DDriftLine::Interpolate3DTableCyl
( 
 Int_t order, 
 Double_t r,   
 Double_t z,
 Double_t phi,   	
 Int_t  nr,    
 Int_t  nz,
 Int_t  nphi,    	
 const Double_t rlist[], 
 const Double_t zlist[],
 const Double_t philist[], 
	
 TMatrixD **arrayofArrays 
  ) 
{
  /// Interpolate table (TMatrix format) - 3D interpolation
  /// Float version (in order to decrease the OCDB size)

  static  Int_t ilow = 0, jlow = 0, klow = 0, m=0;
  Float_t saveArray[5]= {0.,0.,0.,0.,0.};
  Float_t savedArray[5]= {0.,0.,0.,0.,0.} ;

  if (phi < 0.0) phi = TMath::TwoPi() + phi;
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

  Search( nr, rlist, r, ilow   ) ;
  Search( nz, zlist, z, jlow   ) ;  
  Search( nphi, philist, phi, klow   ) ;

	 
	

  if ( ilow < 0 ) {
    ilow = 0 ; 
  }
  if ( jlow < 0 ) {
    jlow = 0 ;
  }
  if (klow < 0) klow = nphi + klow;
  
  if ( ilow + order  >=    nr- 1 ) {
    ilow =   nr - 1 - order ;
  }
  if ( jlow + order  >=    nz - 1 ) {
    jlow =   nz - 1 - order ;
  }

  for ( Int_t k = 0 ; k <  order + 1 ; k++ )
    {
      m = (klow + k) % nphi;			
      TMatrixD &table = *arrayofArrays[m] ;
      
      for ( Int_t i = ilow ; i < ilow + order + 1 ; i++ )
	{
				
	  saveArray[i-ilow] = Interpolate( &zlist[jlow], &table(i,jlow), order, z )   ;
				
	}
			
      savedArray[k] = Interpolate( &rlist[ilow], saveArray, order, r )  ;
      
    }
  return( InterpolatePhi( &philist[0], klow, nphi,  savedArray, order, phi ) )   ;
}


Double_t AliTPCSpaceCharge3DDriftLine::Interpolate3DTableCyl
( 
 Int_t order, 
 Double_t r,   
 Double_t z,
 Double_t phi,   	
 Int_t  nr,    
 Int_t  nz,
 Int_t  nphi,    	
 const Double_t rlist[], 
 const Double_t zlist[],
 const Double_t philist[], 
	
 TMatrixD **arrayofArrays,
 const Int_t zlow
  ) 
{
  /// Interpolate table (TMatrix format) - 3D interpolation
  /// Float version (in order to decrease the OCDB size)

  static  Int_t ilow = 0, jlow = 0, klow = 0, m=0;
  Bool_t bExtrapolate = kFALSE;
  Float_t saveArray[5]= {0.,0.,0.,0.,0.};
  Float_t savedArray[5]= {0.,0.,0.,0.,0.} ;

  if (phi < 0.0) phi = TMath::TwoPi() + phi;
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

  Search( nr, rlist, r, ilow   ) ;
  Search( nz, zlist, z, jlow   ) ;  
  Search( nphi, philist, phi, klow   ) ;

  //if (jlow + order < jhigh)
  //printf("z=%f,jlow=%d,zlow=%d\n",z,jlow,zlow);
  jlow = zlow;
	

  if ( ilow < 0 ) {
    ilow = 0 ; 
  }
  if ( jlow < 0 ) {
    jlow = 0 ;
  }
	
  if (klow < 0) klow = nphi + klow;
  
  if ( ilow + order  >=    nr- 1 ) {
    ilow =   nr - 1 - order ;
  }
  if ( jlow + order  >=    nz - 1 ) {
    jlow =   nz - 1 - order ;
  }

  for ( Int_t k = 0 ; k <  order + 1 ; k++ )
    {
      m = (klow + k) % nphi;			
      TMatrixD &table = *arrayofArrays[m] ;
      
      for ( Int_t i = ilow ; i < ilow + order + 1 ; i++ )
	{
					
	  saveArray[i-ilow] = Interpolate( &zlist[jlow], &table(i,jlow), order, z )   ;
				
	}
			
      savedArray[k] = Interpolate( &rlist[ilow], saveArray, order, r )  ;
      
    }
  return( InterpolatePhi( &philist[0], klow, nphi,  savedArray, order, phi ) )   ;
}


Float_t AliTPCSpaceCharge3DDriftLine::Interpolate3DTableCyl
( 
 Int_t order, 
 Double_t r,   
 Double_t z,
 Double_t phi,   
 Int_t  nr, 
 Int_t  nz,
 Int_t  nphi,    
 const Double_t rlist[], 
 const Double_t zlist[],
 const Double_t philist[], 
 TMatrixF **arrayofArrays 
  ) 
{
  /// Interpolate table (TMatrix format) - 3D interpolation
  /// Float version (in order to decrease the OCDB size)
  //printf("%f,%f,%f\n",r,z,phi);
  static  Int_t ilow = 0, jlow = 0, klow = 0 ,m;
  Float_t saveArray[5]= {0.,0.,0.,0.,0.};
  Float_t savedArray[5]= {0.,0.,0.,0.,0.} ;

  if (phi < 0.0) phi = TMath::TwoPi() + phi;
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

  Search( nr, rlist, r, ilow   ) ;
  Search( nz, zlist, z, jlow   ) ;
  Search( nphi, philist, phi, klow   ) ;
	
  if ( ilow < 0 ) {
    ilow = 0 ;   // check if out of range
	
  }
  if ( jlow < 0 ) {
    jlow = 0 ;
  }
  
  if (klow < 0) klow = nphi + klow;
  
  
  if ( ilow + order  >=    nr- 1 ) {
    ilow =   nr - 1 - order ;
  }
  if ( jlow + order  >=    nz - 1 ) {
    jlow =   nz - 1 - order ;
  }
	
	
  for ( Int_t k = 0 ; k < order + 1 ; k++ )
    {
      m = (klow + k) % nphi;
			
				
			    
      TMatrixF &table = *arrayofArrays[m] ;
      
      for ( Int_t i = ilow ; i < ilow + order + 1 ; i++ )
	{
	  saveArray[i-ilow] = Interpolate( &zlist[jlow], &table(i,jlow), order, z)   ;
	}
      savedArray[k] = Interpolate( &rlist[ilow], saveArray, order, r)  ;
    }
  return( InterpolatePhi( &philist[0], klow, nphi, savedArray, order, phi ) )   ;
}




Double_t AliTPCSpaceCharge3DDriftLine::InterpolatePhi
( 
 const Double_t xArray[], 
 const Int_t ilow,
 const Int_t nx,
 const Float_t yArray[],
 Int_t order, 
 Double_t x 
  )
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.

  Int_t i0 = ilow;
  Double_t xi0 = xArray[ilow];
	
  Int_t i1 = (ilow + 1) % nx;
  Double_t xi1 = xArray[i1];
  Int_t i2 = (ilow + 2) % nx;
  Double_t xi2 = xArray[i2];

  if (x < 0)       x += TMath::TwoPi();
  if (xi0 < 0) xi0 += TMath::TwoPi();
  if (xi1 < 0) xi1 += TMath::TwoPi();
  if (xi2 < 0) xi2 += TMath::TwoPi();

  if (xi0 > xi1) 
    xi1 += TMath::TwoPi();
	
  if (xi1 > xi2)
    xi2 += TMath::TwoPi();

  if (xi0 > x)
    x  += TMath::TwoPi();
  printf("(%f,%f,%f,%f)\n",x,xi0,xi1,xi2);

	
  Double_t y ;
 
  
  if ( order == 2 ) {                // Quadratic Interpolation = 2 Ernst: should be >=2?
		
    y  = (x-xi1) * (x-xi2) * yArray[0] / ( (xi0-xi1) * (xi0-xi2) ) ;
    y += (x-xi2) * (x-xi0) * yArray[1] / ( (xi1-xi2) * (xi1-xi0) ) ;
    y += (x-xi0) * (x-xi1) * yArray[2] / ( (xi2-xi0) * (xi2-xi1) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xi0 ) / (xi1 - xi0 ) ;
  }

  return (y);
}





void AliTPCSpaceCharge3DDriftLine::GetDistCorrFromLookUpTableCyl
(
 const Float_t x[], 
 Short_t roc,
 Float_t dx[],
 TMatrixD**lookUpR,
 TMatrixD**lookUpPhiR,
 TMatrixD**lookUpDz
 ) 
{
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  
  
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

	
  // Get distortion
  //printf("(%f,%f,%f)\n",r,z,phi);
  dR      = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
				  fLookupRList, fLookupZList, fLookupPhiList, lookUpR );
  dPhiR    = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
				   fLookupRList, fLookupZList, fLookupPhiList, lookUpPhiR);
  dZ = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
			     fLookupRList, fLookupZList, fLookupPhiList, lookUpDz   );
  //printf("(%f,%f,%f)\n",dR,dPhiR,dZ);
  
  //printf("correction factor %f\n",fCorrectionFactor);
  // Calculate correction in cartesian coordinates
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}






void AliTPCSpaceCharge3DDriftLine::GetDistCorrFromLookUpTableCyl
(
 const Float_t x[], 
 Short_t roc,
 Float_t dx[],
 AliTPCLookUpTable3DInterpolatorD * lookUp
 ) 
{
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  
 
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

  lookUp->GetValue(r,phi,z,dR,dPhiR,dZ);
  // Get distortion
  //printf("(%f,%f,%f)\n",r,z,phi);
  //dR      = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
  //			  fLookupRList, fLookupZList, fLookupPhiList, lookUpR );
  //dPhiR    = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
  //			  fLookupRList, fLookupZList, fLookupPhiList, lookUpPhiR);
  //dZ = Interpolate3DTableCyl(fInterpolationOrder, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
  //				fLookupRList, fLookupZList, fLookupPhiList, lookUpDz   );
  //printf("(%f,%f,%f)\n",dR,dPhiR,dZ);
  
  //printf("correction factor %f\n",fCorrectionFactor);
  // Calculate correction in cartesian coordinates
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}



TH2F* AliTPCSpaceCharge3DDriftLine::CreateHistoDistDRinXY(Float_t z,Int_t nx,Int_t ny) {
  /// Simple plot functionality.
  /// Returns a 2d hisogram which represents the corrections in radial direction (dr)
  /// in respect to position z within the XY plane.
  /// The histogramm has nx times ny entries.

  // AliTPCParam* tpcparam = new AliTPCParamSR; // removed
  o2::TPC::Mapper &mapper = o2::TPC::Mapper::instance(); // added; mapper for radius calculation

  TH2F *h=CreateTH2F("dr_xy", TString::Format("%s: DRinXY Z=%2.0f", GetTitle(),z).Data(),"x [cm]","y [cm]","dr [cm]",
		     nx,-250.,250.,ny,-250.,250.);
  Float_t x[3],dx[3];
  x[2]=z;
  Int_t roc=z>0.?0:18; // FIXME
  for (Int_t iy=1;iy<=ny;++iy) {
    x[1]=h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      x[0]=h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x,roc,dx);
      
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      
      // if (tpcparam->GetPadRowRadii(0,0)<=r0 && r0<=tpcparam->GetPadRowRadii(36,95)) { // removed
      if (mapper.getPadRegionInfo(0).getRadiusFirstRow()<=r0 && r0<=(mapper.getPadRegionInfo(9).getRadiusFirstRow()+mapper.getPadRegionInfo(9).getNumberOfPadRows()*mapper.getPadRegionInfo(9).getPadHeight())) { // added
	Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
	h->SetBinContent(ix,iy,r1-r0);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  // delete tpcparam; // removed
  return h;
}



TH2F* AliTPCSpaceCharge3DDriftLine::CreateHistoDistDRPhiinXY(Float_t z,Int_t nx,Int_t ny) {
  /// Simple plot functionality.
  /// Returns a 2d hisogram which represents the corrections in rphi direction (drphi)
  /// in respect to position z within the XY plane.
  /// The histogramm has nx times ny entries.

  // AliTPCParam* tpcparam = new AliTPCParamSR; // removed
  o2::TPC::Mapper &mapper = o2::TPC::Mapper::instance(); // added; mapper for radius calculation

  TH2F *h=CreateTH2F("drphi_xy",TString::Format("%s: DRPhiinXY Z=%2.0f", GetTitle(),z).Data(),"x [cm]","y [cm]","drphi [cm]",
		     nx,-250.,250.,ny,-250.,250.);
  Float_t x[3],dx[3];
  x[2]=z;
  Int_t roc=z>0.?0:18; // FIXME
  for (Int_t iy=1;iy<=ny;++iy) {
    x[1]=h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      x[0]=h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x,roc,dx);
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      // if (tpcparam->GetPadRowRadii(0,0)<=r0 && r0<=tpcparam->GetPadRowRadii(36,95)) { // removed 
      if (mapper.getPadRegionInfo(0).getRadiusFirstRow()<=r0 && r0<=(mapper.getPadRegionInfo(9).getRadiusFirstRow()+mapper.getPadRegionInfo(9).getNumberOfPadRows()*mapper.getPadRegionInfo(9).getPadHeight())) { // added
 	Float_t phi0=TMath::ATan2(x[1]      ,x[0]      );
	Float_t phi1=TMath::ATan2(x[1]+dx[1],x[0]+dx[0]);

	Float_t dphi=phi1-phi0;
	if (dphi<TMath::Pi()) dphi+=TMath::TwoPi();
	if (dphi>TMath::Pi()) dphi-=TMath::TwoPi();

	h->SetBinContent(ix,iy,r0*dphi);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  // delete tpcparam; // removed
  return h;
}

TH2F* AliTPCSpaceCharge3DDriftLine::CreateHistoDistDZinXY(Float_t z,Int_t nx,Int_t ny) {
  /// Simple plot functionality.
  /// Returns a 2d hisogram which represents the corrections in longitudinal direction (dz)
  /// in respect to position z within the XY plane.
  /// The histogramm has nx times ny entries.

  // AliTPCParam* tpcparam = new AliTPCParamSR; // removed
  o2::TPC::Mapper &mapper = o2::TPC::Mapper::instance(); // added; mapper for radius calculation

  TH2F *h=CreateTH2F("dz_xy",TString::Format("%s: DZinXY Z=%2.0f", GetTitle(),z).Data(),"x [cm]","y [cm]","dz [cm]",
		     nx,-250.,250.,ny,-250.,250.);
  Float_t x[3],dx[3];
  x[2]=z;
  Int_t roc=z>0.?0:18; // FIXME
  for (Int_t iy=1;iy<=ny;++iy) {
    x[1]=h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      x[0]=h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x,roc,dx);
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      // if (tpcparam->GetPadRowRadii(0,0)<=r0 && r0<=tpcparam->GetPadRowRadii(36,95)) { // removed
      if (mapper.getPadRegionInfo(0).getRadiusFirstRow()<=r0 && r0<=(mapper.getPadRegionInfo(9).getRadiusFirstRow()+mapper.getPadRegionInfo(9).getNumberOfPadRows()*mapper.getPadRegionInfo(9).getPadHeight())) { // added
	h->SetBinContent(ix,iy,dx[2]);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  // delete tpcparam; // removed
  return h;
}


/// Use 3D space charge map as an optional input
/// The layout of the input histogram is assumed to be: (phi,r,z)
/// Density histogram is expreseed is expected to bin in  C/m^3
///
/// Standard histogram interpolation is used in order to use the density at center of voxel
/// Warning: Since  value at phi=0, is discontinued 
void	AliTPCSpaceCharge3DDriftLine::SetInputSpaceCharge
(
 TH3 * hisSpaceCharge3D, 
 Double_t norm
 )
{
  fSpaceChargeHistogram3D = hisSpaceCharge3D;
  TMatrixF * scDensity;

  Double_t  r, phi, z ;
  Double_t rmin=hisSpaceCharge3D->GetYaxis()->GetBinCenter(0);
  Double_t rmax=hisSpaceCharge3D->GetYaxis()->GetBinUpEdge(hisSpaceCharge3D->GetYaxis()->GetNbins());
  Double_t zmin=hisSpaceCharge3D->GetZaxis()->GetBinCenter(0);
  Double_t zmax=hisSpaceCharge3D->GetZaxis()->GetBinCenter(hisSpaceCharge3D->GetZaxis()->GetNbins());

  for ( Int_t k = 0 ; k < fNPhiSlices ; k++ ) {
    phi = fLookupPhiList[k] ;
    scDensity   =  fSCdensityDistribution[k]  ;
    for ( Int_t j = 0 ; j < fNZColumns ; j++ ) {
      z = fLookupZList[j] ;
      for ( Int_t i = 0 ; i < fNRRows ; i++ ) {
	// Full 3D configuration ...
	r = fLookupRList[i] ;
	if (r>rmin && r<rmax && z>zmin && z< zmax){
	  // scDensity(i,j) = norm* InterpolatePhi(fSpaceChargeHistogram3D,phi,r,z);
	  (*scDensity)(i,j) = norm* fSpaceChargeHistogram3D->Interpolate(phi,r,z);
	}
      }
    }
  }
  fInitLookUp = kFALSE;
}


/// for correctness analysis
void AliTPCSpaceCharge3DDriftLine::LocalDistCorrDzExact
( 
 TFormula *intDrDzF,
 TFormula *intDphiDzF,
 TFormula *intDzDzF,
 Double_t *rlist,
 Double_t *philist,
 Double_t *zlist,
 TMatrixD** matricesDistDrDz, 
 TMatrixD** matricesDistDphiRDz, 
 TMatrixD** matricesDistDz, 		
 TMatrixD** matricesCorrDrDz, 
 TMatrixD** matricesCorrDphiRDz, 
 TMatrixD** matricesCorrDz, 
 const Int_t rrow, 
 const Int_t zcolumn, 
 const Int_t phiSlice,  
 const Double_t ezField,
 const Int_t side
  )
{
  Double_t r0,z0,phi0,z1;
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEphiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;

  // pointer declaration
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
		
		
  if (side == 1) {
    for (Int_t j =0;j<zcolumn;j++) zlist[j] = -1 * zlist[j];
  }	
  //printf("c0=%E,c1=%E,ezField=%E,dvdE=%E\n",c0,c1,ezField,dvdE);
	
  // Initialization for j == colomn-1 integration is 0.0
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
    for (Int_t i=0; i< rrow;i++) {
      (*distDrDz)(i,zcolumn-1) =  0.0 ;			
      (*distDphiRDz)(i,zcolumn-1) =  0.0 ;			
      (*distDz)(i,zcolumn-1) =  0.0 ;			
			
      (*corrDrDz)(i,0) =  0.0 ;			
      (*corrDphiRDz)(i,0) =  0.0 ;			
      (*corrDz)(i,0) =  0.0 ;			
    }
		
			
  }
	
  // for this case
  // use trapezoidal rule assume no ROC displacement
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    phi0 = philist[m];
		
    distDrDz  =  matricesDistDrDz[m] ;
    distDphiRDz  =  matricesDistDphiRDz[m] ;
    distDz  =  matricesDistDz[m] ;
		
    corrDrDz  =  matricesCorrDrDz[m] ;
    corrDphiRDz  =  matricesCorrDphiRDz[m] ;
    corrDz  =  matricesCorrDz[m] ;
		
    for (Int_t j=0; j< zcolumn-1;j++) {
      z0 = zlist[j];
      z1 = zlist[j+1];
			
				
      for (Int_t i=0; i< rrow;i++)  {
	r0 = rlist[i];
	localIntErOverEz = (intDrDzF->Eval(r0,phi0,z1) - intDrDzF->Eval(r0,phi0,z0)) /(-1*ezField) ;
	localIntEphiOverEz = (intDphiDzF->Eval(r0,phi0,z1) - intDphiDzF->Eval(r0,phi0,z0)) /(-1*ezField) ;				
	localIntDeltaEz = intDzDzF->Eval(r0,phi0,z1) - intDzDzF->Eval(r0,phi0,z0);				
				
	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
	(*distDz)(i,j)      = localIntDeltaEz  * o2::TPC::AliTPCPoissonSolver::fgkdvdE * o2::TPC::AliTPCPoissonSolver::fgkdvdE ; // two times?				
	(*corrDrDz)(i,j+1) 		= -1* (*distDrDz)(i,j) ;				
	(*corrDphiRDz)(i,j+1) = -1* (*distDphiRDz)(i,j);
	(*corrDz)(i,j+1)      = -1* (*distDz)(i,j);							
      }			
    }				
  }		
}


/*
  TTree* AliTPCSpaceCharge3DDriftLine::CreateDistortionTree
  (
  Double_t step
  )
  {
  /// create the distortion tree on a mesh with granularity given by step
  /// return the tree with distortions at given position
  /// Map is created on the mesh with given step size
  /// type - 0: Call GetDistortion()
  ///        1: Call GetDistortionIntegralDz()


  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("correction%s.root",GetName()));
  Float_t xyz[3];     // current point
  Float_t dist[3];    // distorion
  Float_t corr[3];    // correction
  Float_t xyzdist[3]; // distorted point
  Float_t xyzcorr[3]; // corrected point

  //AliMagF* mag= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  //if (!mag) AliError("Magnetic field - not initialized");

  Int_t roc;
  // AliTPCParam* tpcparam = new AliTPCParamSR; // removed
  o2::TPC::Mapper &mapper = o2::TPC::Mapper::instance(); // added; mapper for radius calculation
  Double_t r,phi,rdist,phidist,drdist,drphidist,rcorr,phicorr,drcorr,drphicorr;
  for (Double_t x= -250; x<250; x+=step){
  for (Double_t y= -250; y<250; y+=step){
			
  r    = TMath::Sqrt(x*x+y*y);
      
  // if (tpcparam->GetPadRowRadii(0,0)>r || r>tpcparam->GetPadRowRadii(36,95))    continue; // removed
  if (mapper.getPadRegionInfo(0).getRadiusFirstRow()>r || r>(mapper.getPadRegionInfo(9).getRadiusFirstRow()+mapper.getPadRegionInfo(9).getNumberOfPadRows()*mapper.getPadRegionInfo(9).getPadHeight())) continue;  // added
      
  //printf("(%f,%f)\n",x,y);
  phi  = TMath::ATan2(y,x);
     
			
  for (Double_t z= -250; z<250; z+=step){
  roc=(z>0)?0:18;
  xyz[0]=x;
  xyz[1]=y;
  xyz[2]=z;

  GetDistortion(xyz, roc, dist);
        
  for (Int_t i=0; i<3; ++i) {
  xyzdist[i]=xyz[i]+dist[i];         
  }


  GetCorrection(xyzdist, roc, corr);

  for (Int_t i=0; i<3; ++i) {
  xyzcorr[i]=xyzdist[i]+corr[i];
  }
        
  // === r, rphi + residuals for the distorted point =========================
  rdist    = TMath::Sqrt(xyzdist[0]*xyzdist[0]+xyzdist[1]*xyzdist[1]);
  phidist  = TMath::ATan2(xyzdist[1],xyzdist[0]);
  //rdist = xyzdist[0];
  //phidist = xyzdist[1];

  while ((phidist-phi)>TMath::Pi()) phidist-=TMath::TwoPi();
  while ((phidist-phi)<-TMath::Pi()) phidist+=TMath::TwoPi();

  drdist=rdist-r;
  drphidist=(phidist-phi)*r;

  // === r, rphi + residuals for the corrected point =========================
  rcorr    = TMath::Sqrt(xyzcorr[0]*xyzcorr[0]+xyzcorr[1]*xyzcorr[1]);
  phicorr  = TMath::ATan2(xyzcorr[1],xyzcorr[0]);
  //rcorr = xyzcorr[0];
  //phicorr = xyzcorr[1];
        
  while ((phicorr-phidist)>TMath::Pi()) phicorr-=TMath::TwoPi();
  while ((phicorr-phidist)<-TMath::Pi()) phicorr+=TMath::TwoPi();

        
  drcorr=rcorr-rdist;
  drphicorr=(phicorr-phidist)*r;

  // === get b field ===============
  // Double_t bxyz[3]={0.,0.,0.};
  // Double_t dblxyz[3] = {Double_t(xyzdist[0]),Double_t(xyzdist[1]),Double_t(xyzdist[2])};
  // Double_t br	= 0.;
  // Double_t brfi = 0.;
  // if (mag) {
  //   mag->Field(dblxyz,bxyz);
  //   if(rdist>0){
  //     br = (bxyz[0]*xyz[0]+bxyz[1]*xyz[1])/rdist;
  //     brfi = (-bxyz[0]*xyz[1]+bxyz[1]*xyz[0])/rdist;
  //   }
  // }
  (*pcstream)<<"distortion"<<
  "x="  << x   <<           // original position
  "y="  << y   <<
  "z="  << z   <<
  "r="  << r   <<
  "phi="<< phi <<
  //
  "x_dist="    << xyzdist[0] <<      // distorted position
  "y_dist="    << xyzdist[1] <<
  "z_dist="    << xyzdist[2] <<
  "r_dist="    << rdist      <<
  "phi_dist="  << phidist    <<
  //
  "dx_dist="   << dist[0]    <<     // distortion
  "dy_dist="   << dist[1]    <<
  "dz_dist="   << dist[2]    <<
  "dr_dist="   << drdist     <<
  "drphi_dist="<< drphidist  <<
  //
  "x_corr="    << xyzcorr[0] <<      // corrected position
  "y_corr="    << xyzcorr[1] <<
  "z_corr="    << xyzcorr[2] <<
  "r_corr="    << rcorr      <<
  "phi_corr="  << phicorr    <<
  //
  "dx_corr="   << corr[0]    <<     // correction
  "dy_corr="   << corr[1]    <<
  "dz_corr="   << corr[2]    <<
  "dr_corr="   << drcorr     <<
  "drphi_corr="<< drphicorr  <<
  "\n";
  }
  }
  }
  delete pcstream;
  TFile f(Form("correction%s.root",GetName()));
  TTree * tree = (TTree*)f.Get("distortion");  
  
  return tree;
  }
*/ // removed due to TTreeSRedirector


Double_t AliTPCSpaceCharge3DDriftLine::InterpolatePhi
( 
 TH3 * h3,
 const Double_t phi,
 const Double_t r,
 const Double_t z
  )
{
	
  Int_t ubx =  h3->GetXaxis()->FindBin(phi);
  if ( phi < h3->GetXaxis()->GetBinCenter(ubx) ) ubx -= 1;
  Int_t obx = ubx + 1;
	
  Int_t uby =  h3->GetYaxis()->FindBin(r);
  if (r < h3->GetYaxis()->GetBinCenter(uby) ) uby -= 1;
  Int_t oby = uby + 1;
	
  Int_t ubz =  h3->GetZaxis()->FindBin(z);
  if ( z < h3->GetZaxis()->GetBinCenter(ubz) ) ubz -= 1;
  Int_t obz = ubz + 1;
	
  if ( uby <=0 || ubz <= 0 ||
       oby > h3->GetYaxis()->GetNbins() || obz > h3->GetZaxis()->GetNbins() ) 
    {
      // 	printf("Phi --> Interpolate Cannot interpolate outside histogram domain. (%f,%f,%f)\n",phi,r,z);
      return 0;
    }
	
  if (ubx <=0) 	ubx = h3->GetXaxis()->GetNbins();
	
  if (obx > h3->GetXaxis()->GetNbins()) obx = 1;
		
  Double_t xw = h3->GetXaxis()->GetBinCenter(obx) - h3->GetXaxis()->GetBinCenter(ubx);
  Double_t yw = h3->GetYaxis()->GetBinCenter(oby) - h3->GetYaxis()->GetBinCenter(uby);
  Double_t zw = h3->GetZaxis()->GetBinCenter(obz) - h3->GetZaxis()->GetBinCenter(ubz);
 
  Double_t xd = (phi - h3->GetXaxis()->GetBinCenter(ubx)) / xw;
  Double_t yd = (r - h3->GetYaxis()->GetBinCenter(uby)) / yw;
  Double_t zd = (z - h3->GetZaxis()->GetBinCenter(ubz)) / zw;


  Double_t v[] = { h3->GetBinContent( ubx, uby, ubz ), h3->GetBinContent( ubx, uby, obz ),
		   h3->GetBinContent( ubx, oby, ubz ), h3->GetBinContent( ubx, oby, obz ),
		   h3->GetBinContent( obx, uby, ubz ), h3->GetBinContent( obx, uby, obz ),
		   h3->GetBinContent( obx, oby, ubz ), h3->GetBinContent( obx, oby, obz ) };
	
  Double_t i1 = v[0] * (1 - zd) + v[1] * zd;
  Double_t i2 = v[2] * (1 - zd) + v[3] * zd;
  Double_t j1 = v[4] * (1 - zd) + v[5] * zd;
  Double_t j2 = v[6] * (1 - zd) + v[7] * zd;
  Double_t w1 = i1 * (1 - yd) + i2 * yd;
  Double_t w2 = j1 * (1 - yd) + j2 * yd;
  Double_t result = w1 * (1 - xd) + w2 * xd;
  return result;
	
}





TH2F * AliTPCSpaceCharge3DDriftLine::CreateHistoSCinXY
(
 Float_t z, 
 Int_t nx, 
 Int_t ny
 ) 
{
  /// return a simple histogramm containing the space charge distribution (input for the calculation)

  TH2F *h=CreateTH2F("spaceCharge",GetTitle(),"x [cm]","y [cm]","#rho_{sc} [C/m^{3}/e_{0}]",
		     nx,-250.,250.,ny,-250.,250.);

  for (Int_t iy=1;iy<=ny;++iy) {
    Double_t yp = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      Double_t xp = h->GetXaxis()->GetBinCenter(ix);

      Float_t r = TMath::Sqrt(xp*xp+yp*yp);
      Float_t phi = TMath::ATan2(yp,xp);

      if (85.<=r && r<=250.) {
	Float_t sc = GetSpaceChargeDensity(r,phi,z)/o2::TPC::AliTPCPoissonSolver::fgke0; // in [C/m^3/e0]
	h->SetBinContent(ix,iy,sc);
      } else {
	h->SetBinContent(ix,iy,0.);
      }
    }
  }

  return h;
}

TH2F * AliTPCSpaceCharge3DDriftLine::CreateHistoSCinZR
(
 Float_t phi, 
 Int_t nz, 
 Int_t nr
 ) 
{
  /// return a simple histogramm containing the space charge distribution (input for the calculation)

  TH2F *h=CreateTH2F("spaceCharge",GetTitle(),"z [cm]","r [cm]","#rho_{sc} [C/m^{3}/e_{0}]",
		     nz,-250.,250.,nr,85.,250.);

  for (Int_t ir=1;ir<=nr;++ir) {
    Float_t r = h->GetYaxis()->GetBinCenter(ir);
    for (Int_t iz=1;iz<=nz;++iz) {
      Float_t z = h->GetXaxis()->GetBinCenter(iz);
      Float_t sc = GetSpaceChargeDensity(r,phi,z)/o2::TPC::AliTPCPoissonSolver::fgke0; // in [C/m^3/e0]
      h->SetBinContent(iz,ir,sc);
    }
  }

  return h;
}



Float_t  AliTPCSpaceCharge3DDriftLine::GetSpaceChargeDensity(Float_t r, Float_t phi, Float_t z) {
  /// returns the (input) space charge density at a given point according
  /// Note: input in [cm], output in [C/m^3/e0] !!

  while (phi<0) phi += TMath::TwoPi();
  while (phi>TMath::TwoPi()) phi -= TMath::TwoPi();


  // Float_t sc =fSCdensityDistribution->Interpolate(r0,phi0,z0);
  const Int_t order = 1; //

  const Float_t  sc = Interpolate3DTableCyl(order, r, z, phi, fNRRows, fNZColumns, fNPhiSlices,
					    fLookupRList, fLookupZList, fLookupPhiList, fSCdensityDistribution );

  return sc;
}


// follow the drift for exact function
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDz 
(
 TFormula * intDrDzF,
 TFormula * intDphiDzF,
 TFormula * intDzDzF,
 const Double_t ezField,
 TMatrixD** matricesGDistDrDz,  
 TMatrixD** matricesGDistDphiRDz, 
 TMatrixD** matricesGDistDz,
 TMatrixD** matricesGCorrDrDz,  
 TMatrixD** matricesGCorrDphiRDz, 
 TMatrixD** matricesGCorrDz, 
 

 TMatrixD** matricesGCorrIrregularDrDz,  
 TMatrixD** matricesGCorrIrregularDphiRDz, 
 TMatrixD** matricesGCorrIrregularDz, 
 TMatrixD** matricesRIrregular,  
 TMatrixD** matricesPhiIrregular, 
 TMatrixD** matricesZIrregular, 


 const Int_t rrow,  
 const Int_t zcolumn, 
 const Int_t phiSlice,
 const Double_t *rlist,
 const Double_t *philist,
 const Double_t *zlist
 )
{

	
	
  Float_t dr,dphir,dz,ddr,ddrphi,ddz;
  Float_t radius0, phi0, z0, radius,phi,z,radiusc,z1,localIntErOverEz,localIntEphiOverEz,localIntDeltaEz,z2,gridzsize;
  radiusc = 0.0;
  radius = 0.0;
  TMatrixD * mDistDrDz;
  TMatrixD * mDistDphiRDz;
  TMatrixD * mDistDz;
	
  TMatrixD * mCorrDrDz;
  TMatrixD * mCorrDphiRDz;
  TMatrixD * mCorrDz;
	


  TMatrixD * mCorrIrregularDrDz;
  TMatrixD * mCorrIrregularDphiRDz;
  TMatrixD * mCorrIrregularDz;

  TMatrixD * mRIrregular;
  TMatrixD * mPhiIrregular;
  TMatrixD * mZIrregular;
	
	
  gridzsize = zlist[zcolumn - 1] - zlist[zcolumn - 2] ;
	

  // initialized for zcolumn - 1
  Int_t j = zcolumn  - 1;
  z0 = zlist[j] ;		
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
    phi0 = philist[m];
			
    mDistDrDz = matricesGDistDrDz[m];
    mDistDphiRDz = matricesGDistDphiRDz[m];
    mDistDz = matricesGDistDz[m];
	
    // 
    mCorrDrDz = matricesGCorrDrDz[m];
    mCorrDphiRDz = matricesGCorrDphiRDz[m];
    mCorrDz = matricesGCorrDz[m];


    mCorrIrregularDrDz  = matricesGCorrIrregularDrDz[m];
    mCorrIrregularDphiRDz = matricesGCorrIrregularDphiRDz[m];
    mCorrIrregularDz = matricesGCorrIrregularDz[m];

    mRIrregular = matricesRIrregular[m];
    mPhiIrregular = matricesPhiIrregular[m];
    mZIrregular= matricesZIrregular[m];

		
    for (Int_t i=0; i< rrow;i++)  {
      // do from j to 0
      // follow the drift
      radius0 = rlist[i];
      phi = phi0;
      radius = radius0;
			
      dr = 0.0;
      dphir = 0.0;
      dz = 0.0;
      ddrphi = 0.0;
			
				
				
				
      (*mDistDrDz)(i,j) = dr;
      (*mDistDphiRDz)(i,j) = dphir;
      (*mDistDz)(i,j) = dz;




      //////////////// use irregular grid look up table for correction
      // set 
      (*mCorrIrregularDrDz)(i,j) = -dr;
      (*mCorrIrregularDphiRDz)(i,j) = -dphir;
      (*mCorrIrregularDz)(i,j) = -dz;


      // distorted point
      (*mRIrregular)(i,j) = radius0 + dr;
      (*mPhiIrregular)(i,j) =	phi0 + (dphir/radius0);
      (*mZIrregular)(i,j) = z0 + dz;
				


      ///////////////				

		
				
    }
  }

  // from end cap to central node
  for (j=zcolumn-2; j>=0;j--) {
    printf("global dist j:%d\n",j);
    z0 = zlist[j] ;		
			
	
    for ( Int_t m = 0 ; m < phiSlice ; m++ ) {
      phi0 = philist[m];
				
      mDistDrDz = matricesGDistDrDz[m];
      mDistDphiRDz = matricesGDistDphiRDz[m];
      mDistDz = matricesGDistDz[m];
      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDphiRDz = matricesGCorrDphiRDz[m];
      mCorrDz = matricesGCorrDz[m];;




      mCorrIrregularDrDz  = matricesGCorrIrregularDrDz[m];
      mCorrIrregularDphiRDz = matricesGCorrIrregularDphiRDz[m];
      mCorrIrregularDz = matricesGCorrIrregularDz[m];

      mRIrregular = matricesRIrregular[m];
      mPhiIrregular = matricesPhiIrregular[m];
      mZIrregular= matricesZIrregular[m];

			
      for (Int_t i=0; i< rrow;i++)  {
	// do from j to 0
	// follow the drift
	radius0 = rlist[i];
	phi = phi0;
	radius = radius0;
				
	dr = 0.0;
	dphir = 0.0;
	dz = 0.0;
	ddrphi = 0.0;
	for (Int_t jj = j; jj < zcolumn;jj++) {
	  // interpolation the local distortion for current position
	  phi += ddrphi/radius;
	  radius = radius0 + dr;
	  z = zlist[jj] + dz;
	  z1 = z + gridzsize;
					
					
	  if (phi < 0.0) phi = TMath::TwoPi() + phi;
	  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

	  //				lookupLocalDist->GetValue(radius,phi,z,ddr,ddrphi,ddz);
	  // get from exact function
					
	  localIntErOverEz = (intDrDzF->Eval(radius,phi,z1) - intDrDzF->Eval(radius,phi,z)) /(-1*ezField) ;
	  localIntEphiOverEz = (intDphiDzF->Eval(radius,phi,z1) - intDphiDzF->Eval(radius,phi,z)) /(-1*ezField) ;				
	  localIntDeltaEz = intDzDzF->Eval(radius,phi,z1) - intDzDzF->Eval(radius,phi,z);				
				
	  ddr    = fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	  ddrphi = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
	  ddz    = localIntDeltaEz  * o2::TPC::AliTPCPoissonSolver::fgkdvdE * o2::TPC::AliTPCPoissonSolver::fgkdvdE ; // two times?				
				
	  dr += ddr;
	  dphir += ddrphi;
	  dz += ddz;
				
	  z0 = z;
	}
				
	(*mDistDrDz)(i,j) = dr;
	(*mDistDphiRDz)(i,j) = dphir;
	(*mDistDz)(i,j) = dz;
			
	/////////// Irregular

	(*mCorrIrregularDrDz)(i,j) = -dr;
	(*mCorrIrregularDphiRDz)(i,j) = -dphir;
	(*mCorrIrregularDz)(i,j) = -dz;


	// distorted point
	(*mRIrregular)(i,j) = radius0 + dr;
	(*mPhiIrregular)(i,j) =	phi0 + (dphir/radius0);
	(*mZIrregular)(i,j) = z0 + dz;





	////////////
			
	if (j == zcolumn - 2) {
	  radiusc = radius0;
	  //z0 = 
	}
			
	if (j < zcolumn - 1) {	
	  dr = (*mCorrDrDz)(i,j+1);
	  dphir = (*mCorrDphiRDz)(i,j+1);
	  dz = (*mCorrDz)(i,j+1) ;
				
				
	  phi = phi0 + dphir/radiusc;
	  radiusc = radius0 + dr;				
	  z = zlist[j + 1] + dz;
				
	  z1 = z - gridzsize;
					
				
	  //lookupLocalCorr->GetValue(radiusc,phi,z,ddr,ddrphi,ddz);	
				
	  localIntErOverEz = (intDrDzF->Eval(radius,phi,z1) - intDrDzF->Eval(radius,phi,z)) /(-1*ezField) ;
	  localIntEphiOverEz = (intDphiDzF->Eval(radius,phi,z1) - intDphiDzF->Eval(radius,phi,z)) /(-1*ezField) ;				
	  localIntDeltaEz = intDzDzF->Eval(radius,phi,z1) - intDzDzF->Eval(radius,phi,z);				
				
	  ddr    = fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
	  ddrphi = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
	  ddz    = localIntDeltaEz  * o2::TPC::AliTPCPoissonSolver::fgkdvdE * o2::TPC::AliTPCPoissonSolver::fgkdvdE ; // two times?				
				
						
	  dr += ddr;
	  dz += ddz;
	  dphir += ddrphi;
					
				
	  (*mCorrDrDz)(i,j) = dr;
	  (*mCorrDphiRDz)(i,j) = dphir;
	  (*mCorrDz)(i,j) = dz;
				
	}
      }
    }
  }
}




/// Inverse from Global to Local Distortion
/// 
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param rrow Int_t number of grid in r direction
/// \param zcolumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortion
(
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice
 )
{
  Double_t z,phi,r,zaft,zprev,zstep,ddr,ddphir,ddz,zl,dr,dphir,dz,ddphi,dphi,deltaz,zm1,zm2,zp1,zp2,gradient_r,gradient_phir,gradient_z;
  Float_t  x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];

  Int_t roc;




  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;


  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    distDrDz  =  matricesDistDrDz[k] ;
    distDphiRDz  =  matricesDistDphiRDz[k] ;
    distDz  =  matricesDistDz[k] ;

    for ( Int_t i = 0 ; i < rrow; i++ ) {
      (*distDrDz)(i,zcolumn-1)    = 0.0;
      (*distDphiRDz)(i,zcolumn-1) = 0.0;
      (*distDz)(i,zcolumn-1)      = 0.0;
    }
  }


  deltaz =  (zList[1] - zList[0]);
  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) {
    zprev = zList[j+1];
    z = zList[j];
    roc = 0; // FIXME

    // calculate neighborhood at z-direction
    if (j > 0)
      zm1 = zList[j - 1];
    if (j > 1)
      zm2 = zList[j - 2];
    if (j < zcolumn -1)
      zp1 = zList[j + 1];
    if (j < zcolumn - 2)
      zp2 = zList[j + 2];


    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
      distDrDz  =  matricesDistDrDz[k] ;
      distDphiRDz  =  matricesDistDphiRDz[k] ;
      distDz  =  matricesDistDz[k] ;
      phi = phiList[k];

      for ( Int_t i = 0 ; i < rrow; i++ ) {
	// get global distortion

	r = rList[i];

	/// CALCULATE  derivatif for inner z

	if (j > 1  && j < zcolumn - 2)
	  {
	    // get coordinate
	    x[0] = r;
	    x[1] = phi;
	    x[2] = zm2;
	    GetDistortionCylAC(x,roc,dxm2);
	    x[2] = zm1;
	    GetDistortionCylAC(x,roc,dxm1);
	    x[2] = zp1;
	    GetDistortionCylAC(x,roc,dxp1);
	    x[2] = zp2;
	    GetDistortionCylAC(x,roc,dxp2);

	    gradient_r   	= (dxm2[0] - 8.0 *  dxm1[0] +  8.0 *  dxp1[0] + dxp2[0]) /  (12.0 * deltaz);
	    gradient_phir   = (dxm2[1] - 8.0 *  dxm1[1] +  8.0 *  dxp1[1] + dxp2[1]) /  (12.0 * deltaz);
	    gradient_z  	= (dxm2[2] - 8.0 *  dxm1[2] +  8.0 *  dxp1[2] + dxp2[2]) /  (12.0 * deltaz);

	    dr = gradient_r * deltaz;
	    dphir = gradient_phir * deltaz;
	    dz = gradient_z * deltaz;
					
	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = z;
				
	    GetDistortionCylAC(x,roc,dx);
					
					
	    x[0] = r  + dr;
	    x[1] = phi + dphir/r;
	    x[2] = zp1 + dz;
					
	    GetDistortionCylAC(x,roc,pdx);

	    //(*distDrDz)(i,j)    = (pdx[0] - dx[0]);
	    //(*distDz)(i,j)      = (pdx[2] - dx[2]);				
	    //(*distDphiRDz)(i,j) = (pdx[1] - dx[1]);
	    (*distDrDz)(i,j)    = dr;
	    (*distDz)(i,j)      = dphir/r;
	    (*distDphiRDz)(i,j) = dz;
					
					
	  }


	else {




	  // in case j > 0
	  // calculate differential at the point of j
	  if (j > 0) {
	    zaft = zList[j-1];	
	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zaft;
				
	    GetDistortionCylAC(x,roc,dx);
					
				
	    x[2] = zprev;
				
				
	    GetDistortionCylAC(x,roc,pdx);
					
				
				
				
	    dr   = (pdx[0] - dx[0]) /  (2.0 * deltaz);
	    dz   = (pdx[2] - dx[2]) / (2.0 * deltaz);				
	    dphi = (pdx[1] - dx[1]) / (2.0 * deltaz * r);
					
	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = z;
				
	    GetDistortionCylAC(x,roc,dx);
					
					
	    x[0] = r  + (dr * deltaz);
	    x[1] = phi + (dphi * deltaz);
	    x[2] = zprev + (dz * deltaz);
				
	    GetDistortionCylAC(x,roc,pdx);
					
					
	    (*distDrDz)(i,j)    = (pdx[0] - dx[0]);
	    (*distDz)(i,j)      = (pdx[2] - dx[2]);				
	    (*distDphiRDz)(i,j) = (pdx[1] - dx[1]);
				
	  } else {
	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = z;
					
	    GetDistortionCylAC(x,roc,dx);
					
				
	    x[2] = zprev ;
				
				
	    GetDistortionCylAC(x,roc,pdx);
					
				
	    dr   = (pdx[0] - dx[0]) ;
	    dz     = (pdx[2] - dx[2]) ;				
	    dphir = (pdx[1] - dx[1]);
					
					
	    x[0] = r  + dr;
	    x[1] = phi + dphir/r;
	    x[2] = zprev + dz;
				
	    GetDistortionCylAC(x,roc,pdx);
					
					
	    (*distDrDz)(i,j)    = (pdx[0] - dx[0]);
	    (*distDz)(i,j)      = (pdx[2] - dx[2]);				
	    (*distDphiRDz)(i,j) = (pdx[1] - dx[1]);
				
	  }
	}
      }
    }
  }

}



/// Inverse from Global to Local Distortion
/// 
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param rrow Int_t number of grid in r direction
/// \param zcolumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nstep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortion
(
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice,
 const Int_t nstep,
 const Bool_t useCylAC,
 Int_t stepR,
 Int_t stepZ,
 Int_t stepPhi,
 Int_t type
 )
{
  Double_t z,phi,r,zaft,zprev,zstep,ddr,ddphir,ddz,zl,dr,dphir,dz,ddphi,dphi,deltaz,zm1,zm2,zp1,zp2,gradient_r,gradient_phir,gradient_z,r0,z0,phi0;
  Float_t  x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];

	
  Int_t roc;
	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
	

	
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;

  // correction build up for inverse flow
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
	
  TMatrixD * listR;
  TMatrixD * listPhi;
  TMatrixD * listZ;


  /// make interpolator for inverse flow
  ///
  TMatrixD *Mat_corrDrDz[phiSlice];
  TMatrixD *Mat_corrDphiRDz[phiSlice];
  TMatrixD *Mat_corrDz[phiSlice];
	
  TMatrixD *Mat_RList[phiSlice];
  TMatrixD *Mat_PhiList[phiSlice];
  TMatrixD *Mat_ZList[phiSlice];
	
  for (Int_t m=0;m<phiSlice;m++) {
    Mat_corrDrDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDphiRDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDz[m] = new TMatrixD(rrow,zcolumn);

    Mat_RList[m] = new TMatrixD(rrow,zcolumn);
    Mat_PhiList[m] = new TMatrixD(rrow,zcolumn);
    Mat_ZList[m] = new TMatrixD(rrow,zcolumn);
  }
	
  AliTPCLookUpTable3DInterpolatorDFull * lookupInverseCorr =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     rrow,
     Mat_corrDrDz,
     Mat_RList,
     rList,
     phiSlice,
     Mat_corrDphiRDz,
     Mat_PhiList,
     phiList,
     zcolumn,
     Mat_corrDz,
     Mat_ZList,
     zList,
     2,
     stepR,
     stepZ,
     stepPhi,
     type
     );


  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    distDrDz  =  matricesDistDrDz[k] ;
    distDphiRDz  =  matricesDistDphiRDz[k] ;
    distDz  =  matricesDistDz[k] ;
		
    listR = Mat_RList[k];
    listPhi = Mat_PhiList[k];
    listZ  = Mat_ZList[k];

		

    for ( Int_t i = 0 ; i < rrow; i++ ) {
      (*distDrDz)(i,zcolumn-1)    = 0.0;
      (*distDphiRDz)(i,zcolumn-1) = 0.0;
      (*distDz)(i,zcolumn-1)      = 0.0;		
			
      for ( Int_t j = 0 ; j < zcolumn; j++ ) {
	(*listR)(i,j)    = rList[i];
	(*listPhi)(i,j) = phiList[k];
	(*listZ)(i,j)      = zList[j];		
				
      }
				
    }
  }

	


  deltaz =  (zList[1] - zList[0]);
  Int_t ianchor,kanchor, zanchor;

  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) { 

    printf("inversion global on z column  = %d\n",j);
		 
    roc = 0; // FIXME
    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
				
      distDrDz  =  matricesDistDrDz[k] ;
      distDphiRDz  =  matricesDistDphiRDz[k] ;
      distDz  =  matricesDistDz[k] ;

			
      corrDrDz  =  Mat_corrDrDz[k] ;
      corrDphiRDz  =  Mat_corrDphiRDz[k] ;
      corrDz  =  Mat_corrDz[k] ;

			
      listR  =  Mat_RList[k] ;
      listPhi  =  Mat_PhiList[k] ;
      listZ  =  Mat_ZList[k] ;
		

      for ( Int_t i = 0 ; i < rrow; i++ ) {
	// get global distortion
				
	r = rList[i];
	phi = phiList[k];			
	z = zList[j];
	zp1 = zList[j+1];
				
	dr    = 0.0;
	dz    = 0.0;				
	dphir	= 0.0;
				
	if (j < zcolumn - 2) 
	  {
	    // get global distortion of this point

	    //printf("inversion global on z column  = %d\n",j);
	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

				  
				  
	    // get position on 0 based on global distortion
	    r0    = r + dx[0];
	    z0    = zList[zcolumn - 1] + dx[2];				
	    phi0  = phi + (dx[1]/r);
				  
	    // follow electron path from correction table
	    for (Int_t jj=zcolumn-1;jj > j+1; jj--) 
	      {
					
		//if (r0 < o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)   
		//  r0 = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius;
		//if (r0 > o2::TPC::AliTPCPoissonSolver::fgkOFCRadius) 
		//    printf("(%d,%d,%d,%f)\n",i,j,k,r0);
		//  r0 = o2::TPC::AliTPCPoissonSolver::fgkOFCRadius;

		if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
		if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

		ianchor = TMath::FloorNint((r0 - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR);
		kanchor = TMath::FloorNint(phi0/gridSizePhi); 
		zanchor = TMath::FloorNint(z0/gridSizeZ);

		if (j > zcolumn - 5)  
		  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,1,j+2);
		else
		  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,3,j+2);

		phi0 = phi0 + ((dphir)/r0);
		r0 = r0 + ( dr );

		z0 = z0 - deltaz + (dz);

	      }
				  
	    x[0] = r0 ;
	    x[1] = phi0;				
	    x[2] = z0;

	    if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);
					
	  } else if (j == (zcolumn -2))
	  {

	    // in case of zcolumn-2, the global distortion is the local distortion
	    //printf("inversion global on z column  = %d\n",j);

	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

	    x[2] = zList[j+1];


	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);
				  

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);
				  
	  } 

	(*distDrDz)(i,j)    = dr;
	(*distDz)(i,j)      = dz;				
	(*distDphiRDz)(i,j) = dphir;

	r = rList[i];
	phi = phiList[k];			
	z = zList[j];

	(*corrDrDz)(i,j+1)    = -dr;
	(*corrDz)(i,j+1)      = -dz;				
	(*corrDphiRDz)(i,j+1) = -dphir; 

	(*listR)(i,j+1)  =  r + dr;
	(*listPhi)(i,j+1)  =  phi + dphir/r;
	(*listZ)(i,j+1) = zp1 + dz;

				
      }

					
    }

    // copy Vals to correction look up table
    // in case RBF should compute weighted to the interpolant
    lookupInverseCorr->CopyVals(j+1);

		

  }
	
  for (Int_t m=0;m<phiSlice;m++) {
    delete Mat_corrDrDz[m];
    delete Mat_corrDphiRDz[m];
    delete Mat_corrDz[m];
    delete Mat_RList[m];
    delete Mat_PhiList[m];
    delete Mat_ZList[m];
  }
  delete lookupInverseCorr;
       
	
  // delete[] Mat_corrDrDz;
  // delete[] Mat_corrDphiRDz;
  // delete[] Mat_corrDz;
       
  // delete[] Mat_RList;
  // delete[] Mat_PhiList;
  // delete[] Mat_ZList;	
}




/// 
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param rrow Int_t number of grid in r direction
/// \param zcolumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nstep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionTwoStages
(
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice,
 const Int_t nstep,
 const Bool_t useCylAC,
 Int_t stepR,
 Int_t stepZ,
 Int_t stepPhi,
 Int_t type // 0 1 grid, 1 two stages
 )
{
  Double_t 	z,phi,r,zaft,zprev,zstep,ddr,ddphir,ddz,zl,dr,dphir,dz,ddphi,dphi,deltaz,zm1,zm2,zp1,zp2,gradient_r,gradient_phir,gradient_z,r0,z0,phi0;
  Float_t  x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];

	
  Int_t roc;
  Int_t ianchor,kanchor, zanchor;

	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
	

	
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;
  /// get 

	



  // correction build up for inverse flow
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
	
  TMatrixD * listR;
  TMatrixD * listPhi;
  TMatrixD * listZ;


  /// make interpolator for inverse flow
  ///
  TMatrixD *Mat_corrDrDz[phiSlice];
  TMatrixD *Mat_corrDphiRDz[phiSlice];
  TMatrixD *Mat_corrDz[phiSlice];
	
  TMatrixD *Mat_RList[phiSlice];
  TMatrixD *Mat_PhiList[phiSlice];
  TMatrixD *Mat_ZList[phiSlice];
	
  for (Int_t m=0;m<phiSlice;m++) {
    Mat_corrDrDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDphiRDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDz[m] = new TMatrixD(rrow,zcolumn);

    Mat_RList[m] = new TMatrixD(rrow,zcolumn);
    Mat_PhiList[m] = new TMatrixD(rrow,zcolumn);
    Mat_ZList[m] = new TMatrixD(rrow,zcolumn);
  }
	
  AliTPCLookUpTable3DInterpolatorDFull * lookupInverseCorr =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     rrow,
     Mat_corrDrDz,
     Mat_RList,
     rList,
     phiSlice,
     Mat_corrDphiRDz,
     Mat_PhiList,
     phiList,
     zcolumn,
     Mat_corrDz,
     Mat_ZList,
     zList,
     2,
     stepR,
     stepZ,
     stepPhi,
     type
     );


  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    distDrDz  =  matricesDistDrDz[k] ;
    distDphiRDz  =  matricesDistDphiRDz[k] ;
    distDz  =  matricesDistDz[k] ;
		
    listR = Mat_RList[k];
    listPhi = Mat_PhiList[k];
    listZ  = Mat_ZList[k];

		

    for ( Int_t i = 0 ; i < rrow; i++ ) {
      (*distDrDz)(i,zcolumn-1)    = 0.0;
      (*distDphiRDz)(i,zcolumn-1) = 0.0;
      (*distDz)(i,zcolumn-1)      = 0.0;		
			
      for ( Int_t j = 0 ; j < zcolumn; j++ ) {
	(*listR)(i,j)    = rList[i];
	(*listPhi)(i,j) = phiList[k];
	(*listZ)(i,j)      = zList[j];		
				
      }
				
    }
  }


  // Do for the lowest (17,17,18)

  /////////////////// get local distortion table for coarser grid /////		
  Int_t rrowCoarser = (rrow / 2) + 1;
  Int_t zcolumnCoarser = (zcolumn / 2) + 1;
  Int_t phiSliceCoarser  = phiSlice/2;
  Double_t * rListCoarser = new Double_t[rrowCoarser];
  Double_t * zListCoarser = new Double_t[zcolumnCoarser];
  Double_t * phiListCoarser = new Double_t[phiSliceCoarser];

  // copy rList
  for (int i=0;i<rrowCoarser;i++) rListCoarser[i] = rList[i*2];
  for (int i=0;i<zcolumnCoarser;i++) zListCoarser[i] = zList[i*2];
  for (int i=0;i<phiSliceCoarser;i++) phiListCoarser[i] = phiList[i*2];

	
  TMatrixD *matricesDistDrDzCoarser[phiSliceCoarser];
  TMatrixD *matricesDistDphiRDzCoarser[phiSliceCoarser];
  TMatrixD *matricesDistDzCoarser[phiSliceCoarser];
	


  for (Int_t m=0;m < phiSliceCoarser;m++)   {
    matricesDistDrDzCoarser[m] =  new TMatrixD(rrowCoarser,zcolumnCoarser);
    matricesDistDphiRDzCoarser[m] =  new TMatrixD(rrowCoarser,zcolumnCoarser);
    matricesDistDzCoarser[m] =  new TMatrixD(rrowCoarser,zcolumnCoarser);	
  }


  // prepare interpolator
  AliTPCLookUpTable3DInterpolatorD *lookupLocalDistCoarser = 
    new AliTPCLookUpTable3DInterpolatorD
    (
     rrowCoarser,
     matricesDistDrDzCoarser,
     rListCoarser,
     phiSliceCoarser,
     matricesDistDphiRDzCoarser,
     phiListCoarser,
     zcolumnCoarser,
     matricesDistDzCoarser,
     zListCoarser,
     fInterpolationOrder
     );

  InverseGlobalToLocalDistortion
    (
     matricesDistDrDzCoarser, 
     matricesDistDphiRDzCoarser,  
     matricesDistDzCoarser,
     rListCoarser,
     zListCoarser,
     phiListCoarser,
     rrowCoarser,
     zcolumnCoarser,
     phiSliceCoarser,
     nstep,
     useCylAC,
     stepR,
     stepZ,
     stepPhi,
     type
     );

  // copy vals from mat** to mat* FIX
  lookupLocalDistCoarser->CopyVals();
  //

  // For current grid we have two cases

  deltaz =  (zList[1] - zList[0]);

  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) { 

    printf("inversion global on z column  = %d\n",j);
		 
    roc = 0; // FIXME
    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
				
      distDrDz  =  matricesDistDrDz[k] ;
      distDphiRDz  =  matricesDistDphiRDz[k] ;
      distDz  =  matricesDistDz[k] ;

			
      corrDrDz  =  Mat_corrDrDz[k] ;
      corrDphiRDz  =  Mat_corrDphiRDz[k] ;
      corrDz  =  Mat_corrDz[k] ;

			
      listR  =  Mat_RList[k] ;
      listPhi  =  Mat_PhiList[k] ;
      listZ  =  Mat_ZList[k] ;
		

      for ( Int_t i = 0 ; i < rrow; i++ ) {
	// get global distortion
				
	r = rList[i];
	phi = phiList[k];			
	z = zList[j];
	zp1 = zList[j+1];
				
	dr    = 0.0;
	dz    = 0.0;				
	dphir	= 0.0;
				
	if (j < zcolumn - 2) 
	  {
	    // get global distortion of this point

	    //printf("inversion global on z column  = %d\n",j);
	    //				if (j  % 2 == 0) {
	    // in this case we know local distortion from coarser grid
	    // find local distortion at


	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

				  
	    //						printf("original point (%f,%f,%f)\n",r,phi,z);

	    lookupLocalDistCoarser->GetValue(r,phi,z,dr,dphir,dz);

	    // distorted point
	    r0    = r + dr;
	    z0    = zList[j+2] + dz;				
	    phi0  = phi + (dphir/r);

	    // find local correction at this point
	    //						printf("distorted point (%f,%f,%f)\n",r0,phi0,z0);

	    if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    ianchor = TMath::FloorNint((r0 - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR);
	    kanchor = TMath::FloorNint(phi0/gridSizePhi); 
	    zanchor = TMath::FloorNint(z0/gridSizeZ);

						
	    if (j > zcolumn - 5)  
	      lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,1,j+2);
	    else
	      lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,3,j+2);

	    // get the last point to calculate local distortion
	    phi0 = phi0 + ((dphir)/r0);
	    r0 = r0 + ( dr );
	    z0 = z0 - deltaz + (dz);


	    // find local correction at this point
	    //					printf("last distorted point (%f,%f,%f)\n",r0,phi0,z0);

	    // get final distortion
	    x[0] = r0 ;
	    x[1] = phi0;				
	    x[2] = z0;

	    if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);



	    // } else {
				
	
	    //     x[0] = r ;
	    //     x[1] = phi;				
	    //     x[2] = zList[j];
	    //     if (useCylAC == kTRUE)
	    //       GetDistortionCylAC(x,roc,dx);
	    //     else
	    //       GetDistortionCyl(x,roc,dx);

				  
				  
	    //     // get position on 0 based on global distortion
	    //     r0    = r + dx[0];
	    //     z0    = zList[zcolumn - 1] + dx[2];				
	    //     phi0  = phi + (dx[1]/r);
					  
	    //     // follow electron path from correction table
	    //     for (Int_t jj=zcolumn-1;jj > j+1; jj--) 
	    //     {
						

	    // 		if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    // 		if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    // 		ianchor = TMath::FloorNint((r0 - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR);
	    // 		kanchor = TMath::FloorNint(phi0/gridSizePhi); 
	    // 		zanchor = TMath::FloorNint(z0/gridSizeZ);

	    // 		//printf("call %f,%f,%f,%d,%d,%d\n",r0,phi0,z0,ianchor,kanchor,zanchor);				
	
	    // 		if (j > zcolumn - 5)  
	    // 		  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,1,j+2);
	    // 		else
	    // 		  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,3,j+2);

	    // 		phi0 = phi0 + ((dphir)/r0);
	    // 		r0 = r0 + ( dr );

	    // 		z0 = z0 - deltaz + (dz);

	    //   	}
				  
	    // x[0] = r0 ;
	    // x[1] = phi0;				
	    // x[2] = z0;

	    // if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    // if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    // if (useCylAC == kTRUE)
	    // 	GetDistortionCylAC(x,roc,pdx);
	    // else
	    // 	GetDistortionCyl(x,roc,pdx);

	    // dr   = (dx[0] - pdx[0]);
	    // dz   = (dx[2] - pdx[2]);				
	    // dphir = (dx[1] - pdx[1]);
	    // }					

	  } else if (j == (zcolumn -2))
	  {

	    // in case of zcolumn-2, the global distortion is the local distortion
	    //printf("inversion global on z column  = %d\n",j);

	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

	    x[2] = zList[j+1];


	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);
				  

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);
				  
	  } 

	(*distDrDz)(i,j)    = dr;
	(*distDz)(i,j)      = dz;				
	(*distDphiRDz)(i,j) = dphir;

	r = rList[i];
	phi = phiList[k];			
	z = zList[j];

	(*corrDrDz)(i,j+1)    = -dr;
	(*corrDz)(i,j+1)      = -dz;				
	(*corrDphiRDz)(i,j+1) = -dphir; 

	(*listR)(i,j+1)  =  r + dr;
	(*listPhi)(i,j+1)  =  phi + dphir/r;
	(*listZ)(i,j+1) = zp1 + dz;

				
      }

					
    }

    // copy Vals to correction look up table
    // in case RBF should compute weighted to the interpolant
    lookupInverseCorr->CopyVals(j+1);

		

  }

	

  // dealocate memory
  delete lookupLocalDistCoarser;

  for (Int_t m=0;m < phiSliceCoarser;m++)   {
    delete matricesDistDrDzCoarser[m];
    delete matricesDistDphiRDzCoarser[m];
    delete matricesDistDzCoarser[m];	
  }


  delete[] rListCoarser;
  delete[] zListCoarser;
  delete[] phiListCoarser;





  for (Int_t m=0;m<phiSlice;m++) {
    delete Mat_corrDrDz[m];
    delete Mat_corrDphiRDz[m];
    delete Mat_corrDz[m];
    delete Mat_RList[m];
    delete Mat_PhiList[m];
    delete Mat_ZList[m];
  }
  delete lookupInverseCorr;
       
	
  // delete[] Mat_corrDrDz;
  // delete[] Mat_corrDphiRDz;
  // delete[] Mat_corrDz;
       
  // delete[] Mat_RList;
  // delete[] Mat_PhiList;
  // delete[] Mat_ZList;	
}




void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionMultigrid
(
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice,
 const Int_t nstep,
 const Bool_t useCylAC,
 Int_t stepR,
 Int_t stepZ,
 Int_t stepPhi,
 Int_t type
 )
{
  Double_t 	z,phi,r,zaft,zprev,zstep,ddr,ddphir,ddz,zl,dr,dphir,dz,ddphi,dphi,deltaz,zm1,zm2,zp1,zp2,gradient_r,gradient_phir,gradient_z,r0,z0,phi0;
  Float_t  x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];

	
  Int_t roc;
  Int_t ianchor,kanchor, zanchor;

	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
	

	
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;
  /// get 

	



  // correction build up for inverse flow
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
	
  TMatrixD * listR;
  TMatrixD * listPhi;
  TMatrixD * listZ;


  /// make interpolator for inverse flow
  ///
  TMatrixD *Mat_corrDrDz[phiSlice];
  TMatrixD *Mat_corrDphiRDz[phiSlice];
  TMatrixD *Mat_corrDz[phiSlice];
	
  TMatrixD *Mat_RList[phiSlice];
  TMatrixD *Mat_PhiList[phiSlice];
  TMatrixD *Mat_ZList[phiSlice];
	
  for (Int_t m=0;m<phiSlice;m++) {
    Mat_corrDrDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDphiRDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDz[m] = new TMatrixD(rrow,zcolumn);

    Mat_RList[m] = new TMatrixD(rrow,zcolumn);
    Mat_PhiList[m] = new TMatrixD(rrow,zcolumn);
    Mat_ZList[m] = new TMatrixD(rrow,zcolumn);
  }
	
  AliTPCLookUpTable3DInterpolatorDFull * lookupInverseCorr =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     rrow,
     Mat_corrDrDz,
     Mat_RList,
     rList,
     phiSlice,
     Mat_corrDphiRDz,
     Mat_PhiList,
     phiList,
     zcolumn,
     Mat_corrDz,
     Mat_ZList,
     zList,
     2,
     stepR,
     stepZ,
     stepPhi,
     type
     );


  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    distDrDz  =  matricesDistDrDz[k] ;
    distDphiRDz  =  matricesDistDphiRDz[k] ;
    distDz  =  matricesDistDz[k] ;
		
    listR = Mat_RList[k];
    listPhi = Mat_PhiList[k];
    listZ  = Mat_ZList[k];

		

    for ( Int_t i = 0 ; i < rrow; i++ ) {
      (*distDrDz)(i,zcolumn-1)    = 0.0;
      (*distDphiRDz)(i,zcolumn-1) = 0.0;
      (*distDz)(i,zcolumn-1)      = 0.0;		
			
      for ( Int_t j = 0 ; j < zcolumn; j++ ) {
	(*listR)(i,j)    = rList[i];
	(*listPhi)(i,j) = phiList[k];
	(*listZ)(i,j)      = zList[j];		
				
      }
				
    }
  }


  // Do for the lowest (17,17,18)
  Int_t rrowCoarsest = 17;
  Int_t zcolumnCoarsest = 17;
  Int_t phiSliceCoarsest = 18;


  /////////////////// get local distortion table for coarser grid /////		
  Int_t rrowCoarser = rrowCoarsest;
  Int_t zcolumnCoarser = 	zcolumnCoarsest;
  Int_t phiSliceCoarser  = phiSliceCoarsest;
  Double_t * rListCoarser = new Double_t[rrowCoarser];
  Double_t * zListCoarser = new Double_t[zcolumnCoarser];
  Double_t * phiListCoarser = new Double_t[phiSliceCoarser];

  // copy rList
  for (int i=0;i<rrowCoarser;i++) rListCoarser[i] = rList[i*2];
  for (int i=0;i<zcolumnCoarser;i++) zListCoarser[i] = zList[i*2];
  for (int i=0;i<phiSliceCoarser;i++) phiListCoarser[i] = phiList[i*2];

	
  TMatrixD *matricesDistDrDzCoarser[phiSliceCoarser];
  TMatrixD *matricesDistDphiRDzCoarser[phiSliceCoarser];
  TMatrixD *matricesDistDzCoarser[phiSliceCoarser];
	


  for (Int_t m=0;m < phiSliceCoarser;m++)   {
    matricesDistDrDzCoarser[m] =  new TMatrixD(rrowCoarser,zcolumnCoarser);
    matricesDistDphiRDzCoarser[m] =  new TMatrixD(rrowCoarser,zcolumnCoarser);
    matricesDistDzCoarser[m] =  new TMatrixD(rrowCoarser,zcolumnCoarser);	
  }


  // prepare interpolator
  AliTPCLookUpTable3DInterpolatorD *lookupLocalDistCoarser = 
    new AliTPCLookUpTable3DInterpolatorD
    (
     rrowCoarser,
     matricesDistDrDzCoarser,
     rListCoarser,
     phiSliceCoarser,
     matricesDistDphiRDzCoarser,
     phiListCoarser,
     zcolumnCoarser,
     matricesDistDzCoarser,
     zListCoarser,
     fInterpolationOrder
     );

  InverseGlobalToLocalDistortion
    (
     matricesDistDrDzCoarser, 
     matricesDistDphiRDzCoarser,  
     matricesDistDzCoarser,
     rListCoarser,
     zListCoarser,
     phiListCoarser,
     rrowCoarser,
     zcolumnCoarser,
     phiSliceCoarser,
     nstep,
     useCylAC,
     stepR,
     stepZ,
     stepPhi,
     type
     );

  // copy vals from mat** to mat* FIX
  lookupLocalDistCoarser->CopyVals();
  //

  // For current grid we have two cases

  deltaz =  (zList[1] - zList[0]);

  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) { 

    printf("inversion global on z column  = %d\n",j);
		 
    roc = 0; // FIXME
    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
				
      distDrDz  =  matricesDistDrDz[k] ;
      distDphiRDz  =  matricesDistDphiRDz[k] ;
      distDz  =  matricesDistDz[k] ;

			
      corrDrDz  =  Mat_corrDrDz[k] ;
      corrDphiRDz  =  Mat_corrDphiRDz[k] ;
      corrDz  =  Mat_corrDz[k] ;

			
      listR  =  Mat_RList[k] ;
      listPhi  =  Mat_PhiList[k] ;
      listZ  =  Mat_ZList[k] ;
		

      for ( Int_t i = 0 ; i < rrow; i++ ) {
	// get global distortion
				
	r = rList[i];
	phi = phiList[k];			
	z = zList[j];
	zp1 = zList[j+1];
				
	dr    = 0.0;
	dz    = 0.0;				
	dphir	= 0.0;
				
	if (j < zcolumn - 2) 
	  {
	    // get global distortion of this point

	    //printf("inversion global on z column  = %d\n",j);
	    //					if (j  % 2 == 0) {
	    // in this case we know local distortion from coarser grid
	    // find local distortion at


	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

				  
	    printf("original point (%f,%f,%f)\n",r,phi,z);

	    lookupLocalDistCoarser->GetValue(r,phi,z,dr,dphir,dz);

	    // distorted point
	    r0    = r + dr;
	    z0    = zList[j+2] + dz;				
	    phi0  = phi + (dphir/r);

	    // find local correction at this point
	    printf("distorted point (%f,%f,%f)\n",r0,phi0,z0);

	    if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    ianchor = TMath::FloorNint((r0 - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR);
	    kanchor = TMath::FloorNint(phi0/gridSizePhi); 
	    zanchor = TMath::FloorNint(z0/gridSizeZ);

						
	    if (j > zcolumn - 5)  
	      lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,1,j+2);
	    else
	      lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,3,j+2);

	    // get the last point to calculate local distortion
	    phi0 = phi0 + ((dphir)/r0);
	    r0 = r0 + ( dr );
	    z0 = z0 - deltaz + (dz);


	    // find local correction at this point
	    printf("last distorted point (%f,%f,%f)\n",r0,phi0,z0);

	    // get final distortion
	    x[0] = r0 ;
	    x[1] = phi0;				
	    x[2] = z0;

	    if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);


	    // } else {
				
	
	    //     x[0] = r ;
	    //     x[1] = phi;				
	    //     x[2] = zList[j];
	    //     if (useCylAC == kTRUE)
	    //       GetDistortionCylAC(x,roc,dx);
	    //     else
	    //       GetDistortionCyl(x,roc,dx);

				  
				  
	    //     // get position on 0 based on global distortion
	    //     r0    = r + dx[0];
	    //     z0    = zList[zcolumn - 1] + dx[2];				
	    //     phi0  = phi + (dx[1]/r);
					  
	    //     // follow electron path from correction table
	    //     for (Int_t jj=zcolumn-1;jj > j+1; jj--) 
	    //     {
						

	    // 		if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    // 		if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    // 		ianchor = TMath::FloorNint((r0 - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR);
	    // 		kanchor = TMath::FloorNint(phi0/gridSizePhi); 
	    // 		zanchor = TMath::FloorNint(z0/gridSizeZ);

	    // 		//printf("call %f,%f,%f,%d,%d,%d\n",r0,phi0,z0,ianchor,kanchor,zanchor);				
	
	    // 		if (j > zcolumn - 5)  
	    // 		  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,1,j+2);
	    // 		else
	    // 		  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,3,j+2);

	    // 		phi0 = phi0 + ((dphir)/r0);
	    // 		r0 = r0 + ( dr );

	    // 		z0 = z0 - deltaz + (dz);

	    //   	}
				  
	    // x[0] = r0 ;
	    // x[1] = phi0;				
	    // x[2] = z0;

	    // if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    // if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    // if (useCylAC == kTRUE)
	    // 	GetDistortionCylAC(x,roc,pdx);
	    // else
	    // 	GetDistortionCyl(x,roc,pdx);

	    // dr   = (dx[0] - pdx[0]);
	    // dz   = (dx[2] - pdx[2]);				
	    // dphir = (dx[1] - pdx[1]);
	    // }					
	  } else if (j == (zcolumn -2))
	  {

	    // in case of zcolumn-2, the global distortion is the local distortion
	    //printf("inversion global on z column  = %d\n",j);

	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

	    x[2] = zList[j+1];


	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);
				  

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);
				  
	  } 

	(*distDrDz)(i,j)    = dr;
	(*distDz)(i,j)      = dz;				
	(*distDphiRDz)(i,j) = dphir;

	r = rList[i];
	phi = phiList[k];			
	z = zList[j];

	(*corrDrDz)(i,j+1)    = -dr;
	(*corrDz)(i,j+1)      = -dz;				
	(*corrDphiRDz)(i,j+1) = -dphir; 

	(*listR)(i,j+1)  =  r + dr;
	(*listPhi)(i,j+1)  =  phi + dphir/r;
	(*listZ)(i,j+1) = zp1 + dz;

				
      }

					
    }

    // copy Vals to correction look up table
    // in case RBF should compute weighted to the interpolant
    lookupInverseCorr->CopyVals(j+1);

		

  }

	

  // dealocate memory
  delete lookupLocalDistCoarser;

  for (Int_t m=0;m < phiSliceCoarser;m++)   {
    delete matricesDistDrDzCoarser[m];
    delete matricesDistDphiRDzCoarser[m];
    delete matricesDistDzCoarser[m];	
  }


  delete[] rListCoarser;
  delete[] zListCoarser;
  delete[] phiListCoarser;





  for (Int_t m=0;m<phiSlice;m++) {
    delete Mat_corrDrDz[m];
    delete Mat_corrDphiRDz[m];
    delete Mat_corrDz[m];
    delete Mat_RList[m];
    delete Mat_PhiList[m];
    delete Mat_ZList[m];
  }
  delete lookupInverseCorr;
       
  // delete[] Mat_corrDrDz;
  // delete[] Mat_corrDphiRDz;
  // delete[] Mat_corrDz;
       
  // delete[] Mat_RList;
  // delete[] Mat_PhiList;
  // delete[] Mat_ZList;
}



/// 
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param rrow Int_t number of grid in r direction
/// \param zcolumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nstep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionGlobalInvTable
(
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice,
 const Int_t nstep,
 const Bool_t useCylAC,
 Int_t stepR,
 Int_t stepZ,
 Int_t stepPhi,
 Int_t type // 0 1 grid, 1 two stages
 )
{

  Double_t z,phi,r,zaft,zprev,zstep,ddr,ddphir,ddz,zl,dr,dphir,dz,ddphi,dphi,deltaz,zm1,zm2,zp1,zp2,gradient_r,gradient_phir,gradient_z,r0,z0,phi0;
  Float_t  x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];

	
  Int_t roc;
	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
	

	
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;

  // correction build up for inverse flow
  TMatrixD * corrDrDz;
  TMatrixD * corrDphiRDz;
  TMatrixD * corrDz;
	
  TMatrixD * listR;
  TMatrixD * listPhi;
  TMatrixD * listZ;


  /// make interpolator for inverse flow
  ///
  TMatrixD *Mat_corrDrDz[phiSlice];
  TMatrixD *Mat_corrDphiRDz[phiSlice];
  TMatrixD *Mat_corrDz[phiSlice];
	
  TMatrixD *Mat_RList[phiSlice];
  TMatrixD *Mat_PhiList[phiSlice];
  TMatrixD *Mat_ZList[phiSlice];
	
  for (Int_t m=0;m<phiSlice;m++) {
    Mat_corrDrDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDphiRDz[m] = new TMatrixD(rrow,zcolumn);
    Mat_corrDz[m] = new TMatrixD(rrow,zcolumn);

    Mat_RList[m] = new TMatrixD(rrow,zcolumn);
    Mat_PhiList[m] = new TMatrixD(rrow,zcolumn);
    Mat_ZList[m] = new TMatrixD(rrow,zcolumn);
  }
	
  AliTPCLookUpTable3DInterpolatorDFull * lookupInverseCorr =
    new AliTPCLookUpTable3DInterpolatorDFull
    (
     rrow,
     Mat_corrDrDz,
     Mat_RList,
     rList,
     phiSlice,
     Mat_corrDphiRDz,
     Mat_PhiList,
     phiList,
     zcolumn,
     Mat_corrDz,
     Mat_ZList,
     zList,
     2,
     stepR,
     stepZ,
     stepPhi,
     type
     );
  lookupInverseCorr->SetKernelType(GetRBFKernelType());

  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    distDrDz  =  matricesDistDrDz[k] ;
    distDphiRDz  =  matricesDistDphiRDz[k] ;
    distDz  =  matricesDistDz[k] ;
		
    listR = Mat_RList[k];
    listPhi = Mat_PhiList[k];
    listZ  = Mat_ZList[k];

		

    for ( Int_t i = 0 ; i < rrow; i++ ) {
      (*distDrDz)(i,zcolumn-1)    = 0.0;
      (*distDphiRDz)(i,zcolumn-1) = 0.0;
      (*distDz)(i,zcolumn-1)      = 0.0;		
			
      for ( Int_t j = 0 ; j < zcolumn; j++ ) {
	(*listR)(i,j)    = rList[i];
	(*listPhi)(i,j) = phiList[k];
	(*listZ)(i,j)      = zList[j];		
				
      }
				
    }
  }

	

  // 1) create global correction 
  deltaz =  (zList[1] - zList[0]);
  Int_t ianchor,kanchor, zanchor;

  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) { 

    printf("create inversion global inversion on z column  = %d\n",j);
		 
    roc = 0; // FIXME
    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
				

			
      corrDrDz  =  Mat_corrDrDz[k] ;
      corrDphiRDz  =  Mat_corrDphiRDz[k] ;
      corrDz  =  Mat_corrDz[k] ;

			
      listR  =  Mat_RList[k] ;
      listPhi  =  Mat_PhiList[k] ;
      listZ  =  Mat_ZList[k] ;
		

      for ( Int_t i = 0 ; i < rrow; i++ ) {
	// get global distortion
				
	r = rList[i];
	phi = phiList[k];			
	z = zList[j];
				
	dr    = 0.0;
	dz    = 0.0;				
	dphir	= 0.0;
				
  	
	// in case of zcolumn-2, the global distortion is the local distortion
	//printf("inversion global on z column  = %d\n",j);

	x[0] = r ;
	x[1] = phi;				
	x[2] = z;

	if (useCylAC == kTRUE)
	  GetDistortionCylAC(x,roc,dx);
	else
	  GetDistortionCyl(x,roc,dx);

				  

	dr   = dx[0];
	dz   = dx[2];				
	dphir = dx[1];
	//	printf("global distortion (%f,%f,%f) => (%f,%f,%f)\n",r,phi,z,dr,dphir,dz);
				
	r = rList[i];
	phi = phiList[k];			
	z = zList[j];

	(*corrDrDz)(i,j+1)    = -dr;
	(*corrDz)(i,j+1)      = -dz;				
	(*corrDphiRDz)(i,j+1) = -dphir; 

	(*listR)(i,j+1)  =  r + dr;
	(*listPhi)(i,j+1)  =  phi + dphir/r;
	(*listZ)(i,j+1) = z + dz;

	//printf("global correction (%f,%f,%f) => (%f,%f,%f)\n",r + dr,phi + dphir/r,z + dz,-dr,-dphir,-dz);


				
      }

					
    }

	
    lookupInverseCorr->CopyVals(j+1);	

		
  }
	


  // 2) calculate local distortion
  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) { 

    printf("calculate local distortion global on z column  = %d\n",j);
		 
    roc = 0; // FIXME
    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
				
      distDrDz  =  matricesDistDrDz[k] ;
      distDphiRDz  =  matricesDistDphiRDz[k] ;
      distDz  =  matricesDistDz[k] ;
			
		

      for ( Int_t i = 0 ; i < rrow; i++ ) {
	// get global distortion
				
	r = rList[i];
	phi = phiList[k];			
	z = zList[j];

				
	dr    = 0.0;
	dz    = 0.0;				
	dphir	= 0.0;
				
	if (j < zcolumn - 2) 
	  {
	    // get global distortion of this point

	    //printf("inversion global on z column  = %d\n",j);

	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = z;
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

	    //					printf("original point: (%f,%f,%f)\n",r,phi,z);				  

				  
	    // get position 
	    r0    = r + dx[0];
	    z0    = zList[j + 1] + dx[2];				
	    phi0  = phi + (dx[1]/r);
				  
	    //				printf("distorted point: (%f,%f,%f)\n",r0,phi0,z0);				  


	    ianchor = TMath::FloorNint((r0 - o2::TPC::AliTPCPoissonSolver::fgkIFCRadius)/gridSizeR);
	    kanchor = TMath::FloorNint(phi0/gridSizePhi); 
	    zanchor = TMath::FloorNint(z0/gridSizeZ);


	    if (j > zcolumn - (GetIrregularGridSize() + 2))  
	      lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,1,0);
	    else
	      lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,rrow/4  + 1, phiSlice/4 + 1,GetIrregularGridSize(),0);
				  
	    phi0 = phi0 + ((dphir)/r0);
	    r0 = r0 + ( dr );

	    z0 += dz;


	    //			printf("correcting: (%f,%f,%f)\n",dr,dphir,dz);				  
	    //				printf("corrected point: (%f,%f,%f)\n",r0,phi0,z0);				  

				  
	    x[0] = r0 ;
	    x[1] = phi0;				
	    x[2] = z0;

	    if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
	    if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);
		
	    //			printf("local distortion: (%f,%f,%f)\n",dr,dz,dphir);				  
	  } else if (j == (zcolumn -2))
	  {

	    // in case of zcolumn-2, the global distortion is the local distortion
	    //printf("inversion global on z column  = %d\n",j);

	    x[0] = r ;
	    x[1] = phi;				
	    x[2] = zList[j];
	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,dx);
	    else
	      GetDistortionCyl(x,roc,dx);

	    x[2] = zList[j+1];


	    if (useCylAC == kTRUE)
	      GetDistortionCylAC(x,roc,pdx);
	    else
	      GetDistortionCyl(x,roc,pdx);
				  

	    dr   = (dx[0] - pdx[0]);
	    dz   = (dx[2] - pdx[2]);				
	    dphir = (dx[1] - pdx[1]);
				  
	  } 

	(*distDrDz)(i,j)    = dr;
	(*distDz)(i,j)      = dz;				
	(*distDphiRDz)(i,j) = dphir;

	//				printf("local distortion: (%d,%d,%d) => (%f,%f,%f)\n",i,j,k,dr,dz,dphir);				  
				
      }					
    }	

  }






  for (Int_t m=0;m<phiSlice;m++) {
    delete Mat_corrDrDz[m];
    delete Mat_corrDphiRDz[m];
    delete Mat_corrDz[m];
    delete Mat_RList[m];
    delete Mat_PhiList[m];
    delete Mat_ZList[m];
  }
  delete lookupInverseCorr;
       
  // delete[] Mat_corrDrDz;
  // delete[] Mat_corrDphiRDz;
  // delete[] Mat_corrDz;
       
  // delete[] Mat_RList;
  // delete[] Mat_PhiList;
  // delete[] Mat_ZList;

}







void AliTPCSpaceCharge3DDriftLine::InverseLocalDistortionToElectricField 
(
 TMatrixD ** matricesEr, 
 TMatrixD ** matricesEphi,  
 TMatrixD ** matricesEz,
 TMatrixD ** matricesInvLocalIntErDz,
 TMatrixD ** matricesInvLocalIntEphiDz,
 TMatrixD ** matricesInvLocalIntEz,
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice
 )
{
  // calculate integral	
  Float_t localIntErOverEz,localIntEphiOverEz,localIntDeltaEz,z2,gridzsize;
  Double_t r;	
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	



  TMatrixD * distDrDz;
  TMatrixD * distDz;
  TMatrixD * distDphiRDz;
	
  TMatrixD * tdistDz;
  TMatrixD * tdistDphiRDz;
  TMatrixD * tdistDrDz;
	
  Float_t c02c12 =  fC0 * fC0 + fC1 * fC1;

  // solve local integration
  for ( Int_t j = 0 ; j < zcolumn; j++ ) { 
    for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
      distDrDz  =  matricesDistDrDz[k] ;
      distDz  =  matricesDistDz[k] ;
      distDphiRDz  =  matricesDistDphiRDz[k] ;
			
      tdistDrDz  =  matricesInvLocalIntErDz[k] ;
      tdistDz  =  matricesInvLocalIntEz[k] ;
      tdistDphiRDz  =  matricesInvLocalIntEphiDz[k] ;
			
      for ( Int_t i = 0 ; i < rrow; i++ ) {
				
				
	localIntErOverEz = fC0 * (*distDrDz)(i,j) - fC1 * (*distDphiRDz)(i,j);
	localIntErOverEz = localIntErOverEz / (fC0 * fC0 + fC1 * fC1);

	localIntEphiOverEz =  ((*distDrDz)(i,j) - (fC0 * localIntErOverEz)) / fC1;
	localIntDeltaEz = (*distDz)(i,j) / ( o2::TPC::AliTPCPoissonSolver::fgkdvdE * o2::TPC::AliTPCPoissonSolver::fgkdvdE ); // two times?				
				
				
	// solve local integral
	//
	(*tdistDrDz)(i,j)    =  localIntErOverEz;
	(*tdistDphiRDz)(i,j) =  localIntEphiOverEz;
	(*tdistDz)(i,j)      = 	localIntDeltaEz;

      }
    }
  }

  TMatrixD * mEphi;
  TMatrixD * mEr;
  TMatrixD * mEz;

  // use central-backward-forward difference for calculating Electric field component 
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {		
    mEphi =  matricesEphi[m] ;
    mEr   =  matricesEr[m] ;
    mEz   =  matricesEz[m] ;
		
    distDrDz  =  matricesInvLocalIntErDz[m] ;
    distDphiRDz  =  matricesInvLocalIntEphiDz[m] ;
    distDz  =  matricesInvLocalIntEz[m] ;
    // for boundary-z
    for ( Int_t i = 0; i < rrow ; i++ )  {
      //		  	(*mEr)(i,0)         =   (( 3.0 * (*distDrDz)(i,0) - 1.0 * (*distDrDz)(i,1)) / (2.0 * gridSizeZ)) * -1 * ezField  ;
      //		  	(*mEphi)(i,0)         =   ((  3.0 * (*distDphiRDz)(i,0) - 1.0 * (*distDphiRDz)(i,1)) / (2.0 * gridSizeZ)) * -1 * ezField  ;
      //		  	(*mEz)(i,0)         =   (   3.0 * (*distDz)(i,0) - 1.0 * (*distDz)(i,1)) / (2.0 * gridSizeZ );
      (*mEr)(i,0)         =   ((*distDrDz)(i,0) / gridSizeZ) * -1 * ezField  ;
      (*mEphi)(i,0)         =   ((*distDphiRDz)(i,0) / gridSizeZ) * -1 * ezField  ;
      (*mEz)(i,0)         =   ((*distDz)(i,0) /  gridSizeZ );

      //			(*mEr)(i,0)         =   (( 0.5 * (*distDrDz)(i,2) - 1.5 * (*distDrDz)(i,1)) / gridSizeZ) * -1 * ezField  ;
      //			(*mEphi)(i,0)         =  (( 0.5 * (*distDphiRDz)(i,2) - 1.5 * (*distDphiRDz)(i,1)) / gridSizeZ) * -1 * ezField  ;
      //			(*mEz)(i,0)         =    (  0.5 * (*distDz)(i,2) - 1.5 * (*distDz)(i,2)) / gridSizeZ ;


			
      (*mEr)(i,zcolumn-1)         =   (( -0.5 * (*distDrDz)(i,zcolumn-3) + 1.5 * (*distDrDz)(i,zcolumn-2)) / gridSizeZ) * -1 * ezField  ;
      (*mEphi)(i,zcolumn-1)         =  (( -0.5 * (*distDphiRDz)(i,zcolumn-3) + 1.5 * (*distDphiRDz)(i,zcolumn-2)) / gridSizeZ) * -1 * ezField  ;
      (*mEz)(i,zcolumn-1)         =    (  -0.5 * (*distDz)(i,zcolumn-3) + 1.5 * (*distDz)(i,zcolumn-2)) / gridSizeZ ;
    }

    for ( Int_t i = 0 ; i < rrow; i++) {
      for ( Int_t j = 1; j  < zcolumn-1  ; j++ ) {
	(*mEr)(i,j) =  (( (*distDrDz)(i,j) + (*distDrDz)(i,j-1) )  / (2*gridSizeZ)) * -1 * ezField ; // z direction				
	(*mEphi)(i,j) =  (((*distDphiRDz)(i,j) + (*distDphiRDz)(i,j-1) ) / (2*gridSizeZ)) * -1 * ezField  ; // z direction				
	(*mEz)(i,j) =  ( (*distDz)(i,j) + (*distDz)(i,j-1) ) / (2*gridSizeZ) ; // z direction				
      }
    }
		

    // for boundary-z
    /**
       for ( Int_t i = 0 ; i < rrow; i++) {
       for ( Int_t j = 1; j  <  zcolumn -1; j++ ) {
			

       }
       }
    **/

  }

}



//
// Inverse Electric Field to Charge
// using partial differential
void AliTPCSpaceCharge3DDriftLine::InverseElectricFieldToCharge
(
 TMatrixD ** matricesCharge,
 TMatrixD ** matricesEr, 
 TMatrixD ** matricesEphi,  
 TMatrixD ** matricesEz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice
 )
{
					
	
  Float_t radius;
  Double_t dr,dz,dphi;

  Int_t mplus, mminus, mplus2, mminus2, signplus, signminus  ;
  Int_t symmetry = 0;
  // iterate over phislices
	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
	
  for ( Int_t m = 0 ; m < phiSlice ; m++ ) {		
    mplus  = m + 1;   signplus  = 1 ;
    mminus = m - 1 ;  signminus = 1 ;
    mplus2  = m + 2;   
    mminus2 = m - 2 ;  
    if (symmetry==1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if ( mplus  > phiSlice-1 ) mplus  = phiSlice - 2 ;
      if ( mminus < 0 )           mminus = 1 ;
			
			
    }
    else if (symmetry==-1) {       // Anti-symmetry in phi
      if ( mplus  > phiSlice-1 ) { mplus  = phiSlice - 2 ;  signplus  = -1 ; }
      if ( mminus < 0 )           { mminus = 1 ;	            signminus = -1 ; }
    }
    else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if ( mplus  > phiSlice-1 ) mplus  = m + 1 - phiSlice ;
      if ( mminus < 0 )           mminus = m - 1 + phiSlice ;
      if ( mplus2  > phiSlice-1 ) mplus2  = m + 2 - phiSlice ;
      if ( mminus2 < 0 )           mminus2 = m - 2 + phiSlice ;
    }
		
    TMatrixD &arrayCharge    =  *matricesCharge[m] ;
    TMatrixD &arrayEr	=  *matricesEr[m] ;
    TMatrixD &arrayEz	=  *matricesEz[m] ;
    TMatrixD &arrayEphi	=  *matricesEphi[m] ;
    TMatrixD &arrayEphiM	=  *matricesEphi[mminus] ;
    TMatrixD &arrayEphiP	=  *matricesEphi[mplus] ;
    TMatrixD &arrayEphiM2	=  *matricesEphi[mminus2] ;
    TMatrixD &arrayEphiP2	=  *matricesEphi[mplus2] ;
	
	
    // for non-boundary V
    for ( Int_t i = 2 ; i < rrow-2; i++) {
      radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius  + i*gridSizeR ;
      for ( Int_t j = 2; j < zcolumn-2 ; j++ ) {				
	dr    =  ( -arrayEr(i+2,j) + 8*arrayEr(i+1,j) - 8* arrayEr(i-1,j) +  arrayEr(i-2,j) ) / (12*gridSizeR); // r direction				
	dz    =  ( -arrayEz(i,j+2) + 8*arrayEz(i,j+1) - 8* arrayEz(i,j-1) +  arrayEz(i,j-2) ) / (12*gridSizeZ); // r direction				
	dphi  =  ( -arrayEphiP2(i,j) + 8 * arrayEphiP(i,j) - 8 * arrayEphiM(i,j) + arrayEphiM2(i,j) ) / (12*gridSizePhi) ; // phi// didrection			      		
			
	arrayCharge(i,j) = -1 * (arrayEr(i,j)/radius + dr + dphi/radius + dz);
      }
    }
		
    // for boundary in r
    for ( Int_t j = 2; j < zcolumn-2 ; j++ )  {

      // r near inner radius
      // for index r[0]
      radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius;
      dr 	=  ( -(11.0/6.0)*arrayEr(0,j) + (3.0*arrayEr(1,j)) - (1.5*arrayEr(2,j)) + ((1.0/3.0)*arrayEr(3,j)) )  / gridSizeR ; // forward difference

      //	dr 	=  ( -(1.5)*arrayEr(0,j) + (2.0*arrayEr(1,j)) - (0.5*arrayEr(2,j)) )  / gridSizeR;

      dz    =  ( -arrayEz(0,j+2) + 8*arrayEz(0,j+1) - 8* arrayEz(0,j-1) +  arrayEz(0,j-2) ) / (12.0*gridSizeZ); ; // z direction 	
      dphi  =  ( -arrayEphiP2(0,j) + 8 * arrayEphiP(0,j) - 8 * arrayEphiM(0,j) + arrayEphiM2(0,j) ) / (12.0*gridSizePhi) ;

      arrayCharge(0,j)  = -1 * (arrayEr(0,j)/radius + dr + dphi/radius + dz);
			

      // index use central difference 3-point center
      radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + gridSizeR;
      //	dr 	=  (-arrayEr(3,j)  +6.0*arrayEr(2,j) - 3.0*arrayEr(1,j) - 2*arrayEr(0,j) ) / (6.0*gridSizeR) ; // forward difference
      dr    =  (arrayEr(2,j) - arrayEr(0,j)) / (2.0 * gridSizeR);

      dz    =  ( -arrayEz(1,j+2) + 8*arrayEz(1,j+1) - 8* arrayEz(1,j-1) +  arrayEz(1,j-2) ) / (12*gridSizeZ);   // z direction				
      dphi  =  ( -arrayEphiP2(1,j) + 8 * arrayEphiP(1,j) - 8 * arrayEphiM(1,j) + arrayEphiM2(1,j) ) / (12*gridSizePhi) ;
      arrayCharge(1,j)  = -1 * (arrayEr(1,j)/radius + dr + dphi/radius + dz);
			
			
      // index use centra difference 3-pont center
      radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius  + (rrow - 2)*gridSizeR ;
      //	dr =   (2.0 * arrayEr(rrow - 1,j)  + 3.0*arrayEr(rrow - 2,j) - 6.0*arrayEr(rrow -3,j) + arrayEr(rrow-4,j) ) / (6.0*gridSizeR) ; 
      dr    =  (arrayEr(rrow-1,j) - arrayEr(rrow-3,j)) / (2.0 * gridSizeR);

      dz    =  ( -arrayEz(rrow-2,j+2) + 8*arrayEz(rrow-2,j+1) - 8* arrayEz(rrow-2,j-1) +  arrayEz(rrow-2,j-2) ) / (12*gridSizeZ);  
      dphi  =  ( -arrayEphiP2(rrow-2,j) + 8 * arrayEphiP(rrow-2,j) - 8 * arrayEphiM(rrow-2,j) + arrayEphiM2(rrow-2,j) ) / (12.0*gridSizePhi) ;			
      arrayCharge(rrow-2,j)  = -1 * (arrayEr(rrow-2,j)/radius + dr + dphi/radius + dz);			
			
      // index r[rrow -1] backward difference
      radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius  + (rrow - 1)*gridSizeR ;
      //dr =  ( 1.5*arrayEr(rrow-1,j) - 2.0*arrayEr(rrow-2,j) + 0.5*arrayEr(rrow-3,j) ) / gridSizeR ; // backward difference
      dr 	=  ( -(11.0/6.0)*arrayEr(rrow-1,j) + (3.0*arrayEr(rrow-2,j)) - (1.5*arrayEr(rrow-3,j)) + ((1.0/3.0)*arrayEr(rrow-4,j)) )  / (-1 *gridSizeR) ;

      //dz    =  ( arrayEz(rrow-1,j+1) - arrayEz(rrow-1,j-1) ) / (2*gridSizeZ) ; // z direction				
      dz    =  ( -arrayEz(rrow-1,j+2) + 8*arrayEz(rrow-1,j+1) - 8* arrayEz(rrow-1,j-1) +  arrayEz(rrow-1,j-2) ) / (12*gridSizeZ);  

			
      dphi  =  ( -arrayEphiP2(rrow-1,j) + 8 * arrayEphiP(rrow-1,j) - 8 * arrayEphiM(rrow-1,j) + arrayEphiM2(rrow-1,j) ) / (12*gridSizePhi) ;	
      arrayCharge(rrow-1,j)  = -1 * (arrayEr(rrow-1,j)/radius + dr + dphi/radius + dz);			
    }
			       
    // boundary z
    for ( Int_t i = 2; i < rrow-2 ; i++ )  {
      // z[0]
      radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
			

      dz 	=  ( -(11.0/6.0)*arrayEz(i,0) + (3.0*arrayEz(i,1)) - (1.5*arrayEz(i,2)) + ((1.0/3.0)*arrayEz(i,3)) )  / (1 *gridSizeZ) ; // forward difference


      dr    =  ( -arrayEr(i+2,0) + 8*arrayEr(i+1,0) - 8* arrayEr(i-1,0) +  arrayEr(i-2,0) ) / (12*gridSizeR);  ; // z direction				
      dphi  =  ( -arrayEphiP2(i,0) + 8 * arrayEphiP(i,0) - 8 * arrayEphiM(i,0) + arrayEphiM2(i,0) ) / (12*gridSizePhi) ;
      arrayCharge(i,0)  = -1 * (arrayEr(i,0)/radius + dr + dphi/radius + dz);
			
      //	z[1] central difference 3 stencil
      //			dz 	=  (-arrayEz(i,3)  +6.0*arrayEz(i,2) - 3.0*arrayEz(i,1) - 2*arrayEz(i,0) ) / (6.0*gridSizeZ) ; // forward difference			
      dz 	=  (arrayEz(i,2) - arrayEz(i,0))  / (2.0*gridSizeZ) ; // forward difference			

      dr    =  ( -arrayEr(i+2,1) + 8*arrayEr(i+1,1) - 8* arrayEr(i-1,1) +  arrayEr(i-2,1) ) / (12*gridSizeR);  ; // z direction				
      dphi  =  ( -arrayEphiP2(i,1) + 8 * arrayEphiP(i,1) - 8 * arrayEphiM(i,1) + arrayEphiM2(i,1) ) / (12*gridSizePhi) ;
      arrayCharge(i,1)  = -1 * (arrayEr(i,1)/radius + dr + dphi/radius + dz);
			
			
      // z[zcolumn - 2]
      //dz =   (2.0 * arrayEz(i,zcolumn - 1)  + 3.0*arrayEz(i,zcolumn - 2) - 6.0*arrayEz(i,zcolumn -3) + arrayEz(i,zcolumn-4) ) / (6.0*gridSizeZ) ; // backward difference			
      dz 	=  (arrayEz(i,zcolumn-1) - arrayEz(i,zcolumn-3))  / (2.0*gridSizeZ) ; // forward difference			

      dr    =  ( -arrayEr(i+2,zcolumn-2) + 8*arrayEr(i+1,zcolumn-2) - 8* arrayEr(i-1,zcolumn-2) +  arrayEr(i-2,zcolumn-2) ) / (12*gridSizeR);  ; // z direction				
      dphi  =  ( -arrayEphiP2(i,zcolumn-2) + 8 * arrayEphiP(i,zcolumn-2) - 8 * arrayEphiM(i,zcolumn-2) + arrayEphiM2(i,zcolumn-2) ) / (12*gridSizePhi) ;
      arrayCharge(i,zcolumn-2)  = -1 * (arrayEr(i,zcolumn-2)/radius + dr + dphi/radius + dz);
			
			
      // z[zcolumn  - 1]
      //dz =    ( 1.5*arrayEz(i,zcolumn-1) - 2.0*arrayEz(i,zcolumn-2) + 0.5*arrayEz(i,zcolumn-3) ) / gridSizeZ ; // backward 	difference			
      dz 	=  ( -(11.0/6.0)*arrayEz(i,zcolumn-1) + (3.0*arrayEz(i,zcolumn-2)) - (1.5*arrayEz(i,zcolumn-3)) + ((1.0/3.0)*arrayEz(i,zcolumn-4)) )  / (-gridSizeZ) ; // backward difference
      dr    =  ( -arrayEr(i+2,zcolumn-1) + 8*arrayEr(i+1,zcolumn-1) - 8* arrayEr(i-1,zcolumn-1) +  arrayEr(i-2,zcolumn-1) ) / (12*gridSizeR);  ; // z direction				
      dphi  =  ( -arrayEphiP2(i,zcolumn-1) + 8 * arrayEphiP(i,zcolumn-1) - 8 * arrayEphiM(i,zcolumn-1) + arrayEphiM2(i,zcolumn-1) ) / (12*gridSizePhi) ;
			
			
      arrayCharge(i,zcolumn-1)  = -1 * (arrayEr(i,zcolumn-1)/radius + dr + dphi/radius + dz);			
    }
	       
		
		
    // for corner points
    // corner points for Ephi
		
    radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius;
    dr 	=  ( -0.5*arrayEr(2,0) + 2.0*arrayEr(1,0) - 1.5*arrayEr(0,0) ) / gridSizeR ; // forward difference
    dz 	=  ( -0.5*arrayEz(0,2) + 2.0*arrayEz(0,1) - 1.5*arrayEz(0,0) ) / gridSizeZ ; // forward difference

    dphi  =  ( -arrayEphiP2(0,0) + 8 * arrayEphiP(0,0) - 8 * arrayEphiM(0,0) + arrayEphiM2(0,0) ) / (12*gridSizePhi) ;
    arrayCharge(0,0)  = -1 * (arrayEr(0,0)/radius + dr + dphi/radius + dz );
		
		
    dr 	=  ( -0.5*arrayEr(2,1) + 2.0*arrayEr(1,1) - 1.5*arrayEr(0,1) ) / gridSizeR ; // forward difference

    dz 	=  (arrayEz(0,2) - arrayEz(0,0))  / (2.0*gridSizeZ) ; // forward difference			

    //		dz 	=  (-arrayEz(0,3)  +6.0*arrayEz(0,2) - 3.0*arrayEz(0,1) - 2*arrayEz(0,0) ) / (6.0*gridSizeZ) ; // forward difference			
    dphi  =  ( -arrayEphiP2(0,1) + 8 * arrayEphiP(0,1) - 8 * arrayEphiM(0,1) + arrayEphiM2(0,1) ) / (12*gridSizePhi) ;
    arrayCharge(0,1)  = -1 * (arrayEr(0,1)/radius + dr + dphi/radius + dz);		
			
    dr 	=  ( -0.5*arrayEr(2,zcolumn-2) + 2.0*arrayEr(1,zcolumn-2) - 1.5*arrayEr(0,zcolumn-2) ) / gridSizeR ; // forward difference

    dz =   (2.0 * arrayEz(0,zcolumn - 1)  + 3.0*arrayEz(0,zcolumn - 2) - 6.0*arrayEz(0,zcolumn -3) + arrayEz(0,zcolumn-4) ) / (6.0*gridSizeZ) ; // backward difference
			
    dphi  =  ( -arrayEphiP2(0,zcolumn-2) + 8 * arrayEphiP(0,zcolumn-2) - 8 * arrayEphiM(0,zcolumn-2) + arrayEphiM2(0,zcolumn-2) ) / (12*gridSizePhi) ;
    arrayCharge(0,zcolumn-2)  = -1 * (arrayEr(0,zcolumn-2)/radius + dr + dphi/radius + dz );
		
    dr 	=  ( -0.5*arrayEr(2,zcolumn-1) + 2.0*arrayEr(1,zcolumn-1) - 1.5*arrayEr(0,zcolumn-1) ) / gridSizeR ; // forward difference
    dz =   ( 1.5*arrayEz(0,zcolumn-1) - 2.0*arrayEz(0,zcolumn-2) + 0.5*arrayEz(0,zcolumn-3) ) / gridSizeZ ; // backward difference
    dphi  =  ( -arrayEphiP2(0,zcolumn-1) + 8 * arrayEphiP(0,zcolumn-1) - 8 * arrayEphiM(0,zcolumn-1) + arrayEphiM2(0,zcolumn-1) ) / (12*gridSizePhi) ;
    arrayCharge(0,zcolumn-1)  = -1 * (arrayEr(0,zcolumn-1)/radius + dr + dphi/radius + dz);	


    radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + gridSizeR;
    dr 	=  (-arrayEr(3,0)  +6.0*arrayEr(2,0) - 3.0*arrayEr(1,0) - 2*arrayEr(0,0) ) / (6.0*gridSizeR) ; // forward difference			
    dz 	=  ( -0.5*arrayEz(1,2) + 2.0*arrayEz(1,1) - 1.5*arrayEz(1,0) ) / gridSizeZ ; // forward difference
    dphi  =  ( -arrayEphiP2(1,0) + 8 * arrayEphiP(1,0) - 8 * arrayEphiM(1,0) + arrayEphiM2(1,0) ) / (12*gridSizePhi) ;
    arrayCharge(1,0)  = -1 * (arrayEr(1,0)/radius + dr + dphi/radius + dz);
		
    dr 	=  (-arrayEr(3,1)  +6.0*arrayEr(2,1) - 3.0*arrayEr(1,1) - 2*arrayEr(0,1) ) / (6.0*gridSizeR) ; // forward difference			
    dz 	=  (-arrayEz(1,3)  +6.0*arrayEz(1,2) - 3.0*arrayEz(1,1) - 2*arrayEz(1,0) ) / (6.0*gridSizeZ) ; // forward difference			
    dphi  =  ( -arrayEphiP2(1,1) + 8 * arrayEphiP(1,1) - 8 * arrayEphiM(1,1) + arrayEphiM2(1,1) ) / (12*gridSizePhi) ;
    arrayCharge(1,1)  = -1 * (arrayEr(1,1)/radius + dr + dphi/radius + dz);		
			
    dr 	=  (-arrayEr(3,zcolumn-2)  +6.0*arrayEr(2,zcolumn-2) - 3.0*arrayEr(1,zcolumn-2) - 2*arrayEr(0,zcolumn-2) ) / (6.0*gridSizeR) ; // forward difference			
    dz =   (2.0 * arrayEz(1,zcolumn - 1)  + 3.0*arrayEz(1,zcolumn - 2) - 6.0*arrayEz(1,zcolumn -3) + arrayEz(1,zcolumn-4) ) / (6.0*gridSizeZ) ; // backward difference			
    dphi  =  ( -arrayEphiP2(1,zcolumn-2) + 8 * arrayEphiP(1,zcolumn-2) - 8 * arrayEphiM(1,zcolumn-2) + arrayEphiM2(1,zcolumn-2) ) / (12*gridSizePhi) ;
    arrayCharge(1,zcolumn-2)  =-1 * ( arrayEr(1,zcolumn-2)/radius + dr + dphi/radius + dz );
		
    dr 	=  (-arrayEr(3,zcolumn-1)  +6.0*arrayEr(2,zcolumn-1) - 3.0*arrayEr(1,zcolumn-1) - 2*arrayEr(0,zcolumn-1) ) / (6.0*gridSizeR) ; // forward difference 				
    dz =   ( 1.5*arrayEz(1,zcolumn-1) - 2.0*arrayEz(1,zcolumn-2) + 0.5*arrayEz(1,zcolumn-3) ) / gridSizeZ ; // backward difference
		
    dphi  =  ( -arrayEphiP2(1,zcolumn-1) + 8 * arrayEphiP(1,zcolumn-1) - 8 * arrayEphiM(1,zcolumn-1) + arrayEphiM2(1,zcolumn-1) ) / (12*gridSizePhi) ;
		
    arrayCharge(1,zcolumn-1)  = -1 * (arrayEr(1,zcolumn-1)/radius + dr + dphi/radius + dz);	
		


    radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + (rrow - 2)*gridSizeR;
		
		
    dr =   (2.0 * arrayEr(rrow - 1,0)  + 3.0*arrayEr(rrow - 2,0) - 6.0*arrayEr(rrow -3,0) + arrayEr(rrow-4,0) ) / (6.0*gridSizeR) ; // backward difference					
		
    dz 	=  ( -0.5*arrayEz(rrow - 2,2) + 2.0*arrayEz(rrow - 2,1) - 1.5*arrayEz(rrow - 2,0) ) / gridSizeZ ; // forward difference
		
    dphi  =  ( -arrayEphiP2(rrow - 2,0) + 8 * arrayEphiP(rrow - 2,0) - 8 * arrayEphiM(rrow - 2,0) + arrayEphiM2(rrow - 2,0) ) / (12*gridSizePhi) ;
		
    arrayCharge(rrow - 2,0)  = -1 * (arrayEr(rrow - 2,0)/radius + dr + dphi/radius + dz);
		
    dr =   (2.0 * arrayEr(rrow - 1,1)  + 3.0*arrayEr(rrow - 2,1) - 6.0*arrayEr(rrow -3,1) + arrayEr(rrow-4,1) ) / (6.0*gridSizeR) ; // backward difference					
    dz 	=  (-arrayEz(rrow - 2,3)  +6.0*arrayEz(rrow - 2,2) - 3.0*arrayEz(rrow - 2,1) - 2*arrayEz(rrow - 2,0) ) / (6.0*gridSizeZ) ; // forward difference			
    dphi  =  ( -arrayEphiP2(rrow - 2,1) + 8 * arrayEphiP(rrow - 2,1) - 8 * arrayEphiM(rrow - 2,1) + arrayEphiM2(rrow - 2,1) ) / (12*gridSizePhi) ;
    arrayCharge(rrow - 2,1)  = -1 * (arrayEr(rrow - 2,1)/radius + dr + dphi/radius + dz);		
			
    dr =   (2.0 * arrayEr(rrow - 1,zcolumn-2)  + 3.0*arrayEr(rrow - 2,zcolumn-2) - 6.0*arrayEr(rrow -3,zcolumn-2) + arrayEr(rrow-4,zcolumn-2) ) / (6.0*gridSizeR) ; // backward difference						
    dz =   (2.0 * arrayEz(rrow - 2,zcolumn - 1)  + 3.0*arrayEz(rrow - 2,zcolumn - 2) - 6.0*arrayEz(rrow - 2,zcolumn -3) + arrayEz(rrow - 2,zcolumn-4) ) / (6.0*gridSizeZ) ; // backward difference			
    dphi  =  ( -arrayEphiP2(rrow - 2,zcolumn-2) + 8 * arrayEphiP(rrow - 2,zcolumn-2) - 8 * arrayEphiM(rrow - 2,zcolumn-2) + arrayEphiM2(rrow - 2,zcolumn-2) ) / (12*gridSizePhi) ;
    arrayCharge(rrow - 2,zcolumn-2)  = -1 * (arrayEr(rrow - 2,zcolumn-2)/radius + dr + dphi/radius + dz );
		
    dr =   (2.0 * arrayEr(rrow - 1,zcolumn-1)  + 3.0*arrayEr(rrow - 2,zcolumn-1) - 6.0*arrayEr(rrow -3,zcolumn-1) + arrayEr(rrow-4,zcolumn-1) ) / (6.0*gridSizeR) ; // backward difference						
    dz =   ( 1.5*arrayEz(0,zcolumn-1) - 2.0*arrayEz(0,zcolumn-2) + 0.5*arrayEz(0,zcolumn-3) ) / gridSizeZ ; // backward difference
    dphi  =  ( -arrayEphiP2(rrow - 2,zcolumn-1) + 8 * arrayEphiP(rrow - 2,zcolumn-1) - 8 * arrayEphiM(rrow - 2,zcolumn-1) + arrayEphiM2(rrow - 2,zcolumn-1) ) / (12*gridSizePhi) ;
		
		
		
    arrayCharge(rrow - 2,zcolumn-1)  = -1 * (arrayEr(rrow - 2,zcolumn-1)/radius + dr + dphi/radius + dz);	
		
    radius = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius  + (rrow - 1)*gridSizeR ;
    dr =     ( 1.5*arrayEr(rrow-1,0) - 2.0*arrayEr(rrow-2,0) + 0.5*arrayEr(rrow-3,0) ) / gridSizeR ; // backward difference
    dz 	=   ( -0.5*arrayEz(rrow-1,2) + 2.0*arrayEz(rrow-1,1) - 1.5*arrayEz(rrow-1,0) ) / gridSizeZ ; // forward difference
		
    dphi  =  ( -arrayEphiP2(rrow-1,0) + 8 * arrayEphiP(rrow-1,0)- 8 * arrayEphiM(rrow-1,0) + arrayEphiM2(rrow-1,0) ) / (12*gridSizePhi) ;
		
		
    arrayCharge(rrow-1,0)  = -1 * (arrayEr(rrow-1,0)/radius + dr + dphi/radius + dz);
		
    dr =     ( 1.5*arrayEr(rrow-1,1) - 2.0*arrayEr(rrow-2,1) + 0.5*arrayEr(rrow-3,1) ) / gridSizeR ; // backward difference
    dz 	=  (-arrayEz(rrow-1,3)  +6.0*arrayEz(rrow-1,2) - 3.0*arrayEz(rrow-1,1) - 2*arrayEz(rrow-1,0) ) / (6.0*gridSizeZ) ; // forward difference			
    dphi  =  ( -arrayEphiP2(rrow-1,1) + 8 * arrayEphiP(rrow-1,1) - 8 * arrayEphiM(rrow-1,1) + arrayEphiM2(rrow-1,1) ) / (12*gridSizePhi) ;
    arrayCharge(rrow-1,1)  = -1 * (arrayEr(rrow-1,1)/radius + dr + dphi/radius + dz);		
			
    dr =     ( 1.5*arrayEr(rrow-1,zcolumn-2) - 2.0*arrayEr(rrow-2,zcolumn-2) + 0.5*arrayEr(rrow-3,zcolumn-2) ) / gridSizeR ; // backward difference
    dz =   (2.0 * arrayEz(rrow-1,zcolumn - 1)  + 3.0*arrayEz(rrow-1,zcolumn - 2) - 6.0*arrayEz(rrow-1,zcolumn -3) + arrayEz(rrow-1,zcolumn-4) ) / (6.0*gridSizeZ) ; // backward difference			
    dphi  =  ( -arrayEphiP2(rrow-1,zcolumn-2) + 8 * arrayEphiP(rrow-1,zcolumn-2) - 8 * arrayEphiM(rrow-1,zcolumn-2) + arrayEphiM2(rrow-1,zcolumn-2) ) / (12*gridSizePhi) ;
    arrayCharge(rrow-1,zcolumn-2)  = -1 * (arrayEr(rrow-1,zcolumn-2)/radius + dr + dphi/radius + dz);
	
		
    dr =   ( 1.5*arrayEr(rrow-1,zcolumn-1) - 2.0*arrayEr(rrow-2,zcolumn-1) + 0.5*arrayEr(rrow-3,zcolumn-1) ) / gridSizeR ; // backward difference
    dz =   ( 1.5*arrayEz(rrow-1,zcolumn-1) - 2.0*arrayEz(rrow-1,zcolumn-2) + 0.5*arrayEz(rrow-1,zcolumn-3) ) / gridSizeZ ; // backward difference
		
    dphi  =  ( -arrayEphiP2(rrow-1,zcolumn-1) + 8 * arrayEphiP(rrow-1,zcolumn-1)- 8 * arrayEphiM(rrow-1,zcolumn-1) + arrayEphiM2(rrow-1,zcolumn-1) ) / (12*gridSizePhi) ;
		
    arrayCharge(rrow-1,zcolumn-1)  = -1 * (arrayEr(rrow-1,zcolumn-1)/radius + dr + dphi/radius + dz );
		
		
  }
	
}


void AliTPCSpaceCharge3DDriftLine::InverseDistortionMaps(
							 TMatrixD ** matricesCharge,
							 TMatrixD ** matricesEr, 
							 TMatrixD ** matricesEphi,  
							 TMatrixD ** matricesEz,	
							 TMatrixD ** matricesInvLocalIntErDz,
							 TMatrixD ** matricesInvLocalIntEphiDz,
							 TMatrixD ** matricesInvLocalEz,
							 TMatrixD ** matricesDistDrDz, 
							 TMatrixD ** matricesDistDphiRDz,  
							 TMatrixD ** matricesDistDz,
							 const Int_t rrow,	
							 const Int_t zcolumn,	
							 const Int_t phiSlice,
							 const Int_t nsize,
							 const Bool_t useCylAC,
							 Int_t stepR,
							 Int_t stepZ,
							 Int_t stepPhi,
							 Int_t interpType,
							 Int_t inverseType
							 ) 
{
  // can inverse after lookup table for global distortion been calculated
  Double_t * rList = new Double_t[rrow];
  Double_t * zList = new Double_t[zcolumn];
  Double_t * phiList = new Double_t[phiSlice];
	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
		
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) phiList[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow; i++ ) rList[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn; j++ ) zList[j]  = (j * gridSizeZ); 
  // memory alocation
  if (fInitLookUp) {
    // 1)  get local distortion		
    if (inverseType == 2)
      InverseGlobalToLocalDistortionGlobalInvTable(matricesDistDrDz, matricesDistDphiRDz,  matricesDistDz,rList, zList, phiList, rrow,zcolumn,phiSlice,nsize,useCylAC,stepR,stepZ,stepPhi,interpType);
    else if (inverseType == 1)
      InverseGlobalToLocalDistortionTwoStages(matricesDistDrDz, matricesDistDphiRDz,  matricesDistDz,rList, zList, phiList, rrow,zcolumn,phiSlice,nsize,useCylAC,stepR,stepZ,stepPhi,interpType);
    else
      InverseGlobalToLocalDistortion(matricesDistDrDz, matricesDistDphiRDz,  matricesDistDz,rList, zList, phiList, rrow,zcolumn,phiSlice,nsize,useCylAC,stepR,stepZ,stepPhi,interpType);
		

    fLookupInverseDistA->SetLookUpR(matricesDistDrDz);
    fLookupInverseDistA->SetLookUpPhi(matricesDistDphiRDz);
    fLookupInverseDistA->SetLookUpZ(matricesDistDz);
    fLookupInverseDistA->CopyVals();

    // 2)  calculate local integral
    InverseLocalDistortionToElectricField(matricesEr, matricesEphi,  matricesEz, matricesInvLocalIntErDz, matricesInvLocalIntEphiDz, matricesInvLocalEz, 
					  matricesDistDrDz, matricesDistDphiRDz,  matricesDistDz, rList, zList, phiList, rrow,zcolumn,phiSlice);
    // 3)  get potential from electric field assuming zero boundaries
    InverseElectricFieldToCharge(matricesCharge, matricesEr, matricesEphi,  matricesEz,  rList, zList, phiList, rrow,zcolumn,phiSlice);
  }

  // copy charge inverse here just for side A (TODO: do for side C) 
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) 	*(fLookUpChargeInverseA[k]) = *(matricesCharge[k]);
  fInterpolatorInverseChargeA->SetVals(fLookUpChargeInverseA);
  fInterpolatorInverseChargeA->InitCubicSpline();
  // 



  delete[] zList;
  delete[] rList;
  delete[] phiList;
}	
	



/// CalculateEField (New Version: with reorganization of modules)
/// Calculate E field based on look-up table created by Poisson Solver
/// * Differentiate V(r) and solve for E(r) using special equations for the first and last row
/// * Integrate E(r)/E(z) from point of origin to pad plane
/// * Differentiate V(r) and solve for E(phi)
/// * Integrate E(phi)/E(z) from point of origin to pad plane
/// * Differentiate V(r) and solve for E(z) using special equations for the first and last row
/// * Integrate (E(z)-Ezstd) from point of origin to pad plane
///
/// \param arrayofArrayV TMatrixD** 3D matrix representing calculated potential 
/// \param arrayofEroverEz TMatrix** 3D matrix representing e-field at Er/Ez
/// \param arrayofEPhioverEz TMatrix** 3D matrix representing e-field at Ephi/Ez
/// \param arrayofDeltaZ TMatrix** 3D matrix representing e-field at DeltaZ
/// \param rows Int_t number of rows of discritization (in R direction)
/// \param columns Int_t number of columns  of discritization (in Z direction)
/// \param phislices Int_t number of (phislices in phi direction) 
/// \param symmetry Int_t symmetry?
/// \param rocDisplace rocDisplacement
/// 
/// \pre   Matrix arrayofArrayV is assumed had been calculated  by Poisson solver
/// \post  Results of Integration and Derivations for E-field calculation are stored in arrayofEroverEz, arrayofEPhioverEz, arrayofDeltaZ 
///
void AliTPCSpaceCharge3DDriftLine::CalculateEField
( 
 TMatrixD**arrayofArrayV,  
 TMatrixD**arrayofEroverEz, 
 TMatrixD**arrayofEPhioverEz, 
 TMatrixD**arrayofDeltaEz,
 const Int_t rows, 
 const  Int_t columns, 
 const Int_t phislices,
 const   Int_t symmetry, 
 Bool_t rocDisplacement  
  ) 
{
	
	
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rows - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (columns -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phislices;
	
	
  TMatrixD *arrayOfArrayEr[phislices], *arrayOfArrayEz[phislices], *arrayOfArrayEphi[phislices];
  
  			
  
  // AliSysInfo::AddStamp("CalcField", 100,0,0);
	
  //Allocate memory for electric field r,z, phi direction
  for (Int_t k=0;k < phislices ;k++) {
    arrayOfArrayEr[k]   =   new TMatrixD(rows,columns) ;
    arrayOfArrayEz[k] =   new TMatrixD(rows,columns) ;
    arrayOfArrayEphi[k]    =   new TMatrixD(rows,columns) ;
  }
		
  //Differentiate V(r) and solve for E(r) using special equations for the first and last row
  TStopwatch w; 
  w.Start();
		
		
				
  ElectricField
    (
     arrayofArrayV, 
     arrayOfArrayEr, 
     arrayOfArrayEphi, 
     arrayOfArrayEz,
     rows, 
     columns, 
     phislices, 
     gridSizeR, 
     gridSizePhi ,
     gridSizeZ,
     symmetry, 
     o2::TPC::AliTPCPoissonSolver::fgkIFCRadius
     ); 
	
  w.Stop();
  // AliInfo(Form("Time for calculation E-field CPU = %f s\n",w.CpuTime()));
  LOG(INFO) << "Time for calculation E-field CPU = %f s\n" << FairLogger::endl;

  // AliSysInfo::AddStamp("Electron Drift Calc", 120,0,0);
	
  //Integrate E(r)/E(z) from point of origin to pad plane	
	
  IntegrateEz(arrayofEroverEz, arrayOfArrayEr, rows, columns,  phislices, ezField);
  IntegrateEz(arrayofEPhioverEz, arrayOfArrayEphi, rows, columns,  phislices, ezField);
  IntegrateEz(arrayofDeltaEz, arrayOfArrayEz, rows, columns,  phislices, -1.0);


  // calculate z distortion from the integrated Delta Ez residuals
  // and include the aquivalence (Volt to cm) of the ROC shift !!
  for (Int_t m=0; m < phislices;m++) {
    TMatrixD& arrayV    =  *arrayofArrayV[m] ;
    TMatrixD& deltaEz  =  *arrayofDeltaEz[m] ;
		
    for ( Int_t j = 0 ; j < columns ; j++ )  {
      for ( Int_t i = 0 ; i < rows ; i++ ) {
	// Scale the Ez distortions with the drift velocity pertubation -> delivers cm
	deltaEz(i,j) = deltaEz(i,j)*o2::TPC::AliTPCPoissonSolver::fgkdvdE;
	// ROC Potential in cm aquivalent
	Double_t dzROCShift =  arrayV(i, columns -1)/ezField;
	if ( rocDisplacement ) deltaEz(i,j) = deltaEz(i,j) + dzROCShift;  // add the ROC misaligment
      }
    }
  }
  // clear the temporary arrays lists
	
  for ( Int_t k = 0 ; k < phislices ; k++ )  {
    delete arrayOfArrayEr[k];
    delete arrayOfArrayEz[k];
    delete arrayOfArrayEphi[k];
  }			
}


///
/// Integrate at z direction Ez for electron drift calculation
///
///
/// \param arrayofArrayExoverEz TMatrixD** 3D matrix representing ExoverEz
/// \param arrayofArrayEx TMatrix** 3D matrix representing e-field at x direction
/// \param rows const Int_t number of rows of discritization (in R direction)
/// \param columns const Int_t number of columns  of discritization (in Z direction)
/// \param phislices const Int_t number of (phislices in phi direction) 
/// \param ezField const Double_t Electric field in z direction
/// 
/// \pre   arrayofArrayEx is assumed already been calculated by ElectricFieldCalculation
/// \post  Matrix arrayofArrayExoverEz is calculated by integration of arrayofArrayEx
///
void AliTPCSpaceCharge3DDriftLine::IntegrateEz
(
 TMatrixD**arrayofArrayExoverEz, 
 TMatrixD**arrayofArrayEx, 	
 const Int_t rows, 
 const Int_t columns, 
 const Int_t phislices, 
 const Double_t ezField
 ) 
{
		
  const Float_t  gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (columns-1) ;
			
  for ( Int_t m = 0 ; m < phislices ; m++ ) {
    TMatrixD& eXoverEz  =  *arrayofArrayExoverEz[m] ;
    TMatrixD& arrayEx	  =  *arrayofArrayEx[m] ;
		
    for ( Int_t j = columns-1 ; j >= 0 ; j-- ) {  
      for ( Int_t i = 0 ; i < rows ; i++ ) {
				
	/// Calculate integration from int^{0}_{j} (TODO: Split the integration)
	if (j < columns - 3) {
	  eXoverEz(i,j) = eXoverEz(i,j+2) + 	(gridSizeZ/3.0)*(arrayEx(i,j) +  4 * arrayEx(i,j+1) +  	arrayEx(i,j+2))/(-1*ezField);					
	}
	else {
	  if  (j == columns-3)	{
	    eXoverEz(i,j) = (gridSizeZ/3.0)*(arrayEx(i,columns-3)	+ 4*arrayEx(i,columns-2)	+ arrayEx(i,columns-1))/(-1*ezField) ;
	  }
	  if  (j == columns-2) {
	    eXoverEz(i,j) =  (gridSizeZ/3.0)*(1.5*arrayEx(i,columns-2)+ 1.5*arrayEx(i,columns-1))/(-1*ezField) ;
	  }
	  if ( j == columns-1) {
	    eXoverEz(i,j) =  0.0 ;
	  }
	}
      }
    }
  }	
			
}


void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylNoDrift(const Float_t x[],const Short_t roc,Float_t dx[]) {
  /// Calculates the correction due the Space Charge effect within the TPC drift volume

  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    //    InitSpaceCharge3DDistortion();
    return;
  }

  Float_t intEr, intEphi, intdEz ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1] ;
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi()) phi -= TMath::TwoPi();

  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

  // Get the Er and Ephi field integrals plus the integral over DeltaEz
  if (sign == -1 && z < 0.0) {
    printf("call C side\n");
    fLookupIntENoDriftC->GetValue(r,phi,z,intEr,intEphi,intdEz);
  }
  else
    fLookupIntENoDriftA->GetValue(r,phi,z,intEr,intEphi,intdEz);

  // Calculate distorted position
  if ( r > 0.0 ) {
    phi =  phi + fCorrectionFactor *( fC0*intEphi - fC1*intEr ) / r;
    r   =  r   + fCorrectionFactor *( fC0*intEr   + fC1*intEphi );
  }
  Double_t dz = intdEz * fCorrectionFactor * o2::TPC::AliTPCPoissonSolver::fgkdvdE;

  // Calculate correction in cartesian coordinates
  dx[0] = - (r - x[0]);
  dx[1] = - (phi - x[1]);
  dx[2] = - dz;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}

void AliTPCSpaceCharge3DDriftLine::GetDistortionCylNoDrift(const Float_t x[], Short_t roc,Float_t dx[]) {
  /// This function delivers the distortion values dx in respect to the inital coordinates x
  /// roc represents the TPC read out chamber (offline numbering convention)

  GetCorrectionCylNoDrift(x,roc,dx);
  for (Int_t j=0;j<3;++j) dx[j]=-dx[j];
}

// inverse for no drift
// inverse from global distortion to local distortion
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionNoDrift
(
 TMatrixD ** matricesDistDrDz, 
 TMatrixD ** matricesDistDphiRDz,  
 TMatrixD ** matricesDistDz,
 Double_t * rList,
 Double_t * zList,
 Double_t * phiList,
 const Int_t rrow,
 const Int_t zcolumn,
 const Int_t phiSlice
 )
{
  Double_t z,phi,r,zaft,zprev,zstep,ddr,ddphir,ddz,zl,dr,dphir,dz;
  Float_t  x[3], dx[3], pdx[3], dxp1[3], dxp2[3];

  Int_t roc;
  Int_t maxzstep = 1;
  //Int_t izstep;

	
  TMatrixD * distDrDz;
  TMatrixD * distDphiRDz;
  TMatrixD * distDz;

	

  for ( Int_t k = 0 ; k < phiSlice ; k++ ) {
    distDrDz  =  matricesDistDrDz[k] ;
    distDphiRDz  =  matricesDistDphiRDz[k] ;
    distDz  =  matricesDistDz[k] ;

    for ( Int_t i = 0 ; i < rrow; i++ ) {
      (*distDrDz)(i,zcolumn-1)    = 0.0;
      (*distDphiRDz)(i,zcolumn-1) = 0.0;
      (*distDz)(i,zcolumn-1)      = 0.0;
			
      //(*distDrDz)(i,0)    = 0.0;
      //(*distDphiRDz)(i,0) = 0.0;
      //(*distDz)(i,0)      = 0.0;	
				
    }
  }

	
  zstep =  (zList[1] - zList[0]) / maxzstep;
	
  for ( Int_t j = zcolumn-2 ; j >= 0; j-- ) 
    {
	  
      roc = 0; // FIXME
		
      for ( Int_t k = 0 ; k < phiSlice ; k++ ) 
	{
				
	  distDrDz  =  matricesDistDrDz[k] ;
	  distDphiRDz  =  matricesDistDphiRDz[k] ;
	  distDz  =  matricesDistDz[k] ;
	  for ( Int_t i = 0 ; i < rrow; i++ ) 
	    {
	      // get global distortion
				
	      r = rList[i];
	      phi = phiList[k];			
	      z = zList[j];	
	      zprev = zList[j+1];	
	      //zaft = zList[j-1];	
	      zl = fLookupZListA[j] ;  // Symmetric solution in Z that depends only on ABS(Z)
			
	      //printf("(%f,%f)\n",z,zl);
			
				
				
	      (*distDrDz)(i,j)    = 0.0;
	      (*distDphiRDz)(i,j) = 0.0;
	      (*distDz)(i,j)      = 0.0;
	      dr = 0.0;
	      dphir = 0.0;
	      dz = 0.0;
				
	      r = rList[i];
	      phi = phiList[k];			
	      z = zList[j];	
				
	      x[0] = r ;
	      x[1] = phi;				
	      x[2] = z;
				
	      GetDistortionCylNoDrift(x,roc,dx);
				
				
	      //x[0] = x[0] + dr;
	      //x[1] = x[1] + dphir/r;				
	      x[2] = zprev ;
				
				
	      GetDistortionCylNoDrift(x,roc,pdx);
					
				
	      (*distDrDz)(i,j)    = ( dx[0] - pdx[0]);
	      (*distDphiRDz)(i,j) = ( dx[1] - pdx[1]) * r;
	      (*distDz)(i,j)      = ( dx[2] - pdx[2]);
				
				
	      //(*distDrDz)(i,j)    = dr;
	      //(*distDphiRDz)(i,j) = dphir;
	      //(*distDFz)(i,j)      = dz;
				
	      //if (j==1) {
	      //		printf("(step  ) %d => (%f,%f,%f)\n",j,dr,dphir,dz);
	      //		printf("(direct) %d => (%f,%f,%f)\n",j,(*distDrDz)(i,j),(*distDphiRDz)(i,j),(*distDz)(i,j));
	      //}
			   
	    }
			  
			
		      
	}
    }
		
}


void AliTPCSpaceCharge3DDriftLine::InverseDistortionMapsNoDrift(
								TMatrixD ** matricesCharge,
								TMatrixD ** matricesEr, 
								TMatrixD ** matricesEphi,  
								TMatrixD ** matricesEz,	
								TMatrixD ** matricesInvLocalIntErDz,
								TMatrixD ** matricesInvLocalIntEphiDz,
								TMatrixD ** matricesInvLocalEz,
								TMatrixD ** matricesDistDrDz, 
								TMatrixD ** matricesDistDphiRDz,  
								TMatrixD ** matricesDistDz,
								const Int_t rrow,	
								const Int_t zcolumn,	
								const Int_t phiSlice 
								) 
{
  // can inverse after lookup table for global distortion been calculated
  Double_t * rList = new Double_t[rrow];
  Double_t * zList = new Double_t[zcolumn];
  Double_t * phiList = new Double_t[phiSlice];
	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) phiList[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow; i++ ) rList[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn; j++ ) zList[j]  = (j * gridSizeZ); 
  // memory alocation
  if (fInitLookUp) {
    // 1)  get local distortion		
    InverseGlobalToLocalDistortionNoDrift(matricesDistDrDz, matricesDistDphiRDz,  matricesDistDz,rList, zList, phiList, rrow,zcolumn,phiSlice);
    // 2)  calculate local integral
    InverseLocalDistortionToElectricField(matricesEr, matricesEphi,  matricesEz, matricesInvLocalIntErDz, matricesInvLocalIntEphiDz, matricesInvLocalEz, 
					  matricesDistDrDz, matricesDistDphiRDz,  matricesDistDz, rList, zList, phiList, rrow,zcolumn,phiSlice);
    // 3)  get potential from electric field assuming zero boundaries
    InverseElectricFieldToCharge(matricesCharge, matricesEr, matricesEphi,  matricesEz,  rList, zList, phiList, rrow,zcolumn,phiSlice);
  }
  delete[] zList;
  delete[] rList;
  delete[] phiList;
}	
	


void AliTPCSpaceCharge3DDriftLine::GetChargeDensity
(
 TMatrixD** matricesChargeA,		
 TMatrixD** matricesChargeC,		
 TH3 *  spaceChargeHistogram3D,
 const Int_t rrow,		
 const Int_t zcolumn,		
 const Int_t phiSlice
 )
{
  Int_t phiSlicesPerSector = phiSlice / kNumSector;	
  const Float_t gridSizeR   =  (o2::TPC::AliTPCPoissonSolver::fgkOFCRadius-o2::TPC::AliTPCPoissonSolver::fgkIFCRadius) / (rrow - 1) ;
  const Float_t gridSizeZ   =  o2::TPC::AliTPCPoissonSolver::fgkTPCZ0 / (zcolumn -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / phiSlice;
  const Double_t ezField = (o2::TPC::AliTPCPoissonSolver::fgkCathodeV-o2::TPC::AliTPCPoissonSolver::fgkGG)/o2::TPC::AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;	

	
  // local variables
  Float_t radius0, phi0, z0;
	
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t  rlist[rrow], zedlist[zcolumn] , philist[phiSlice];
	
  for ( Int_t k = 0 ; k < phiSlice ; k++ ) philist[k] =  gridSizePhi * k;
  for ( Int_t i = 0 ; i < rrow ; i++ )  rlist[i] = o2::TPC::AliTPCPoissonSolver::fgkIFCRadius + i*gridSizeR ;
  for ( Int_t j = 0 ; j < zcolumn ; j++ ) zedlist[j]  = j * gridSizeZ ;

	
  TMatrixD *mCharge;
  for ( Int_t side = 0 ; side < 2 ; side++ ) {  					
			
    for ( Int_t k = 0 ; k < phiSlice ; k++ )  {				
      if (side == 0)
	mCharge    =  matricesChargeA[k] ;				
      else
	mCharge    =  matricesChargeC[k] ;				
					
      phi0    = philist[k];	
      for ( Int_t i = 0 ; i < rrow ; i++ ) {
	radius0 = rlist[i];					
	for ( Int_t j = 0 ; j < zcolumn ; j++ ) {						
	  z0 = zedlist[j];						
	  if (side==1) z0= -TMath::Abs(zedlist[j]);							
	  if (spaceChargeHistogram3D != NULL) {
	    (*mCharge)(i,j) = InterpolatePhi(spaceChargeHistogram3D,phi0,radius0,z0);
	    //InterpolatePhi(spaceChargeHistogram3D,phi0,radius0,z0);
	  }
	}
      }				
    }	
  }
}



// charge
Double_t AliTPCSpaceCharge3DDriftLine::GetChargeCylAC
(
 const Float_t x[], 
 Short_t roc
 )
{
	
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
	
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

    
  if (z > -1e-6) 
    return fInterpolatorChargeA->GetValue(r,phi,z);
  else 
    return fInterpolatorChargeC->GetValue(r,phi,z);
}


// chargeInverse
Double_t AliTPCSpaceCharge3DDriftLine::GetInverseChargeCylAC
(
 const Float_t x[], 
 Short_t roc
 )
{
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
	
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

    
  if (z > -1e-6) 
    return fInterpolatorInverseChargeA->GetValue(r,phi,z);
  else 
    return fInterpolatorInverseChargeC->GetValue(r,phi,z);
}


void AliTPCSpaceCharge3DDriftLine::GetLocalDistortionCylAC(const Float_t x[], Short_t roc,Float_t dx[]) {
	
	
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
	
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

  if (z > -1e-6) 
    fLookupDistA->GetValue(r,phi,z,dR,dPhiR,dZ);
  else 
    fLookupDistC->GetValue(r,phi,-z,dR,dPhiR,dZ);
	
		
	
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)
	
}




void AliTPCSpaceCharge3DDriftLine::GetInverseLocalDistortionCylAC(const Float_t x[], Short_t roc,Float_t dx[]) {
  if (!fInitLookUp) {
    // AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    LOG(INFO) << "Lookup table was not initialized! Performing the inizialisation now ..." << FairLogger::endl;
    InitSpaceCharge3DPoissonIntegralDz(129,129,144,100,1e-8);
  }
	
	
  Float_t dR, dPhiR, dZ ;
  Double_t r, phi, z ;
  Int_t    sign;

  r      =  x[0];
  phi    =  x[1];
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  if ( phi > TMath::TwoPi() ) phi = phi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
	
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( (roc%36) < 18 ) {
    sign =  1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if ( sign==1  && z <  o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z =  o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE
  if ( sign==-1 && z > -o2::TPC::AliTPCPoissonSolver::fgkZOffSet ) z = -o2::TPC::AliTPCPoissonSolver::fgkZOffSet;    // Protect against discontinuity at CE


  if ( (sign==1 && z<0) || (sign==-1 && z>0) ) // just a consistency check
    // AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!"); // removed
    LOG(ERROR) << "ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!" << FairLogger::endl;

  if (z > -1e-6) 
    fLookupInverseDistA->GetValue(r,phi,z,dR,dPhiR,dZ);
  else 
    fLookupInverseDistC->GetValue(r,phi,-z,dR,dPhiR,dZ);
	
		
	
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)
	
}


// added; needed since there is no AliTPCCorrection inheritance anymore
void AliTPCSpaceCharge3DDriftLine::Search
( 
 Int_t n, 
 const Double_t xArray[], 
 Double_t x, 
 Int_t &low 
  ) 
{
  /// Search an ordered table by starting at the most recently used point

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[n-1] >= xArray[0] ) ascend = 1 ;  // Ascending ordered table if true

  if ( low < 0 || low > n-1 ) {
    low = -1 ; high = n ;
  } else {                                            // Ordered Search phase
    if ( (Int_t)( x >= xArray[low] ) == ascend )  {
      if ( low == n-1 ) return ;
      high = low + 1 ;
      while ( (Int_t)( x >= xArray[high] ) == ascend ) {
	low = high ;
	increment *= 2 ;
	high = low + increment ;
	if ( high > n-1 )  {  high = n ; break ;  }
      }
    } else {
      if ( low == 0 )  {  low = -1 ;  return ;  }
      high = low - 1 ;
      while ( (Int_t)( x < xArray[low] ) == ascend ) {
	high = low ;
	increment *= 2 ;
	if ( increment >= high )  {  low = -1 ;  break ;  }
	else  low = high - increment ;
      }
    }
  }

  while ( (high-low) != 1 ) {                     // Binary Search Phase
    middle = ( high + low ) / 2 ;
    if ( (Int_t)( x >= xArray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x == xArray[n-1] ) low = n-2 ;
  if ( x == xArray[0]   ) low = 0 ;

}


Float_t AliTPCSpaceCharge3DDriftLine::Interpolate
( 
 const Double_t xArray[], 
 const Float_t yArray[],
 Int_t order, 
 Double_t x 
  ) 
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
  /// Float version (in order to decrease the OCDB size)

  Float_t y ;
  if ( order == 2 ) {                // Quadratic Interpolation = 2
    y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ;
    y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ;
    y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}


Double_t AliTPCSpaceCharge3DDriftLine::Interpolate
( 
 const Double_t xArray[], 
 const Double_t yArray[],
 Int_t order, 
 Double_t x 
  ) 
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.

  Double_t y ;
  if ( order == 2 ) {                // Quadratic Interpolation = 2
    y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ;
    y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ;
    y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}


TH2F* AliTPCSpaceCharge3DDriftLine::CreateTH2F(const char *name,const char *title,
					       const char *xlabel,const char *ylabel,const char *zlabel,
					       Int_t nbinsx,Double_t xlow,Double_t xup,
					       Int_t nbinsy,Double_t ylow,Double_t yup) {
  /// Helper function to create a 2d histogramm of given size

  TString hname=name;
  Int_t i=0;
  if (gDirectory) {
    while (gDirectory->FindObject(hname.Data())) {
      hname =name;
      hname+="_";
      hname+=i;
      ++i;
    }
  }
  TH2F *h=new TH2F(hname.Data(),title,
		   nbinsx,xlow,xup,
		   nbinsy,ylow,yup);
  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  h->GetZaxis()->SetTitle(zlabel);
  h->SetStats(0);
  return h;
}
