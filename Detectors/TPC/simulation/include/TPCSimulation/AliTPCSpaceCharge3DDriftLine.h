#ifndef O2_TPC_ALITPCSPACECHARGE3DDRIFTLINE_H
#define O2_TPC_ALITPCSPACECHARGE3DDRIFTLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
 * See cxx source for full Copyright notice                                */

/// \class AliTPCSpaceCharge3DDriftLine
/// \brief This class provides correction-distortion following the drift line          
///  
/// Example usage:
///
/// ~~~{.cpp}
///	 
/// ~~~
///
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date November 11, 2015///

#include "TNamed.h"
#include "TVectorD.h"
#include "TFormula.h"
#include "TPCSimulation/AliTPCPoissonSolver.h"
#include "TPCSimulation/AliTPCLookUpTable3DInterpolatorD.h"
#include "TPCSimulation/AliTPC3DCylindricalInterpolator.h"
#include "TPCSimulation/AliTPCLookUpTable3DInterpolatorDFull.h" 
#include "TPCSimulation/AliTPC3DCylindricalInterpolatorFull.h"

// removed
/*
  #include "AliTPCSpaceCharge3D.h" 
  #include "AliTPCCorrection.h"
*/

class TCollection;
class TTimeStamp;
class TFormula;
class TH3F;
class TH3;
class TH2F;
class TH2;

namespace o2 {
  namespace TPC {

    // class AliTPCSpaceCharge3DDriftLine : public AliTPCCorrection { // removed
    class AliTPCSpaceCharge3DDriftLine  : public TNamed { // added TNamed
    public:
      AliTPCSpaceCharge3DDriftLine();
      AliTPCSpaceCharge3DDriftLine	
      (
       Int_t nrrows, 
       Int_t nzcolumns, 
       Int_t nphislices,
       Int_t interpolationorder,
       Int_t irregularGridSize,
       Int_t strategy,
       Int_t rbfKernelType
       );
	
      virtual ~AliTPCSpaceCharge3DDriftLine();


      void Recreate	(
			 Int_t nrrows, 
			 Int_t nzcolumns, 
			 Int_t nphislices,
			 Int_t interpolationorder,
			 Int_t strategy
			 );



      void InitSpaceCharge3DPoissonIntegralDz
      (
       Int_t rrow, 
       Int_t zcolumn, 
       Int_t phiSlice,
       Int_t maxIteration,
       Double_t stoppingConv
       );
	
      void InitSpaceCharge3DPoisson
      (
       Int_t rrow, 
       Int_t zcolumn, 
       Int_t phiSlice,
       Int_t maxIteration,
       Double_t stoppingConv
       );
	
      void ForceInitSpaceCharge3DPoissonIntegralDz
      (
       Int_t rrow, 
       Int_t zcolumn, 
       Int_t phiSlice,
       Int_t maxIteration,
       Double_t stoppingConv
       );

      void InitSpaceCharge3DPoissonIntegralDz
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
       );
	
      void InitSpaceCharge3DPoissonIntegralDz	
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
       );
	
      // not following the drift
      void InitSpaceCharge3DPoisson
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
       );
	
      void GetDistortionCyl
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[]
       );
	
      void GetDistortionCylAC
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[]
       );
	
      void GetCorrectionCyl
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[]
       );


      void GetCorrectionCylAC
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[]
       );

      // charge
      Double_t GetChargeCylAC
      (
       const Float_t x[], 
       Short_t roc
       );


      // chargeInverse
      Double_t GetInverseChargeCylAC
      (
       const Float_t x[], 
       Short_t roc
       );


      void GetDistortion
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[]
       );
	
      void GetCorrection
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[]
       );


      void SetStrategyType(Int_t strategy) {
	fStrategy = strategy;
      }
	
      TH2F * CreateHistoDistDRinXY(Float_t z,Int_t nx,Int_t ny);
      TH2F * CreateHistoDistDRPhiinXY(Float_t z,Int_t nx,Int_t ny);
      TH2F * CreateHistoDistDZinXY(Float_t z,Int_t nx,Int_t ny);
	
	
      enum {
	kNumSector = 18
      };
	
      ///< Enumeration of strategy in building lookup table
      enum StrategyType {
	kNaive = 0, ///< Naive, calculate each drift line point
	kUseInterpolator = 1,  ///< Use interpolation for 
      };

      TH3 *GetInputSpaceChargeHistogram() {return fSpaceChargeHistogram3D;}
      // override the set input space charge
      void    SetInputSpaceCharge(TH3 * hisSpaceCharge3D, Double_t norm);
      void    SetInputSpaceCharge3D(TH3 * hisSpaceCharge3D){SetInputSpaceCharge(hisSpaceCharge3D, 1);}
      void    SetInputSpaceChargeA(TMatrixD** matricesLookUpCharge) 
      {
	fInterpolatorChargeA->SetVals(matricesLookUpCharge);
      }

      void    SetInputSpaceChargeC(TMatrixD** matricesLookUpCharge) 
      {
	fInterpolatorChargeC->SetVals(matricesLookUpCharge);
      }
	


      TTree* CreateDistortionTree
      (	
       Double_t step
	);

      void SetNRRows(Int_t nrrows) {fNRRows = nrrows;}
      void SetNPhiSlices(Int_t nphislices) {fNPhiSlices = nphislices;}
      void SetNZColumns(Int_t nzcolumns) {fNZColumns = nzcolumns;}
	
      Int_t GetNRRows() {return fNRRows;}
      Int_t GetNPhiSlices() {return fNPhiSlices;}
      Int_t GetNZColumns() {return fNZColumns;}

      TH2F * CreateHistoSCinXY
      (
       Float_t z, 
       Int_t nx, 
       Int_t ny
       );
	
      TH2F * CreateHistoSCinZR
      (
       Float_t phi, 
       Int_t nz, 
       Int_t nr
       );
	
      // setter and getter
      void SetPoissonSolver(AliTPCPoissonSolver *poissonSolver) {fPoissonSolver = poissonSolver;}
      AliTPCPoissonSolver * GetPoissonSolver() { return fPoissonSolver;}
      void SetInterpolationOrder(Int_t order) {fInterpolationOrder = order;}
      Int_t GetInterpolationOrder(){ return fInterpolationOrder;}

      // common setters and getters for tangled ExB effect
      virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
	fT1=t1; fT2=t2;
	const Double_t wt0=t2*omegaTau;     fC0=1./(1.+wt0*wt0);
	const Double_t wt1=t1*omegaTau;     fC1=wt1/(1.+wt1*wt1);
      };
      void SetC0C1(Float_t c0,Float_t c1) {fC0=c0;fC1=c1;} // CAUTION: USE WITH CARE
      Float_t GetC0() const {return fC0;}
      Float_t GetC1() const {return fC1;}

      // setters and getters
      void SetCorrectionFactor(Float_t correctionFactor) {fCorrectionFactor=correctionFactor;}
      Float_t GetCorrectionFactor() const {return fCorrectionFactor;}


      // inverse the calculated distortion map into charge distribution
      //void InverseDistortionMaps(const Int_t rrow,	const Int_t zcolumn,	const Int_t phiSlice );
		
      void InverseDistortionMaps(
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
				 const Int_t phiSliceconst,
				 const Int_t nstep,
				 const Bool_t useCylAC,
				 Int_t stepR,
				 Int_t stepZ,
				 Int_t stepPhi,
				 Int_t interpType,
				 Int_t inverseType
 
				 );
	
	
      void InverseDistortionMapsNoDrift(
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
					);
      void GetCorrectionCylNoDrift
      (
       const Float_t x[],
       const Short_t roc,
       Float_t dx[]
       );

      void GetDistortionCylNoDrift(const Float_t x[], Short_t roc,Float_t dx[]);
      void InverseGlobalToLocalDistortionNoDrift
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
       );
	
      void GetChargeDensity(
			    TMatrixD** matricesChargeA,		
			    TMatrixD** matricesChargeC, 
			    TH3 *  spaceChargeHistogram3D,
			    const Int_t rrow,
			    const Int_t zcolumn,
			    const Int_t phiSlice
			    );
	
	
      void GetInverseLocalDistortionCylAC(const Float_t x[], Short_t roc,Float_t dx[]);
      void GetLocalDistortionCylAC(const Float_t x[], Short_t roc,Float_t dx[]);

      // setter and getter
      void SetIrregularGridSize(Int_t size) {fIrregularGridSize = size;}
      Int_t GetIrregularGridSize() {return fIrregularGridSize;}

      Int_t GetRBFKernelType() {return fRBFKernelType;}

    private:
      static const Int_t kNMaxPhi = 360;
  
      Int_t fNRRows;     ///< the maximum on row-slices so far ~ 2cm slicing
      Int_t fNPhiSlices; ///< the maximum of phi-slices so far = (8 per sector)
      Int_t fNZColumns;  ///< the maximum on column-slices so  ~ 2cm slicing
	
      // added tensor terms from AliTPCCorrection class
      Double_t fT1;         ///< tensor term of wt - T1
      Double_t fT2;         ///< tensor term of wt - T2
      //
      Float_t fC0; ///< coefficient C0                 (compare Jim Thomas's notes for definitions)
      Float_t fC1; ///< coefficient C1                 (compare Jim Thomas's notes for definitions)
      Float_t fCorrectionFactor;       ///< Space Charge Correction factor in comparison to initialized
  
      Bool_t fInitLookUp;                 ///< flag to check if the Look Up table was created

 
      /// list of r-coordinate of grids
      Double_t * fLookupRList; //[fNRRows]
      /// list of \f$ \phi\f$ -coordinate of grids
      Double_t * fLookupPhiList; //[fNPhiSlices]
      /// list of z-coordinate of grids
      Double_t * fLookupZList; //[fNZColumns]
  
  
      /// list of z-coordinate of grids
      Double_t * fLookupZListA; //[fNZColumns]
  
      /// list of z-coordinate of grids
      Double_t * fLookupZListC; //[fNZColumns]
  
      Int_t fStrategy; ///>  Strategy for building  (1-> naive algorithm (O(n^4)), 2->use interpolation (O(n^3))	
      Int_t  fInterpolationOrder; ///>  Order of interpolation (1-> tri linear, 2->Lagrange interpolation order 2)
      Int_t  fIrregularGridSize; ///>  Size of irregular grid cubes for interpolation (min 3(

      Int_t fRBFKernelType;


	
      /// Matrices for storing Global  \f$ R \f$ direction
      TMatrixD *fLookUpIntDistDrEz[kNMaxPhi];    
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fLookUpIntDistDphiREz[kNMaxPhi];    
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fLookUpIntDistDz[kNMaxPhi];   		 



      /// Matrices for storing Global  \f$ R \f$ direction
      TMatrixD *fLookUpIntDistDrEzA[kNMaxPhi];   	 
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fLookUpIntDistDphiREzA[kNMaxPhi];   
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fLookUpIntDistDzA[kNMaxPhi];   		
	

      /// Matrices for storing Global  \f$ R \f$ direction
      TMatrixD *fLookUpIntDistDrEzC[kNMaxPhi];   		  
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fLookUpIntDistDphiREzC[kNMaxPhi];    
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fLookUpIntDistDzC[kNMaxPhi];   			
	

      TMatrixD *fLookUpErOverEz[kNMaxPhi]  ;   
      TMatrixD *fLookUpEphiOverEz[kNMaxPhi];   
      TMatrixD *fLookUpDeltaEz[kNMaxPhi]   ;   

      TMatrixD *fLookUpErOverEzA[kNMaxPhi]  ;   
      TMatrixD *fLookUpEphiOverEzA[kNMaxPhi];   
      TMatrixD *fLookUpDeltaEzA[kNMaxPhi]   ;   


      TMatrixD *fLookUpErOverEzC[kNMaxPhi]  ;   
      TMatrixD *fLookUpEphiOverEzC[kNMaxPhi];   
      TMatrixD *fLookUpDeltaEzC[kNMaxPhi]   ;   


	
	 
      /// Matrices for storing Global  \f$  R \f$ correction
      TMatrixD *fLookUpIntCorrDrEz[kNMaxPhi];   		   
      /// Matrices for storing Global  \f$ \phi R \f$  correction
      TMatrixD *fLookUpIntCorrDphiREz[kNMaxPhi];   
      /// Matrices for storing Global  \f$ X \f$ correctiona
      TMatrixD *fLookUpIntCorrDz[kNMaxPhi];        



	 
      /// Matrices for storing Global  \f$  R \f$ correction
      TMatrixD *fLookUpIntCorrDrEzA[kNMaxPhi];   	  
      /// Matrices for storing Global  \f$ \phi R \f$  correction
      TMatrixD *fLookUpIntCorrDphiREzA[kNMaxPhi];   
      /// Matrices for storing Global  \f$ X \f$ correctiona
      TMatrixD *fLookUpIntCorrDzA[kNMaxPhi];        


		
	
      // Correction with irregular interpolation
      /// Matrices for storing Global  \f$ R \f$ direction
      TMatrixD *fLookUpIntCorrDrEzIrregularA[kNMaxPhi];   	 
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fLookUpIntCorrDphiREzIrregularA[kNMaxPhi];  
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fLookUpIntCorrDzIrregularA[kNMaxPhi];   		



      TMatrixD *fRListIrregularA[kNMaxPhi];   		
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fPhiListIrregularA[kNMaxPhi];   
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fZListIrregularA[kNMaxPhi];   			



      // Correction with irregular interpolation
      /// Matrices for storing Global  \f$ R \f$ direction
      TMatrixD *fLookUpIntCorrDrEzIrregularC[kNMaxPhi];   		
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fLookUpIntCorrDphiREzIrregularC[kNMaxPhi];   
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fLookUpIntCorrDzIrregularC[kNMaxPhi];   			



      TMatrixD *fRListIrregularC[kNMaxPhi];   		 
      /// Matrices for storing Globar \f$ \phi R \f$ Distortion
      TMatrixD *fPhiListIrregularC[kNMaxPhi];   
      /// Matrices for storing Globar \f$ z \f$ Distortion
      TMatrixD *fZListIrregularC[kNMaxPhi];   			


      /// look up for charge densities
      TMatrixD *fLookUpChargeA[kNMaxPhi]; 
      TMatrixD *fLookUpChargeC[kNMaxPhi]; 
      TMatrixD *fLookUpChargeInverseA[kNMaxPhi]; 
      TMatrixD *fLookUpChargeInverseC[kNMaxPhi]; 



      // lookup table for charge densities, just interpolator
      AliTPC3DCylindricalInterpolator *fInterpolatorChargeA;
      AliTPC3DCylindricalInterpolator *fInterpolatorChargeC;
      AliTPC3DCylindricalInterpolator *fInterpolatorInverseChargeA;
      AliTPC3DCylindricalInterpolator *fInterpolatorInverseChargeC;

      //	AliTPCLookUpTable3DInterpolatorDFull *fLookupIntCorrIrregular;
      AliTPCLookUpTable3DInterpolatorDFull *fLookupIntCorrIrregularA;
      AliTPCLookUpTable3DInterpolatorDFull *fLookupIntCorrIrregularC;
	
	
	



	
      AliTPCLookUpTable3DInterpolatorD *fLookupIntDist;
      AliTPCLookUpTable3DInterpolatorD *fLookupIntCorr;

	
	
      /// split A and C for loookup table
      AliTPCLookUpTable3DInterpolatorD *fLookupIntDistA;
      AliTPCLookUpTable3DInterpolatorD *fLookupIntCorrA;
	
	
      AliTPCLookUpTable3DInterpolatorD *fLookupIntDistC;
      AliTPCLookUpTable3DInterpolatorD *fLookupIntCorrC;
	
      //
      AliTPCLookUpTable3DInterpolatorD *fLookupIntENoDriftA;
      AliTPCLookUpTable3DInterpolatorD *fLookupIntENoDriftC;
	
      AliTPCLookUpTable3DInterpolatorD *fLookupIntENoDrift;

      // look up for local distortion
      AliTPCLookUpTable3DInterpolatorD *fLookupDistA;
      AliTPCLookUpTable3DInterpolatorD *fLookupInverseDistA;
      AliTPCLookUpTable3DInterpolatorD *fLookupDistC;
      AliTPCLookUpTable3DInterpolatorD *fLookupInverseDistC;



	
      /// Matrices for storing interpolated space charge distribution for TH3
      TMatrixF *fSCdensityDistribution[kNMaxPhi]; 
	
      TH3 *    fSpaceChargeHistogram3D;      ///< Histogram with the input space charge histogram - used as an optional input

	
      /// Pointer to a poisson solver
      AliTPCPoissonSolver * fPoissonSolver; //-> Pointer to a poisson solver

	
      void ElectricField(
			 TMatrixD** matricesV,  
			 TMatrixD** matricesEr,  
			 TMatrixD** matricesEphi, 
			 TMatrixD** matricesEz,
			 const Int_t rrow, 
			 const Int_t zcolumn, 
			 const Int_t phiSlices,
			 const Float_t gridSizeR, 
			 const Float_t  gridSizePhi ,
			 const Float_t  gridSizeZ,
			 const  Int_t symmetry, 
			 const Float_t  innerRadius 
			 ); 
	
      void LocalDistCorrDz
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
       );



      void LocalDistCorrDz
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
       );
      void LocalDistDz
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
       );
	
      void LocalCorrDz
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
       );
	
	
	
      void IntegrateDistCorrDriftLineDz 
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
       );

      void IntegrateDistCorrDriftLineDz 
      (
		
       TFormula * intErDzTestFunction,
       TFormula * intEphiRDzTestFunction,
       TFormula * intDzTestFunction,
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
       );



      void IntegrateDistCorrDriftLineDz 
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
       );
	
	
	
      void IntegrateDistCorrDriftLineDzOpt2
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
       );

	
      void IntegrateDistDriftLineDz 
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
       );
	
	
      void IntegrateCorrDriftLineDz 
      (
       TMatrixD** matricesDistDrDz,  
       TMatrixD** matricesDistDphiRDz, 
       TMatrixD** matricesDistDz, 
       TMatrixD** matricesGCorrDrDz,  
       TMatrixD** matricesGCorrDphiRDz, 
       TMatrixD** matricesGCorrDz, 
       const Int_t rrow,  
       const Int_t zcolumn, 
       const Int_t phiSlice,
       const Double_t *rlist,
       const Double_t *philist,
       const Double_t *zlist
       );
	
      void FillLookUpTable
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
       );
	
      void FillLookUpTable
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
       );


      void FillLookUpTableA
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
       );
	
	
      void FillLookUpTableC
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
       );
	
      void GetDistCorrFromLookUpTable
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[],
       TMatrixD**lookUpR,
       TMatrixD**lookUpPhiR,
       TMatrixD**lookUpDz
       );
	
      void GetDistCorrFromLookUpTableCyl
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[],
       TMatrixD**lookUpR,
       TMatrixD**lookUpPhiR,
       TMatrixD**lookUpDz
       );
	
	
	
      void GetDistCorrFromLookUpTable
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[],
       AliTPCLookUpTable3DInterpolatorD * lookUp
       );
	
      void GetDistCorrFromLookUpTableCyl
      (
       const Float_t x[], 
       Short_t roc,
       Float_t dx[],
       AliTPCLookUpTable3DInterpolatorD * lookUp
       );
	
	
      Double_t Interpolate3DTableCyl
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
	);
	
      Double_t Interpolate3DTableCyl
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
	);
	
      Double_t InterpolatePhi
      ( 
       const Double_t xArray[], 
       const Int_t ilow,
       const Int_t nx,
       const Float_t yArray[],
       Int_t order, 
       Double_t x 
	);


      Float_t Interpolate3DTableCyl
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
	);
      /// Function for testing
      /// for correctness analysis
      void LocalDistCorrDzExact
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
	 );
	
	
      void IntegrateDistCorrDriftLineDzOpt2 
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
       );
	
      Double_t InterpolatePhi
      ( 
       TH3 * h3,
       const Double_t r,
       const Double_t phi,
       const Double_t z
	);
	
	
      Float_t  GetSpaceChargeDensity
      (
       Float_t r, 
       Float_t phi, 
       Float_t z
       );


	
      void InverseGlobalToLocalDistortion
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
       );


      void InverseGlobalToLocalDistortion
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
       );
    

 
      void InverseGlobalToLocalDistortionMultigrid
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
       );
	


      void InverseGlobalToLocalDistortionTwoStages
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
       );


      void InverseGlobalToLocalDistortionGlobalInvTable
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
       );


      void InverseLocalDistortionToElectricField
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
       );
	
	
      void InverseElectricFieldToCharge
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
       );
	
	
      // no drift line
      void CalculateEField
      ( 
       TMatrixD**arrayofArrayV,  
       TMatrixD**arrayofEroverEz, 
       TMatrixD**arrayofEPhioverEz, 
       TMatrixD**arrayofDeltaEz,
       const Int_t rows, 
       const  Int_t columns, 
       const Int_t phislices,
       const   Int_t symmetry, 
       Bool_t rocDisplacement  = kFALSE
	);
	
	
      void IntegrateEz
      (
       TMatrixD**arrayofArrayExoverEz, 
       TMatrixD**arrayofArrayEx, 	
       const Int_t rows, 
       const Int_t columns, 
       const Int_t phislices, 
       const Double_t ezField
       );

      // added; needed since there is no AliTPCCorrection inheritance anymore
      void Search
      ( 
       Int_t n, 
       const Double_t xArray[], 
       Double_t x, 
       Int_t &low 
	);
	
	
      Float_t Interpolate
      ( 
       const Double_t xArray[], 
       const Float_t yArray[],
       Int_t order, 
       Double_t x 
	);
	
      Double_t Interpolate
      ( 
       const Double_t xArray[], 
       const Double_t yArray[],
       Int_t order, 
       Double_t x 
	);

      TH2F* CreateTH2F(const char *name,const char *title,
		       const char *xlabel,const char *ylabel,const char *zlabel,
		       Int_t nbinsx,Double_t xlow,Double_t xup,
		       Int_t nbinsy,Double_t ylow,Double_t yup);
    

      /// \cond CLASSIMP
      ClassDef(AliTPCSpaceCharge3DDriftLine,1);
      /// \endcond
    };

  }
}

#endif	// O2_TPC_ALITPCSPACECHARGE3DDRIFTLINE_H
