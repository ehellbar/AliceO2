#ifndef O2_TPC_ALITPCPOISSONSOLVER_H
#define O2_TPC_ALITPCPOISSONSOLVER_H

/************************************************************************
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

#include <TNamed.h>
#include "TMatrixD.h"
#include "TMatrixF.h"
#include "TObjArray.h"
// #include "AliTPCParam.h" // removed
#include "TVectorD.h"

// library kernel
// #include "../cuda_imp/relaxation.h" // removed

class TH2F;
class TTimeStamp;
class TCollection;
// class TTreeSRedirector;  // removed
// class AliExternalTrackParam;  // removed
class TTree;
class THnSparse;
// class AliESDVertex;  // removed

namespace o2 {
namespace TPC {

/// \class AliPoissonSolver
/// \brief Provides poisson solver with many strategies (relaxation, multigrid, fft, cpu and gpu)
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Mar 4, 2015


class AliTPCPoissonSolver : public TNamed {
public:
  
  ///< Enumeration of Poisson Solver Strategy Type
  enum StrategyType {
		kRelaxation = 0, ///< S.O.R Cascaded Multigrid
		kMultigrid = 1,  ///< Geometric MG
		kFFT = 2	       ///< Spectral (TODO)
  };

	///< Enumeration of Cycles Type
	enum CycleType {
		kVCycle = 0,		///< V Cycle
		kWCycle = 1,   	///< W Cycle (TODO)
		kFCycle = 2			///< Full Cycle 
	};
	
	///< Fine -> Coarse Grid transfer operator types
	enum GridTransferType {
		kHalf = 0, ///< Half weighting
		kFull = 1, ///< Full weighting
	};
	
	
	///< Smoothing (Relax) operator types
	enum RelaxType {
		kJacobi = 0,			///< Jacobi (5 Stencil 2D, 7 Stencil 3D_
		kWeightedJacobi = 1, ///< (TODO)
		kGaussSeidel = 2  ///< Gauss Seidel 2D (2 Color, 5 Stensil), 3D (7 Stencil)
	};

	
	///< Coarse -> fine  operator types (TODO: Interp and Restrict in one packet, just one enumeration)
	enum InterpType {
		kHalfInterp = 0,		//< Half biliniear interporlation
		kFullInterp = 1   	//< Full biliniear interporlation
	};
	

	///< Parameters choice for Multigrid 	algorithm
	struct MGParameters {
		Bool_t isFull3D; //< True: coarsening in phi direction
		CycleType cycleType;
		GridTransferType gtType;
		RelaxType relaxType;
		Int_t gamma;
		Int_t NPRE;
		Int_t NPOST;
		Int_t NMGCYCLE;
		Int_t MAXLOOP;		
		
		
		// default values
		MGParameters() {
			isFull3D = kFALSE;
			cycleType = kFCycle;
			gtType = kFull; // default full
			relaxType = kGaussSeidel; // default relaxation method 
			NPRE = 2; 
			NPOST = 2;
			NMGCYCLE = 200;
			MAXLOOP = 6;		
			
		}
	};
    	
  AliTPCPoissonSolver();
  AliTPCPoissonSolver(const char *name,const char *title);
  virtual ~AliTPCPoissonSolver();


  void PoissonSolver2D(TMatrixD &arrayV, TMatrixD &chargeDensity,
					   Int_t rows, Int_t columns, Int_t maxIterations);


  void PoissonSolver3D( TMatrixD **arrayofArrayV, TMatrixD **arrayofChargeDensities,
			    Int_t rows, Int_t columns,  Int_t phislices, Int_t maxIterations,  Int_t symmetry);


  // setter and getter
  void SetStrategy(StrategyType strategy) {fStrategy = strategy;}
  StrategyType GetStrategy() {return fStrategy;}
  
  
  
  static const Double_t fgkTPCZ0;       ///< nominal gating grid position
  static const Double_t fgkIFCRadius;   ///< Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static const Double_t fgkOFCRadius;   ///< Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static const Double_t fgkZOffSet;     ///< Offset from CE: calculate all distortions closer to CE as if at this point
  static const Double_t fgkCathodeV;    ///< Cathode Voltage (volts)
  static const Double_t fgkGG;          ///< Gating Grid voltage (volts)
  static const Double_t fgkdvdE;        ///< [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)
  static const Double_t fgkEM;          ///< charge/mass in [C/kg]
  static const Double_t fgke0;          ///< vacuum permittivity [A·s/(V·m)]
  
	static Double_t fgExactErr; ///< Error tolerated 
	static Double_t fgConvErr; ///< Error tolerated 

	void SetExactSolution(TMatrixD **exactSolution,const Int_t fPhiSlices);
	
	
	Int_t fIterations;
  	MGParameters fMgParameters;	
	void SetCycleType(AliTPCPoissonSolver::CycleType cycleType) {
		fMgParameters.cycleType = cycleType;
	}	
private:
  AliTPCPoissonSolver(const AliTPCPoissonSolver &);               // not implemented
  AliTPCPoissonSolver &operator=(const AliTPCPoissonSolver &);    // not implemented

  StrategyType fStrategy;	///< strategy used
  
  
  TMatrixD **fExactSolution;
	TVectorD *fErrorConvergenceNorm2;
	TVectorD *fErrorConvergenceNormInf;

	TVectorD *fError;
	
	Double_t GetMaxExact() {return fMaxExact;};
	

	
  // AliTPCParam *fTPCParam;	///< pointer to AliTPCParam // removed
  

  void PoissonRelaxation2D(TMatrixD &arrayV, TMatrixD &chargeDensity, Int_t rows, Int_t columns, Int_t maxIterations);
  void PoissonRelaxation3D(TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities,Int_t rows, Int_t columns,  Int_t phislices, Int_t maxIterations,  Int_t symmetry);
  
  void PoissonMultigrid2D(TMatrixD &arrayV, TMatrixD &chargeDensity, Int_t rows, Int_t columns);
  void PoissonMultigrid3D2D(TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities,Int_t rows, Int_t columns,  Int_t phislices,   Int_t symmetry);


	void PoissonMultigrid3D(TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities,Int_t rows, Int_t columns,  Int_t phislices,  Int_t symmetry);
	
//  void PoissonFFT2D(TMatrixD &arrayV, TMatrixD &chargeDensity, Int_t rows, Int_t columns, Int_t maxIterations);
//  void PoissonFFT3D(TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities,Int_t rows, Int_t columns,  Int_t phislices, Int_t maxIterations,  Int_t symmetry);

  Int_t IsPowerOfTwo ( Int_t i ) const  ;
  
  
  

  void Relax2D(TMatrixD &curArrayV,TMatrixD &curCharge,const Int_t trows,const Int_t tcolumns,const Float_t h2, const Float_t tempFourth, const Float_t tempRatio,std::vector<float> &coef1,std::vector<float> &coef2);
	
	void Relax3D( TMatrixD **curArrayV,TMatrixD **curCharge,const Int_t trows,const Int_t tcolumns, const Int_t phislices, const Int_t symmetry,const Float_t h2, const Float_t tempRatioZ,std::vector<float> &coef1,std::vector<float> &coef2, std::vector<float> &coef3,std::vector<float> &coef4);
	
	void Residu2D(TMatrixD&residu, TMatrixD &curArrayV,TMatrixD &curCharge,const Int_t trows,const Int_t tcolumns,const Float_t ih2, const Float_t itempFourth, const Float_t tempRatio,std::vector<float> &coef1,std::vector<float> &coef2);
	void Residu3D(TMatrixD **residu,TMatrixD **curArrayV,TMatrixD **curCharge,const Int_t trows,const Int_t tcolumns,const Int_t phislices, const Int_t symmetry, const Float_t ih2, const Float_t tempRatio,std::vector<float> &coef1,std::vector<float> &coef2,std::vector<float> &coef3,std::vector<float> &icoef4);
	
	
	void Restrict2D(TMatrixD &curCharge,TMatrixD &residu,const Int_t trows,const Int_t tcolumns);
	void Restrict3D(TMatrixD **curCharge,TMatrixD **residu,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices);
	
	void RestrictBoundary2D(TMatrixD &curCharge,TMatrixD &residu,const Int_t trows,const Int_t tcolumns);
	void RestrictBoundary3D(TMatrixD **curCharge,TMatrixD **residu,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices);
	
	void AddInterp2D(TMatrixD &curArrayV,TMatrixD &curArrayVC,const Int_t trows,const Int_t tcolumns);
	void AddInterp3D(TMatrixD **curArrayV,TMatrixD **curArrayVC,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices);
	void Interp2D(TMatrixD &curArrayV,TMatrixD &curArrayVC,const Int_t trows,const Int_t tcolumns);
	
	void Interp3D(TMatrixD **curArrayV,TMatrixD **curArrayVC,const Int_t trows,const Int_t tcolumns, const Int_t newphislices, const Int_t oldphislices);
	void VCycle2D(const Int_t rows, const Int_t columns,const Int_t gridfrom,const  Int_t gridto,const Int_t NPRE, const Int_t NPOST,const  Float_t gridSizeR,const  Float_t ratio,std::vector<TMatrixD *> & tvArrayV,std::vector<TMatrixD *> & tvCharge,std::vector<TMatrixD *> & tvResidu);
  void WCycle2D(const Int_t rows, const Int_t columns, const Int_t gridfrom,const  Int_t gridto, const Int_t gamma,
				const Int_t NPRE, const Int_t NPOST, const  Float_t gridSizeR,const  Float_t ratio,std::vector<TMatrixD *> & tvArrayV,
				std::vector<TMatrixD *> & tvCharge,std::vector<TMatrixD *> & tvResidu);


	void VCycle3D(const Int_t rows, const Int_t columns, const Int_t phislices, const Int_t symmetry,const Int_t gridfrom,const  Int_t gridto,const Int_t NPRE, const Int_t NPOST,const Float_t gridSizeR, const Float_t ratioZ,
		std::vector<TMatrixD **> & tvArrayV,	std::vector<TMatrixD **> & tvCharge,std::vector<TMatrixD **> & tvResidu,std::vector<float> &coef1,std::vector<float> &coef2, std::vector<float> &coef3,std::vector<float> &coef4,std::vector<float> &icoef4);
	
	// V Cycle 3D2D CPU
	void VCycle3D2D
	(
		const Int_t rows, 
		const Int_t columns, 
		const Int_t phislices, 
		const Int_t symmetry,
		const Int_t gridfrom,
		const Int_t gridto,
		const Int_t NPRE, 
		const Int_t NPOST,
		const Float_t gridSizeR, 
		const Float_t ratioZ, 
		const Float_t ratioPhi,
		std::vector<TMatrixD **> & tvArrayV,	
		std::vector<TMatrixD **> & tvCharge,
		std::vector<TMatrixD **> & tvResidu,
		std::vector<float> &coef1,
		std::vector<float> &coef2, 
		std::vector<float> &coef3,
		std::vector<float> &coef4,
		std::vector<float> &icoef4
	);


	// V Cycle 3D2D GPU Version
	void VCycle3D2DGPU
	(
		const Int_t rows, 
		const Int_t columns, 
		const Int_t phislices, 
		const Int_t symmetry,
		const Int_t gridfrom,
		const Int_t gridto,
		const Int_t NPRE, 
		const Int_t NPOST,
		const Float_t gridSizeR, 
		const Float_t ratioZ, 
		const Float_t ratioPhi,
		std::vector<TMatrixD **> & tvArrayV,	
		std::vector<TMatrixD **> & tvCharge,
		std::vector<TMatrixD **> & tvResidu,
		std::vector<float> &coef1,
		std::vector<float> &coef2, 
		std::vector<float> &coef3,
		std::vector<float> &coef4,
		std::vector<float> &icoef4
	);

		
	Double_t fMaxExact;
	
	
	Double_t GetErrorExact(TMatrixD **curArrayV, TMatrixD **tempArrayV,   const Int_t phislices);
	Double_t GetErrorConv(TMatrixD **curArrayV, TMatrixD **prevArrayV, const Int_t phislices);
	
	
	Bool_t fExactPresent;
/// \cond CLASSIMP
	ClassDef(AliTPCPoissonSolver,5);
/// \endcond
};

}
}
#endif	// O2_TPC_ALITPCPOISSONSOLVER_H
