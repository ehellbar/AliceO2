#ifndef O2_TPC_AliTPC3DCylindricalInterpolatorFull_H
#define O2_TPC_AliTPC3DCylindricalInterpolatorFull_H


/// \class AliTPC3DInterpolatorFull
/// \brief Interpolator for cylindrical coordinate with r,phi,z different coordinates
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Jan 5, 2016

#include <TMath.h>
#include <TMatrixD.h>

namespace o2 {
namespace TPC {
    
class AliTPC3DCylindricalInterpolatorFull {
public:

	void SetNR(Int_t nr) {fNR = nr;}
	void SetNPhi(Int_t nphi) {fNPhi = nphi;}
	void SetNZ(Int_t nz) {fNZ = nz;}
	
	Int_t GetNR() {return fNR;}
	Int_t GetNPhi() {return fNPhi;}
	Int_t GetNZ() {return fNZ;}
	
	
	void SetRList(Double_t *rlist) {fRlist = rlist; }
	void SetPhiList(Double_t *philist) {fPhilist = philist; }
	void SetZList(Double_t *zlist) {fZlist = zlist; }

	void SetRListNormalized(Double_t *rlist) {fRlistNormalized = rlist; }
	void SetPhiListNormalized(Double_t *philist) {fPhilistNormalized = philist; }
	void SetZListNormalized(Double_t *zlist) {fZlistNormalized = zlist; }
	
	void SetVals(Double_t *vals) {fVals = vals;}
	
	void SetVals(TMatrixD **mvals, TMatrixD **rlist, TMatrixD **philist, TMatrixD **zlist);
	void SetValsCartesian(TMatrixD **mvals, TMatrixD **rlist, TMatrixD **philist, TMatrixD **zlist);
	void SetVals(TMatrixD **mvals, TMatrixD **rlist, TMatrixD **philist, TMatrixD **zlist, Int_t jy);

	
	AliTPC3DCylindricalInterpolatorFull(Int_t nr, Int_t nz, Int_t nphi, Int_t stepr, Int_t stepz, Int_t stepphi, Int_t intType);	
	AliTPC3DCylindricalInterpolatorFull();
	virtual ~AliTPC3DCylindricalInterpolatorFull();
  
	Double_t GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ
	);



	Double_t GetValueCartesian
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ
	);
	

	Double_t GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Int_t minzindex
	);



	Double_t GetValueCartesian
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Int_t minzindex
	);
	Double_t GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z
	);
       
	
	void SetOrder(Int_t order) {
		fOrder = order;
	}
  
  
	
	void InitCubicSpline();
	void InitRBFWeight();
	void InitRBFWeightCartesian();

	Int_t GetIrregularGridSize() {return fIrregularGridSize;}
	void SetIrregularGridSize(Int_t size) {fIrregularGridSize = size;}

	void SetKernelType(Int_t kernelType)
	{
		fKernelType = kernelType;
	}

	Int_t GetKernelType() {return fKernelType;}

	///< Enumeration of Poisson Solver Strategy Type
	enum RBFKernelType {
		kRBFMultiQuadratic = 0, 
		kRBFInverseMultiQuadratic = 1,  
		kRBFThinPlateSpline = 2,
		kRBFGaussian = 3	       
  	};

	
private:

 	Int_t fOrder;			///< Order of interpolation, 1 - linear, 2 - quadratic, 3 - cubic
	Int_t fType; // 0 INVERSE WEIGHT, 1 RBF FULL, 2 RBF Half
	Int_t fKernelType; /// < type kernel RBF 1--5

	Int_t fIrregularGridSize; // size when interpolating for irregular grid

	Int_t fNR;				///< Grid size in direction of R
	Int_t fNPhi;			///< Grid size in direction of Phi
	Int_t fNZ;				///< Grid size in direction of Z
	
	Int_t fMinZIndex; ///<index z minimal as lower bound


	Int_t fStepR;
	Int_t fStepZ;
	Int_t fStepPhi; 	

	Double_t *fVals;  ///< 3D for storing known values interpolation should be in size fNR*fNPhi*fNZ	
    Double_t *fValsNormalized;

	Double_t fRadiusRBF0;

	Int_t * fRBFWeightLookUp;


	Double_t *fRlist; ///< coordinate in R (cm) (should be increasing) in 3D
	Double_t *fPhilist; ///< coordinate in philist (rad) (should be increasing) 0 <= < 2 pi (cyclic) in 3D
	Double_t *fZlist; ///< coordinate in z list (cm) (should be increasing) in 3D
	Double_t *fRlistNormalized; ///< coordinate in R (cm) (should be increasing) in 3D
	Double_t *fPhilistNormalized; ///< coordinate in philist (rad) (should be increasing) 0 <= < 2 pi (cyclic) in 3D
	Double_t *fZlistNormalized; ///< coordinate in z list (cm) (should be increasing) in 3D
	
	//Double_t *fSecondDerR; ///< store second derivative of cubic interpolation in r direction
	Double_t *fSecondDerZ; ///< store second derivative of cubic interpolation in z direction
	//Double_t *fSecondDerPhi; ///< store second derivative of cubic interpolation in phi direction

	//Double_t *fWeigthRBF; ///< weight for RBF
	Double_t *fRBFWeight; 
	
	Bool_t fIsAllocatingLookUp; ///< is allocating memory?
	Bool_t fIsInitCubic;
	
	Double_t InterpolatePhi
	( 
		Double_t xArray[], 
		const Int_t ilow,
		const Int_t nx,
		Double_t yArray[],
		Double_t x 
	);

	Double_t InterpolatePhi
	( 
		Double_t xArray[], 
		Int_t offsetX,
		const Int_t ilow,
		const Int_t nx,
		Double_t yArray[],
		Double_t x 
	);
	
	
	Double_t Interpolate3DTableCyl
	( 
		Double_t r,   
		Double_t z,
		Double_t phi,
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepZ,
		Int_t stepPhi
	);

	Double_t Interpolate3DTableCylIDW
	( 
		Double_t r,   
		Double_t z,
		Double_t phi,
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepZ,
		Int_t stepPhi
	);


	Double_t Interpolate3DTableCylRBF
	( 
		Double_t r,   
		Double_t z,
		Double_t phi,
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepZ,
		Int_t stepPhi,
		Double_t radiusRBF0
	);


	Double_t Interpolate3DTableCylRBFCartesian
	( 
		Double_t r,   
		Double_t z,
		Double_t phi,
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepZ,
		Int_t stepPhi,
		Double_t radiusRBF0
	);

	Double_t Interpolate3DTableCyl
	( 
		Double_t r,   
		Double_t z,
		Double_t phi
	);
	
	void Search
	( 
		Int_t n, 
		const Double_t xArray[], 
		Double_t x, 
		Int_t &low 
	);

	void Search
	( 
		Int_t n, 
		Double_t  * xArray, 
		Int_t offset,
		Double_t x, 
		Int_t &low 
	);
	
	
	
	Double_t Interpolate
	( 
		Double_t xArray[], 
		Double_t yArray[],
		Double_t x
	);


	Double_t Interpolate
	( 
		Double_t xArray[], 
		Int_t offsetX,
		Double_t yArray[],
		Int_t offsetY,
		Double_t x
	);

	
	void InitCubicSpline
	(	
		Double_t * xArray,
		Double_t * yArray,
		const Int_t n,
		Double_t * y2Array,
		const Int_t skip
	);
	
	
	
	//void InitCubicSplinePhi
	//(	
	//	Double_t * xArray,
	//	Double_t * yArray,
	//	const Int_t n,
	//	Double_t * y2Array,
	//	const Int_t skip
	//);
	
	
	
	Double_t InterpolateCubicSpline
	(	
	
		Double_t * xArray, 
		Double_t * yArray,
		Double_t * y2Array,	
		const Int_t nxArray,
		const Int_t nyArray,
		const Int_t ny2Array,
		Double_t x,
		const Int_t skip
	);
	
	
	Double_t InterpolateCubicSplinePhi
	(	
		Double_t * xArray, 
		Double_t * yArray,
		Double_t * y2Array,	
		const Int_t n,
		Double_t x,
		const Int_t skip
	);

	Double_t Distance
	(
	   Double_t r0,
	   Double_t phi0,
	   Double_t z0,
	   Double_t r,
	   Double_t phi,
	   Double_t z
	);
// RBF

	void RBFWeight 
	( 
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0, 
		Int_t kernelType,
		Double_t * weight
	);



	void RBFWeightCartesian
	( 
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0, 
		Int_t kernelType,
		Double_t * weight
	);

	void GetRBFWeight 
	( 
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0, 
		Int_t kernelType,
		Double_t * weight
	);


	void GetRBFWeightCartesian
	( 
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0, 
		Int_t kernelType,
		Double_t * weight
	);

	void Phi
	( 
		Int_t n, 
		Double_t r[], 
		Double_t r0, 
		Double_t v[] 
	);

// RBF functions
	void rbf1 
	( 
		Int_t n, 
		Double_t r[], 
		Double_t r0, 
		Double_t v[] 
	);

	void rbf2 
	( 
		Int_t n, 
		Double_t r[], 
		Double_t r0, 
		Double_t v[] 
	);

	void rbf3 
	( 
		Int_t n, 
		Double_t r[], 
		Double_t r0, 
		Double_t v[] 
	);

	void rbf4
	( 
		Int_t n, 
		Double_t r[], 
		Double_t r0, 
		Double_t v[] 
	);


	Double_t InterpRBF
	( 
		Double_t r,
		Double_t phi,
		Double_t z,
		Int_t startR,
		Int_t startPhi,
		Int_t startZ,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0,
		Int_t kernelType,
		Double_t * weight
	);



	Double_t InterpRBFCartesian
	( 
		Double_t r,
		Double_t phi,
		Double_t z,
		Int_t startR,
		Int_t startPhi,
		Int_t startZ,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0,
		Int_t kernelType,
		Double_t * weight
	);

	void InitRBFWeight
	(
		Int_t jy
	);


	// this thin spline plate use half information of cubes
	void RBFWeightHalf 
	( 
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0, 
		Int_t kernelType,
		Double_t * weight
	);


	Double_t InterpRBFHalf
	( 
		Double_t r,
		Double_t phi,
		Double_t z,
		Int_t startR,
		Int_t startPhi,
		Int_t startZ,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0,
		Int_t kernelType,
		Double_t * weight
	);



	void GetRBFWeightHalf
	( 
		Int_t rindex,
		Int_t zindex,
		Int_t phiindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Double_t radius0, 
		Int_t kernelType,
		Double_t * weight
	);


	Double_t GetRadius0RBF
	(
		const Int_t indexr,
		const Int_t indexphi,
		const Int_t indexz
	);
	//void SetRadiusRBF0();

/// \cond CLASSIMP
	ClassDef(AliTPC3DCylindricalInterpolatorFull,1);
/// \endcond
};

}
}

#endif
