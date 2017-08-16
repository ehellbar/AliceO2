#ifndef O2_TPC_AliTPC3DCylindricalInterpolator_H
#define O2_TPC_AliTPC3DCylindricalInterpolator_H


/// \class AliTPC3DInterpolator
/// \brief Interpolator for cylindrical coordinate
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Jan 5, 2016

#include <TMath.h>
#include <TMatrixD.h>

namespace o2 {
namespace TPC {

class AliTPC3DCylindricalInterpolator {
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
	
  void SetVals(Double_t *vals) {fVals = vals;}
	
  void SetVals(TMatrixD **mvals);
	
  AliTPC3DCylindricalInterpolator();
  virtual ~AliTPC3DCylindricalInterpolator();
  
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
	
 private:

  Int_t fOrder;			///< Order of interpolation, 1 - linear, 2 - quadratic, 3 - cubic
  Int_t fNR;				///< Grid size in direction of R
  Int_t fNPhi;			///< Grid size in direction of Phi
  Int_t fNZ;				///< Grid size in direction of Z
	
	
  Double_t *fVals;  ///< 3D for storing known values interpolation should be in size fNR*fNPhi*fNZ	
  Double_t *fRlist; ///< coordinate in R (cm) (should be increasing)
  Double_t *fPhilist; ///< coordinate in philist (rad) (should be increasing) 0 <= < 2 pi (cyclic)
  Double_t *fZlist; ///< coordinate in z list (cm) (should be increasing)
	
  //Double_t *fSecondDerR; ///< store second derivative of cubic interpolation in r direction
  Double_t *fSecondDerZ; ///< store second derivative of cubic interpolation in z direction
  //Double_t *fSecondDerPhi; ///< store second derivative of cubic interpolation in phi direction
	
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
	
	
	
  Double_t Interpolate
    ( 
     Double_t xArray[], 
     Double_t yArray[],
     Double_t x
      );
	

  // natural cubic spline
  void InitCubicSpline
    (	
     Double_t * xArray,
     Double_t * yArray,
     const Int_t n,
     Double_t * y2Array,
     const Int_t skip
	);
	
  // clamped cubic spline
  void InitCubicSpline
    (	
     Double_t * xArray,
     Double_t * yArray,
     const Int_t n,
     Double_t * y2Array,
     const Int_t skip,
     Double_t yp0,
     Double_t ypn1
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
  /// \cond CLASSIMP
  ClassDef(AliTPC3DCylindricalInterpolator,1);
  /// \endcond
};
 
}
}

#endif	// O2_TPC_AliTPC3DCylindricalInterpolator_H
